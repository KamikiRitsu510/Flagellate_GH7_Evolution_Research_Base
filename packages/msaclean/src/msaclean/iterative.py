"""Iterative alignment cleaning with external realignment.

Implements the iterative curation loop originally prototyped in
``iterative_clean.sh``: realign → trim → remove worst sequences → repeat.
Requires MAFFT (and optionally trimAl) to be installed and on PATH.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Collection, List, Optional, Sequence, Union

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

from .core import _accession, clean_alignment, gappiest_sequences

logger = logging.getLogger(__name__)


# ── Data types ────────────────────────────────────────────────────────────────

@dataclass
class IterationResult:
    """Outcome of a single cleaning round."""
    round: int
    seq_count: int
    alignment_length: int
    removed: list[str] = field(default_factory=list)
    stopped_reason: str = ""


@dataclass
class IterativeCleanResult:
    """Full result of :func:`iterative_clean`."""
    rounds: list[IterationResult] = field(default_factory=list)
    final_alignment: Optional[MultipleSeqAlignment] = None
    output_file: Optional[Path] = None

    @property
    def log_table(self) -> str:
        """Human-readable table of round-by-round statistics."""
        lines = [f"{'Round':>6}  {'Seqs':>6}  {'Aln len':>8}  Removed"]
        for r in self.rounds:
            removed_str = ", ".join(r.removed) if r.removed else "—"
            lines.append(f"{r.round:>6}  {r.seq_count:>6}  {r.alignment_length:>8}  {removed_str}")
        if self.rounds and self.rounds[-1].stopped_reason:
            lines.append(f"\nStopped: {self.rounds[-1].stopped_reason}")
        return "\n".join(lines)


# ── External tool runners ─────────────────────────────────────────────────────

def _run_mafft(
    input_fasta: Path,
    output_fasta: Path,
    threads: int = 4,
    mode: str = "localpair",
    maxiterate: int = 1000,
) -> None:
    """Run MAFFT to realign sequences."""
    cmd = [
        "mafft",
        f"--{mode}",
        "--maxiterate", str(maxiterate),
        "--thread", str(threads),
        "--quiet",
        str(input_fasta),
    ]
    with open(output_fasta, "w") as out_fh:
        subprocess.run(cmd, stdout=out_fh, check=True)


def _run_trimal(
    input_fasta: Path,
    output_fasta: Path,
    method: str = "automated1",
) -> None:
    """Run trimAl to trim alignment columns."""
    cmd = [
        "trimal",
        "-in", str(input_fasta),
        "-out", str(output_fasta),
        f"-{method}",
    ]
    subprocess.run(cmd, check=True, capture_output=True)


# ── Core iterative loop ───────────────────────────────────────────────────────

def iterative_clean(
    input_fasta: Union[str, Path],
    output_fasta: Optional[Union[str, Path]] = None,
    whitelist: Collection[str] = (),
    remove_per_round: int = 2,
    max_rounds: int = 20,
    stop_at_n: Optional[int] = None,
    max_gap_fraction: float = 1.0,
    mafft_threads: int = 4,
    mafft_mode: str = "localpair",
    mafft_maxiterate: int = 1000,
    use_trimal: bool = True,
    trimal_method: str = "automated1",
    verbose: bool = True,
) -> IterativeCleanResult:
    """Iteratively realign and remove the gappiest sequences.

    Each round:
    1. Realign remaining sequences with MAFFT.
    2. (Optional) Trim alignment columns with trimAl.
    3. Identify the *remove_per_round* gappiest non-whitelisted sequences.
    4. Remove them and proceed to the next round.

    Stops when the target sequence count is reached, no removable sequences
    remain, or *max_rounds* is exceeded.

    Parameters
    ----------
    input_fasta : str or Path
        Input FASTA file (unaligned or pre-aligned sequences).
    output_fasta : str or Path, optional
        Where to write the final cleaned alignment.
    whitelist : collection of str, optional
        Accessions that are never removed.
    remove_per_round : int
        How many sequences to remove per round (default 2).
    max_rounds : int
        Maximum number of rounds (default 20).
    stop_at_n : int, optional
        Stop when this many sequences remain. If None, runs until no
        more non-whitelisted sequences can be removed.
    max_gap_fraction : float
        Only remove sequences whose gap fraction exceeds this threshold
        (default 1.0, i.e. always remove the top *remove_per_round* gappiest).
    mafft_threads : int
        Threads for MAFFT (default 4).
    mafft_mode : str
        MAFFT alignment mode: ``"localpair"`` (L-INS-i), ``"globalpair"``,
        ``"genafpair"``, ``"auto"`` (default: ``"localpair"``).
    mafft_maxiterate : int
        MAFFT ``--maxiterate`` value (default 1000).
    use_trimal : bool
        Run trimAl after each MAFFT alignment (default True).
    trimal_method : str
        trimAl trimming method (default ``"automated1"``).
    verbose : bool
        Print round-by-round progress (default True).

    Returns
    -------
    IterativeCleanResult
        Contains per-round statistics and the final alignment.

    Raises
    ------
    FileNotFoundError
        If MAFFT (or trimAl when ``use_trimal=True``) is not on PATH.

    Examples
    --------
    >>> from msaclean import iterative_clean
    >>> WHITELIST = {
    ...     "BAB64565.3", "AAY83390.3", "BAB64553.1",
    ... }
    >>> result = iterative_clean(
    ...     "GH7_clustered.fasta",
    ...     output_fasta="GH7_curated.fasta",
    ...     whitelist=WHITELIST,
    ...     remove_per_round=2,
    ...     stop_at_n=20,
    ... )
    >>> print(result.log_table)
    """
    _check_executables(use_trimal)

    input_fasta = Path(input_fasta)
    result = IterativeCleanResult()
    wl_set = set(whitelist)

    # Read initial sequences (unaligned)
    sequences = list(SeqIO.parse(str(input_fasta), "fasta"))
    if verbose:
        print(f"Starting with {len(sequences)} sequences "
              f"({len(wl_set)} whitelisted).")

    with tempfile.TemporaryDirectory(prefix="msaclean_") as tmpdir:
        tmpdir = Path(tmpdir)
        current_fasta = tmpdir / "current.fasta"

        # Write initial unaligned FASTA
        SeqIO.write(sequences, str(current_fasta), "fasta")

        for rnd in range(max_rounds):
            # ── Align ──
            aligned_fasta = tmpdir / f"round{rnd}_aligned.fasta"
            _run_mafft(current_fasta, aligned_fasta,
                       threads=mafft_threads, mode=mafft_mode,
                       maxiterate=mafft_maxiterate)

            # ── Trim (optional) ──
            if use_trimal:
                trimmed_fasta = tmpdir / f"round{rnd}_trimmed.fasta"
                try:
                    _run_trimal(aligned_fasta, trimmed_fasta, method=trimal_method)
                    eval_fasta = trimmed_fasta
                except subprocess.CalledProcessError:
                    logger.warning("trimAl failed on round %d; skipping trim.", rnd)
                    eval_fasta = aligned_fasta
            else:
                eval_fasta = aligned_fasta

            # ── Read alignment ──
            aln = AlignIO.read(str(eval_fasta), "fasta")
            aln_len = aln.get_alignment_length()
            seq_count = len(aln)

            if verbose:
                print(f"Round {rnd:>2}: {seq_count:>4} seqs, "
                      f"{aln_len:>6} columns", end="")

            # ── Stop condition ──
            if stop_at_n is not None and seq_count <= stop_at_n:
                ir = IterationResult(rnd, seq_count, aln_len,
                                     stopped_reason=f"reached stop_at_n={stop_at_n}")
                result.rounds.append(ir)
                if verbose:
                    print(f"  → STOP (target {stop_at_n} reached)")
                # Save the aligned/trimmed result
                _finalize(aln, output_fasta, result)
                return result

            # ── Identify worst sequences ──
            worst = gappiest_sequences(aln, n=remove_per_round, whitelist=wl_set)
            # Optionally filter by gap threshold
            worst = [(acc, gf) for acc, gf in worst if gf > max_gap_fraction
                     ] if max_gap_fraction < 1.0 else worst

            if not worst:
                ir = IterationResult(rnd, seq_count, aln_len,
                                     stopped_reason="no removable sequences remaining")
                result.rounds.append(ir)
                if verbose:
                    print("  → STOP (nothing to remove)")
                _finalize(aln, output_fasta, result)
                return result

            removed_accs = [acc for acc, _ in worst]
            ir = IterationResult(rnd, seq_count, aln_len, removed=removed_accs)
            result.rounds.append(ir)
            if verbose:
                print(f"  → removing: {', '.join(removed_accs)}")

            # ── Remove and write next round input (unaligned) ──
            remove_set = set(removed_accs)
            sequences = [r for r in sequences if _accession(r) not in remove_set]
            SeqIO.write(sequences, str(current_fasta), "fasta")

        # Ran out of rounds — do a final alignment
        if verbose:
            print(f"Reached max_rounds={max_rounds}; performing final alignment.")
        aligned_fasta = tmpdir / "final_aligned.fasta"
        _run_mafft(current_fasta, aligned_fasta,
                   threads=mafft_threads, mode=mafft_mode,
                   maxiterate=mafft_maxiterate)
        aln = AlignIO.read(str(aligned_fasta), "fasta")
        result.rounds[-1].stopped_reason = f"max_rounds={max_rounds} reached"
        _finalize(aln, output_fasta, result)

    return result


def _finalize(
    aln: MultipleSeqAlignment,
    output_fasta: Optional[Union[str, Path]],
    result: IterativeCleanResult,
) -> None:
    result.final_alignment = aln
    if output_fasta:
        AlignIO.write(aln, str(output_fasta), "fasta")
        result.output_file = Path(output_fasta)


def _check_executables(use_trimal: bool) -> None:
    if not shutil.which("mafft"):
        raise FileNotFoundError(
            "MAFFT is not installed or not on PATH. "
            "Install it with: conda install -c bioconda mafft"
        )
    if use_trimal and not shutil.which("trimal"):
        raise FileNotFoundError(
            "trimAl is not installed or not on PATH. "
            "Install it with: conda install -c bioconda trimal  "
            "— or pass use_trimal=False to skip column trimming."
        )
