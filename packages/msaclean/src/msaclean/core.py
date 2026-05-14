"""Core alignment-cleaning functions for msaclean.

Provides gap-based filtering, whitelist protection, and column trimming
for multiple sequence alignments in any Biopython-supported format.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Collection, Iterable, Optional, Sequence, Union

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _accession(record: SeqRecord) -> str:
    """Extract the bare accession from a SeqRecord id (first whitespace-delimited token)."""
    return record.id.strip().split()[0]


def _gap_fraction(record: SeqRecord) -> float:
    """Fraction of gap/dot characters in the aligned sequence."""
    seq = str(record.seq)
    if not seq:
        return 1.0
    gaps = seq.count("-") + seq.count(".")
    return gaps / len(seq)


def _residue_count(record: SeqRecord) -> int:
    """Number of non-gap residues."""
    return len(str(record.seq).replace("-", "").replace(".", ""))


# ── Public API ────────────────────────────────────────────────────────────────

def filter_empty_seqs(
    alignment: MultipleSeqAlignment,
    whitelist: Collection[str] = (),
    min_residues: int = 1,
) -> tuple[MultipleSeqAlignment, list[str]]:
    """Remove sequences that contain fewer than *min_residues* non-gap characters.

    Sequences whose accession appears in *whitelist* are always kept regardless
    of their gap content.

    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Input alignment (Biopython object).
    whitelist : collection of str, optional
        Accession strings that are never removed.
    min_residues : int
        Minimum number of non-gap residues a sequence must have to be kept.
        Default 1 (removes all-gap sequences).

    Returns
    -------
    clean_alignment : MultipleSeqAlignment
        Alignment with empty sequences removed.
    removed : list of str
        Accessions that were removed.
    """
    kept: list[SeqRecord] = []
    removed: list[str] = []

    for rec in alignment:
        acc = _accession(rec)
        if acc in whitelist:
            kept.append(rec)
            continue
        if _residue_count(rec) >= min_residues:
            kept.append(rec)
        else:
            removed.append(acc)
            logger.debug("Removed empty sequence: %s", acc)

    return MultipleSeqAlignment(kept), removed


def filter_gappy_seqs(
    alignment: MultipleSeqAlignment,
    whitelist: Collection[str] = (),
    max_gap_fraction: float = 0.9,
) -> tuple[MultipleSeqAlignment, list[str]]:
    """Remove sequences whose gap content exceeds *max_gap_fraction*.

    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Input alignment.
    whitelist : collection of str, optional
        Accessions that are never removed, regardless of gap content.
    max_gap_fraction : float
        Sequences with a gap fraction above this threshold are removed.
        Range [0, 1]. Default 0.9.

    Returns
    -------
    clean_alignment : MultipleSeqAlignment
        Filtered alignment.
    removed : list of str
        Accessions that were removed.
    """
    kept: list[SeqRecord] = []
    removed: list[str] = []

    for rec in alignment:
        acc = _accession(rec)
        if acc in whitelist:
            kept.append(rec)
            continue
        gf = _gap_fraction(rec)
        if gf <= max_gap_fraction:
            kept.append(rec)
        else:
            removed.append(acc)
            logger.debug("Removed gappy sequence %s (gap=%.1f%%)", acc, gf * 100)

    return MultipleSeqAlignment(kept), removed


def protect_sequences(
    alignment: MultipleSeqAlignment,
    whitelist: Collection[str],
) -> tuple[MultipleSeqAlignment, MultipleSeqAlignment]:
    """Split alignment into protected and unprotected sub-alignments.

    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Input alignment.
    whitelist : collection of str
        Accessions to protect.

    Returns
    -------
    protected : MultipleSeqAlignment
        Sequences whose accession is in *whitelist*.
    unprotected : MultipleSeqAlignment
        All other sequences.
    """
    protected: list[SeqRecord] = []
    unprotected: list[SeqRecord] = []

    wl_set = set(whitelist)
    for rec in alignment:
        if _accession(rec) in wl_set:
            protected.append(rec)
        else:
            unprotected.append(rec)

    return MultipleSeqAlignment(protected), MultipleSeqAlignment(unprotected)


def clean_alignment(
    alignment: MultipleSeqAlignment,
    whitelist: Collection[str] = (),
    max_gap_fraction: float = 1.0,
    min_residues: int = 1,
) -> tuple[MultipleSeqAlignment, dict]:
    """Clean an alignment by removing empty or overly gappy sequences.

    Combines :func:`filter_empty_seqs` and :func:`filter_gappy_seqs` in one
    call. Whitelisted sequences are always retained.

    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Input alignment.
    whitelist : collection of str, optional
        Accessions that are never removed.
    max_gap_fraction : float
        Remove sequences with gap content above this threshold (default 1.0,
        i.e. only all-gap sequences).
    min_residues : int
        Minimum non-gap residues to retain a sequence (default 1).

    Returns
    -------
    clean : MultipleSeqAlignment
        Cleaned alignment.
    report : dict
        Summary with keys ``"original"``, ``"kept"``, ``"removed_empty"``,
        ``"removed_gappy"``, ``"protected"``.

    Examples
    --------
    >>> from Bio import AlignIO
    >>> from msaclean import clean_alignment
    >>> aln = AlignIO.read("GH7_aligned.fasta", "fasta")
    >>> WHITELIST = {"BAB64565.3", "AAY83390.3"}
    >>> clean, report = clean_alignment(aln, whitelist=WHITELIST, max_gap_fraction=0.95)
    >>> print(report)
    """
    original_count = len(alignment)
    wl_set = set(whitelist)
    protected_count = sum(1 for r in alignment if _accession(r) in wl_set)

    # Step 1: remove all-gap / near-empty sequences
    step1, removed_empty = filter_empty_seqs(alignment, whitelist=wl_set,
                                              min_residues=min_residues)
    # Step 2: remove sequences that are too gappy
    step2, removed_gappy = filter_gappy_seqs(step1, whitelist=wl_set,
                                              max_gap_fraction=max_gap_fraction)

    report = {
        "original": original_count,
        "kept": len(step2),
        "removed_empty": len(removed_empty),
        "removed_gappy": len(removed_gappy),
        "removed_empty_ids": removed_empty,
        "removed_gappy_ids": removed_gappy,
        "protected": protected_count,
    }

    return step2, report


# ── File-level convenience wrappers ───────────────────────────────────────────

def clean_alignment_file(
    input_file: Union[str, Path],
    output_file: Union[str, Path],
    whitelist: Collection[str] = (),
    max_gap_fraction: float = 1.0,
    min_residues: int = 1,
    fmt: str = "fasta",
) -> dict:
    """Read an alignment file, clean it, and write the result.

    Parameters
    ----------
    input_file : str or Path
        Path to the input alignment (FASTA, PHYLIP, etc.).
    output_file : str or Path
        Path for the cleaned output alignment.
    whitelist : collection of str, optional
        Accessions that are never removed.
    max_gap_fraction : float
        Sequences above this gap fraction are removed (default 1.0).
    min_residues : int
        Minimum non-gap residues to keep a sequence (default 1).
    fmt : str
        Biopython alignment format string (default ``"fasta"``).

    Returns
    -------
    dict
        Cleaning report (same as :func:`clean_alignment`).

    Examples
    --------
    >>> from msaclean import clean_alignment_file
    >>> report = clean_alignment_file(
    ...     "GH7_aligned.fasta", "GH7_clean.fasta",
    ...     whitelist={"BAB64565.3"}, max_gap_fraction=0.95
    ... )
    >>> print(f"Removed {report['removed_empty']} empty, "
    ...       f"{report['removed_gappy']} gappy sequences")
    """
    aln = AlignIO.read(str(input_file), fmt)
    clean, report = clean_alignment(
        aln,
        whitelist=whitelist,
        max_gap_fraction=max_gap_fraction,
        min_residues=min_residues,
    )
    AlignIO.write(clean, str(output_file), fmt)
    return report


def gappiest_sequences(
    alignment: MultipleSeqAlignment,
    n: int = 5,
    whitelist: Collection[str] = (),
) -> list[tuple[str, float]]:
    """Return the *n* sequences with the highest gap fraction.

    Useful for identifying candidates to remove during iterative curation.

    Parameters
    ----------
    alignment : MultipleSeqAlignment
    n : int
        Number of sequences to return (default 5).
    whitelist : collection of str, optional
        Accessions to exclude from consideration.

    Returns
    -------
    list of (accession, gap_fraction)
        Sorted highest-gap-first.
    """
    wl_set = set(whitelist)
    candidates = [
        (_accession(r), _gap_fraction(r))
        for r in alignment
        if _accession(r) not in wl_set
    ]
    candidates.sort(key=lambda x: x[1], reverse=True)
    return candidates[:n]
