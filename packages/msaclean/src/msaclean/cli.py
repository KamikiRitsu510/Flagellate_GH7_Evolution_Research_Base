"""Command-line interface for msaclean."""

import argparse
import logging
import sys
from pathlib import Path

from .core import clean_alignment_file
from .iterative import iterative_clean

logger = logging.getLogger(__name__)


def _read_whitelist(path: str) -> list[str]:
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="msaclean",
        description=(
            "msaclean — Clean and curate multiple sequence alignments for phylogenetics.\n"
            "\n"
            "Three subcommands:\n"
            "  clean       One-pass gap filtering on an existing alignment\n"
            "  iterative   Iterative MAFFT realignment + gap-based removal (requires mafft)\n"
            "  stats       Inspect gap statistics without modifying the alignment\n"
            "\n"
            "Whitelist (-w) protects specific accessions from removal regardless of gap content.\n"
            "Accession is extracted as the first whitespace-delimited token in each FASTA header."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # Inspect gap distribution first
  msaclean stats GH7_aligned.fasta --top 20 -w whitelist.txt

  # Remove sequences with >90% gaps, keep whitelist
  msaclean clean GH7_aligned.fasta -o GH7_clean.fasta \\
      -w whitelist.txt --max-gap 0.9

  # Iterative realignment: remove 2 worst per round, stop at 50 sequences
  msaclean iterative GH7_clustered.fasta -o GH7_curated.fasta \\
      -w whitelist.txt --stop-at 50 --rounds 2 --mafft-threads 8

  # Faster MAFFT mode, skip trimAl
  msaclean iterative GH7_clustered.fasta -o GH7_curated.fasta \\
      --mafft-mode auto --no-trimal --stop-at 30 -w whitelist.txt

Whitelist file format
---------------------
  One accession per line. Lines starting with # are ignored.
  Example:
    # Core flagellate sequences — never remove
    BAB64565.3
    AAY83390.3

MAFFT modes
-----------
  localpair   L-INS-i, highest accuracy, slow   (default)
  globalpair  G-INS-i, global pairwise
  genafpair   E-INS-i, for sequences with long unalignable regions
  auto        Automatic selection, fastest

See also: EXAMPLE.md for a full worked walkthrough.
""",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Verbose output."
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="Suppress progress messages."
    )

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    # ── clean ────────────────────────────────────────────────────────────────
    clean_p = sub.add_parser(
        "clean",
        help="Remove empty/gappy sequences from an existing alignment.",
    )
    clean_p.add_argument("input", help="Input alignment file (FASTA).")
    clean_p.add_argument("--output", "-o", required=True, help="Output FASTA file.")
    clean_p.add_argument(
        "--whitelist", "-w", metavar="FILE",
        help="Text file with accessions to always keep (one per line).",
    )
    clean_p.add_argument(
        "--max-gap", type=float, default=1.0, metavar="FRAC",
        help="Remove sequences with gap fraction above this (0–1, default 1.0).",
    )
    clean_p.add_argument(
        "--min-residues", type=int, default=1, metavar="N",
        help="Minimum non-gap residues to keep a sequence (default 1).",
    )
    clean_p.add_argument(
        "--format", "-f", default="fasta", metavar="FMT",
        help="Biopython format string (default: fasta).",
    )

    # ── iterative ────────────────────────────────────────────────────────────
    iter_p = sub.add_parser(
        "iterative",
        help="Iteratively realign (MAFFT) and remove gappiest sequences.",
    )
    iter_p.add_argument("input", help="Input FASTA (unaligned or pre-aligned).")
    iter_p.add_argument("--output", "-o", required=True, help="Output alignment FASTA.")
    iter_p.add_argument(
        "--whitelist", "-w", metavar="FILE",
        help="Text file with accessions to never remove (one per line).",
    )
    iter_p.add_argument(
        "--stop-at", type=int, default=None, metavar="N",
        help="Stop when this many sequences remain.",
    )
    iter_p.add_argument(
        "--rounds", type=int, default=2, metavar="N",
        help="Sequences to remove per round (default 2).",
    )
    iter_p.add_argument(
        "--max-rounds", type=int, default=20, metavar="N",
        help="Maximum number of rounds (default 20).",
    )
    iter_p.add_argument(
        "--max-gap", type=float, default=1.0, metavar="FRAC",
        help="Only remove sequences above this gap fraction (default 1.0).",
    )
    iter_p.add_argument(
        "--mafft-threads", type=int, default=4, metavar="N",
        help="Threads for MAFFT (default 4).",
    )
    iter_p.add_argument(
        "--mafft-mode", default="localpair",
        choices=["localpair", "globalpair", "genafpair", "auto"],
        help="MAFFT alignment mode (default: localpair = L-INS-i).",
    )
    iter_p.add_argument(
        "--no-trimal", action="store_true",
        help="Skip trimAl column-trimming step.",
    )

    # ── stats ─────────────────────────────────────────────────────────────────
    stats_p = sub.add_parser(
        "stats",
        help="Report gap statistics for each sequence in an alignment.",
    )
    stats_p.add_argument("input", help="Input alignment FASTA.")
    stats_p.add_argument(
        "--top", type=int, default=20, metavar="N",
        help="Show the N gappiest sequences (default 20).",
    )
    stats_p.add_argument(
        "--whitelist", "-w", metavar="FILE",
        help="Text file of accessions to flag as protected.",
    )

    return parser


def _cmd_clean(args) -> int:
    whitelist = _read_whitelist(args.whitelist) if args.whitelist else []
    verbose = args.verbose and not args.quiet

    report = clean_alignment_file(
        args.input,
        args.output,
        whitelist=whitelist,
        max_gap_fraction=args.max_gap,
        min_residues=args.min_residues,
        fmt=args.format,
    )

    if verbose or not args.quiet:
        print(f"\n=== msaclean report ===")
        print(f"Input sequences : {report['original']}")
        print(f"Protected       : {report['protected']}")
        print(f"Removed (empty) : {report['removed_empty']}")
        print(f"Removed (gappy) : {report['removed_gappy']}")
        print(f"Kept            : {report['kept']}")
        print(f"Output          : {args.output}")
        if report["removed_empty_ids"]:
            print(f"\nEmpty sequences : {', '.join(report['removed_empty_ids'])}")
        if report["removed_gappy_ids"]:
            print(f"Gappy sequences : {', '.join(report['removed_gappy_ids'])}")
    return 0


def _cmd_iterative(args) -> int:
    whitelist = _read_whitelist(args.whitelist) if args.whitelist else []
    verbose = not args.quiet

    result = iterative_clean(
        args.input,
        output_fasta=args.output,
        whitelist=whitelist,
        remove_per_round=args.rounds,
        max_rounds=args.max_rounds,
        stop_at_n=args.stop_at,
        max_gap_fraction=args.max_gap,
        mafft_threads=args.mafft_threads,
        mafft_mode=args.mafft_mode,
        use_trimal=not args.no_trimal,
        verbose=verbose,
    )

    if verbose:
        print(f"\n{result.log_table}")
        if result.final_alignment:
            print(f"\nFinal: {len(result.final_alignment)} sequences → {args.output}")
    return 0


def _cmd_stats(args) -> int:
    from Bio import AlignIO
    from .core import gappiest_sequences, _accession, _gap_fraction

    aln = AlignIO.read(args.input, "fasta")
    whitelist = set(_read_whitelist(args.whitelist)) if args.whitelist else set()

    total = len(aln)
    aln_len = aln.get_alignment_length()
    print(f"Alignment: {total} sequences × {aln_len} columns")
    print(f"\n{'Accession':<30} {'Gap%':>7}  {'Residues':>9}  {'Protected':>10}")
    print("-" * 62)

    rows = []
    for rec in aln:
        acc = _accession(rec)
        seq = str(rec.seq)
        gaps = seq.count("-") + seq.count(".")
        res = len(seq) - gaps
        gf = gaps / len(seq) if seq else 1.0
        rows.append((acc, gf, res, acc in whitelist))

    rows.sort(key=lambda x: x[1], reverse=True)
    shown = 0
    for acc, gf, res, prot in rows:
        if shown >= args.top and not prot:
            continue
        prot_str = "✓" if prot else ""
        print(f"{acc:<30} {gf*100:>6.1f}%  {res:>9}  {prot_str:>10}")
        shown += 1

    return 0


def main(argv=None) -> int:
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "clean":
        return _cmd_clean(args)
    if args.command == "iterative":
        return _cmd_iterative(args)
    if args.command == "stats":
        return _cmd_stats(args)
    return 1


if __name__ == "__main__":
    sys.exit(main())
