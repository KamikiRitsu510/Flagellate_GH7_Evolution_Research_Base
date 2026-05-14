"""Command-line interface for phylofetch."""

import argparse
import sys
import logging

from .utils import set_email
from .cds import batch_fetch_cds
from .taxonomy import batch_fetch_taxonomy

logger = logging.getLogger(__name__)


def _read_accessions(path: str) -> list:
    """Read one accession per line from a plain text file."""
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phylofetch",
        description=(
            "phylofetch — Batch fetch CDS sequences and taxonomy lineages from NCBI.\n"
            "\n"
            "Two subcommands are available:\n"
            "  cds        Download coding DNA sequences (FASTA output)\n"
            "  taxonomy   Download full taxonomic lineages (Excel / TSV output)\n"
            "\n"
            "NCBI requires a valid email with every request (--email / -e).\n"
            "Rate limiting (≤3 req/s) and retry logic are applied automatically."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # Download CDS for accessions in a file
  phylofetch cds --input accs.txt --output sequences.fasta --email you@example.com

  # Download CDS for accessions passed directly
  phylofetch cds BAB64565.3 AAY83390.3 --email you@example.com

  # Download taxonomy → Excel + TSV
  phylofetch taxonomy --input accs.txt \\
      --excel taxonomy.xlsx --tsv taxonomy.tsv --email you@example.com

  # Include extra taxonomic ranks (e.g. subphylum)
  phylofetch taxonomy --input accs.txt --extra-ranks subphylum tribe \\
      --tsv full_taxonomy.tsv --email you@example.com

  # Suppress progress output
  phylofetch cds --input accs.txt -o out.fasta -e you@example.com --quiet

Input file format (accs.txt)
-----------------------------
  One accession per line. Lines starting with # are ignored.
  Example:
    # Flagellate GH7 sequences
    BAB64565.3
    AAY83390.3
    WP_165634028.1

Output columns (taxonomy)
--------------------------
  accession | taxid | lineage_string | superkingdom | kingdom |
  phylum | class | order | family | genus | species [+ extras]

See also: EXAMPLE.md for a full worked walkthrough.
""",
    )
    parser.add_argument(
        "--email", "-e",
        required=True,
        help="Valid email address for NCBI Entrez (required by NCBI policy).",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        default=True,
        help="Print per-accession progress (default: on).",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress per-accession progress messages.",
    )
    parser.add_argument(
        "--retries", "-r",
        type=int,
        default=3,
        metavar="N",
        help="Retry attempts on network failure (default: 3).",
    )

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    # ── cds sub-command ──────────────────────────────────────────────────────
    cds_p = sub.add_parser(
        "cds",
        help="Fetch CDS nucleotide sequences.",
        description="Fetch CDS (coding DNA sequences) for a list of NCBI accessions.",
    )
    cds_p.add_argument(
        "accessions",
        nargs="*",
        metavar="ACCESSION",
        help="Accession(s) to fetch (positional).",
    )
    cds_p.add_argument(
        "--input", "-i",
        metavar="FILE",
        help="Text file with one accession per line.",
    )
    cds_p.add_argument(
        "--output", "-o",
        metavar="FILE",
        help="Output FASTA file (default: print to stdout).",
    )

    # ── taxonomy sub-command ─────────────────────────────────────────────────
    tax_p = sub.add_parser(
        "taxonomy",
        help="Fetch taxonomy lineages.",
        description="Fetch taxonomy lineages for a list of NCBI accessions.",
    )
    tax_p.add_argument(
        "accessions",
        nargs="*",
        metavar="ACCESSION",
        help="Accession(s) to fetch (positional).",
    )
    tax_p.add_argument(
        "--input", "-i",
        metavar="FILE",
        help="Text file with one accession per line.",
    )
    tax_p.add_argument(
        "--excel", "-x",
        metavar="FILE",
        help="Save results to an Excel (.xlsx) file.",
    )
    tax_p.add_argument(
        "--tsv", "-t",
        metavar="FILE",
        help="Save results to a tab-separated (.tsv) file.",
    )
    tax_p.add_argument(
        "--extra-ranks",
        nargs="+",
        metavar="RANK",
        default=[],
        help="Extra NCBI rank names to include (e.g. subphylum tribe).",
    )

    return parser


def main(argv=None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    # Collect accessions
    accessions = list(args.accessions)
    if hasattr(args, "input") and args.input:
        accessions += _read_accessions(args.input)

    if not accessions:
        parser.error("No accessions provided. Use positional args or --input FILE.")

    set_email(args.email)
    verbose = args.verbose and not args.quiet

    # ── cds ──────────────────────────────────────────────────────────────────
    if args.command == "cds":
        results = batch_fetch_cds(
            accessions,
            output_fasta=args.output,
            retries=args.retries,
            verbose=verbose,
        )
        if not args.output:
            # Print FASTA to stdout
            for acc, seq in results.items():
                print(f">{acc}")
                # Wrap at 70 chars
                for i in range(0, len(seq), 70):
                    print(seq[i:i + 70])
        found = len(results)
        total = len(accessions)
        if verbose:
            print(f"\nDone: {found}/{total} sequences retrieved.")
        return 0 if found > 0 else 1

    # ── taxonomy ─────────────────────────────────────────────────────────────
    if args.command == "taxonomy":
        if not args.excel and not args.tsv:
            # Default: print TSV to stdout
            import io
            buf = io.StringIO()
            df = batch_fetch_taxonomy(
                accessions,
                extra_ranks=args.extra_ranks,
                retries=args.retries,
                verbose=verbose,
            )
            print(df.to_csv(sep="\t", index=False), end="")
        else:
            df = batch_fetch_taxonomy(
                accessions,
                output_excel=args.excel,
                output_tsv=args.tsv,
                extra_ranks=args.extra_ranks,
                retries=args.retries,
                verbose=verbose,
            )
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
