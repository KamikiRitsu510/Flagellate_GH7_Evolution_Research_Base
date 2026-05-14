# msaclean

**Clean and curate multiple sequence alignments for phylogenetic analysis.**

`msaclean` generalises the whitelist-based gap-filtering and iterative MAFFT curation workflow originally developed for GH7 cellulase phylogeny into a reusable pip-installable package. It handles the two most common alignment-quality problems:

- **All-gap / near-empty sequences** — sequences that survived CD-HIT clustering but are nearly entirely gaps after alignment.
- **Iterative contamination** — outlier sequences that distort alignment quality and need to be removed round-by-round with realignment in between.

---

## Installation

```bash
pip install msaclean
```

**Requirements:** Python ≥ 3.8, Biopython ≥ 1.80

For iterative cleaning you also need [MAFFT](https://mafft.cbrc.jp/alignment/software/) and optionally [trimAl](http://trimal.cgenomics.org/) on your PATH:

```bash
conda install -c bioconda mafft trimal
```

---

## Quick start

```python
from msaclean import clean_alignment_file, clean_alignment, iterative_clean
from Bio import AlignIO

# Sequences to always keep, regardless of gap content
WHITELIST = {
    "BAB64565.3", "AAY83390.3", "BAB64553.1", "BAB64554.1",
    "BAB64555.1", "BAB64556.1", "BAB64557.1", "BAB64558.1",
}

# ── One-call file cleaning ─────────────────────────────────────────────────

report = clean_alignment_file(
    "GH7_aligned.fasta",
    "GH7_clean.fasta",
    whitelist=WHITELIST,
    max_gap_fraction=0.95,   # remove sequences with >95% gaps
)
print(f"Removed {report['removed_empty']} empty + "
      f"{report['removed_gappy']} gappy sequences")

# ── Object-level API ───────────────────────────────────────────────────────

aln = AlignIO.read("GH7_aligned.fasta", "fasta")
clean, report = clean_alignment(aln, whitelist=WHITELIST, max_gap_fraction=0.9)

# ── Inspect what's gappy ──────────────────────────────────────────────────

from msaclean import gappiest_sequences
worst = gappiest_sequences(aln, n=10, whitelist=WHITELIST)
for acc, gap_frac in worst:
    print(f"{acc}  {gap_frac:.1%} gaps")

# ── Iterative realignment (requires MAFFT) ─────────────────────────────────

result = iterative_clean(
    "GH7_clustered.fasta",
    output_fasta="GH7_curated.fasta",
    whitelist=WHITELIST,
    remove_per_round=2,
    stop_at_n=20,
)
print(result.log_table)
```

---

## Command-line interface

```bash
# Remove all-gap sequences, protect a whitelist
msaclean clean GH7_aligned.fasta -o GH7_clean.fasta -w whitelist.txt

# Remove sequences with >90% gaps
msaclean clean GH7_aligned.fasta -o GH7_clean.fasta --max-gap 0.9

# Show gap statistics (top 20 gappiest)
msaclean stats GH7_aligned.fasta --top 20 -w whitelist.txt

# Iterative realignment: remove 2 per round, stop at 20 sequences
msaclean iterative GH7_raw.fasta -o GH7_curated.fasta \
    -w whitelist.txt --stop-at 20 --rounds 2

# Use faster MAFFT mode and skip trimAl
msaclean iterative GH7_raw.fasta -o GH7_curated.fasta \
    --mafft-mode auto --no-trimal --stop-at 30
```

`whitelist.txt` is a plain text file with one accession per line (lines starting with `#` are ignored).

---

## API reference

### `clean_alignment_file(input_file, output_file, whitelist=(), max_gap_fraction=1.0, min_residues=1, fmt="fasta") → dict`

Read an alignment file, remove empty/gappy sequences, and write the result. Returns a cleaning report dict.

---

### `clean_alignment(alignment, whitelist=(), max_gap_fraction=1.0, min_residues=1) → (MultipleSeqAlignment, dict)`

Clean a Biopython `MultipleSeqAlignment` object in memory.

Report dict keys: `original`, `kept`, `removed_empty`, `removed_gappy`, `removed_empty_ids`, `removed_gappy_ids`, `protected`.

---

### `filter_empty_seqs(alignment, whitelist=(), min_residues=1) → (alignment, removed_list)`

Remove sequences with fewer than `min_residues` non-gap characters.

---

### `filter_gappy_seqs(alignment, whitelist=(), max_gap_fraction=0.9) → (alignment, removed_list)`

Remove sequences whose gap fraction exceeds `max_gap_fraction`.

---

### `gappiest_sequences(alignment, n=5, whitelist=()) → list[(accession, gap_fraction)]`

Return the `n` sequences with the highest gap fraction, sorted highest-first. Useful for deciding what to remove before running iterative cleaning.

---

### `iterative_clean(input_fasta, output_fasta=None, whitelist=(), remove_per_round=2, max_rounds=20, stop_at_n=None, ...) → IterativeCleanResult`

Iteratively realign (MAFFT) and remove the gappiest sequences each round.

Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `whitelist` | `()` | Accessions that are never removed |
| `remove_per_round` | `2` | Sequences removed each round |
| `stop_at_n` | `None` | Stop when this many sequences remain |
| `max_rounds` | `20` | Hard cap on number of rounds |
| `max_gap_fraction` | `1.0` | Only remove seqs above this gap fraction |
| `mafft_mode` | `"localpair"` | MAFFT mode (L-INS-i by default) |
| `use_trimal` | `True` | Trim columns with trimAl after each alignment |

`IterativeCleanResult.log_table` returns a formatted table of round-by-round statistics.

---

## Design notes

The whitelist mechanism was central to the original research workflow — the 17 manually verified flagellate GH7 sequences (BAB64565.3, AAY83390.3, etc.) had to be preserved through all cleaning rounds regardless of their gap content, because they were the biological signal of interest. The `protect_sequences()` function and `whitelist=` parameter throughout the API reflect this requirement.

The iterative cleaning algorithm mirrors `iterative_clean.sh`: align → trim → find worst 2 non-whitelisted sequences → remove → repeat. Translating this into Python makes it reproducible, testable, and independent of the user's shell environment.

---

## License

MIT
