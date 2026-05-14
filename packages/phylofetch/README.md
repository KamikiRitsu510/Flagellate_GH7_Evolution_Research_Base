# phylofetch

**Fetch CDS sequences and taxonomy lineages from NCBI — for phylogenetic studies.**

`phylofetch` wraps Biopython's Entrez API into a clean, rate-limited, retry-safe interface for batch retrieval of coding sequences and full taxonomy lineages. It was originally built to support GH7 cellulase evolution studies in termite gut flagellate symbionts, but works for any NCBI protein or nucleotide accession.

---

## Installation

```bash
pip install phylofetch
```

With optional `ete3` support (for tree handling):

```bash
pip install "phylofetch[full]"
```

**Requirements:** Python ≥ 3.8, Biopython ≥ 1.80, pandas ≥ 1.3

---

## Quick start

```python
from phylofetch import set_email, fetch_cds, fetch_taxonomy
from phylofetch import batch_fetch_cds, batch_fetch_taxonomy

# NCBI requires a valid email for Entrez access
set_email("you@example.com")

# ── Single accession ──────────────────────────────────────────────────────

seq = fetch_cds("BAB64565.3")
print(f"CDS length: {len(seq)} bp")

tax = fetch_taxonomy("BAB64565.3")
print(tax["phylum"], "/", tax["class"])

# ── Batch download ────────────────────────────────────────────────────────

accs = ["BAB64565.3", "AAY83390.3", "XP_006965674.1", "WP_165634028.1"]

# Save all CDS to a FASTA file
seqs = batch_fetch_cds(accs, output_fasta="my_sequences.fasta")
print(f"Downloaded {len(seqs)}/{len(accs)} sequences")

# Save taxonomy to Excel and TSV
df = batch_fetch_taxonomy(accs, output_excel="taxonomy.xlsx", output_tsv="taxonomy.tsv")
print(df[["accession", "phylum", "class", "genus", "species"]])
```

---

## Command-line interface

After installation, the `phylofetch` command is available:

```bash
# Fetch CDS for accessions listed in a file
phylofetch cds --input accs.txt --output sequences.fasta --email you@example.com

# Pass accessions directly
phylofetch cds BAB64565.3 AAY83390.3 --email you@example.com

# Fetch taxonomy → Excel + TSV
phylofetch taxonomy --input accs.txt \
    --excel taxonomy.xlsx --tsv taxonomy.tsv \
    --email you@example.com

# Include extra taxonomic ranks
phylofetch taxonomy --input accs.txt --extra-ranks subphylum tribe \
    --tsv full_taxonomy.tsv --email you@example.com
```

`accs.txt` is a plain text file with one accession per line (lines starting with `#` are ignored).

---

## API reference

### `set_email(email)`

Configure the email address sent to NCBI with every request (required by NCBI policy). Call once at the start of your script.

---

### `fetch_cds(accession, retries=3) → str | None`

Fetch the CDS nucleotide sequence for one accession.

The function tries the nucleotide database first. If nothing is found it falls back to the protein database and extracts CDS coordinates from the `coded_by` qualifier (format `NM_012345.1:100..500`).

Returns the sequence as a plain string (`ACGT…`), or `None` if not found.

---

### `batch_fetch_cds(accessions, output_fasta=None, retries=3, verbose=True) → dict`

Fetch CDS for a list of accessions.

| Parameter | Description |
|-----------|-------------|
| `accessions` | List of NCBI accession strings |
| `output_fasta` | Path to write a FASTA file (optional) |
| `retries` | Per-accession retry attempts (default 3) |
| `verbose` | Print progress to stdout |

Returns `{accession: sequence_string}` for successful fetches.

---

### `fetch_taxonomy(accession, extra_ranks=(), retries=3) → dict`

Fetch the full taxonomy lineage for one accession.

The returned dict contains:

| Key | Value |
|-----|-------|
| `accession` | The input accession |
| `taxid` | NCBI Taxonomy ID |
| `lineage_string` | Free-text lineage from GenBank |
| `superkingdom` | e.g. `"Eukaryota"` |
| `kingdom` | e.g. `"Fungi"` (or `None`) |
| `phylum` | e.g. `"Metamonada"` |
| `class` | e.g. `"Parabasalia"` |
| `order` | e.g. `"Trichomonadida"` |
| `family` | — |
| `genus` | — |
| `species` | Full binomial name |
| *extra ranks* | Any names passed via `extra_ranks` |

`None` for any rank absent in that organism's lineage.

---

### `batch_fetch_taxonomy(accessions, output_excel=None, output_tsv=None, extra_ranks=(), retries=3, verbose=True) → pd.DataFrame`

Fetch taxonomy for a list of accessions and return a tidy DataFrame (one row per accession). Optionally saves to Excel and/or TSV.

---

## Rate limiting and retries

`phylofetch` automatically respects NCBI's rate limits:

- **Without API key:** ≤ 3 requests/second (enforced via a 0.4 s minimum interval).
- **With API key:** set `Entrez.api_key = "your_key"` before calling any function; NCBI allows up to 10 req/s in that case, but the default limiter stays conservative.

All network calls are wrapped with exponential back-off retry logic (3 attempts, 2× back-off by default).

---

## Project background

Built during research on horizontal gene transfer (HGT) of GH7 cellulases from bacteria into termite gut flagellate symbionts (*Pseudotrichonympha*, *Holomastigotoides*, etc.). The original one-off scripts for bulk accession fetching were generalised here into a reusable package.

---

## License

MIT
