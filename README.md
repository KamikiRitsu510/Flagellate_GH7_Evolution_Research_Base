# Flagellate GH7 Evolution Research Base

Study of horizontal gene transfer and adaptive evolution of GH7 cellulases in flagellate gut symbionts of termites.

## Dataset Summary

| Dataset | Sequences | File |
|---------|-----------|------|
| **GH7 + GH12 outgroup (current)** | **539** (534 GH7 + 5 GH12) | `data/GH7_GH12_aligned_Pure.fasta` |
| GH7 + GH12 relaxed alignment | 539 | `data/GH7_GH12_relaxed.fasta` |
| GH7 + GH12 trimmed alignment | 539 | `data/GH7_GH12_trimmed.fasta` |
| GH7 only (prior run) | 534 | `data/final_matrix.fasta` |
| Flagellates only | 17 | `data/sequences/flagellate_only.fasta` |
| Publication set | 525 | `data/sequences/GH7_FINAL_525.fasta` |

> **Current main alignment**: `data/GH7_GH12_aligned_Pure.fasta` — 534 GH7 cellulase sequences plus 5 GH12 sequences used as outgroup, LINSi-aligned and manually curated.

## Directory Structure

```
Flagellate_GH7_Evolution_Research_Base/
├── data/
│   ├── GH7_GH12_aligned_Pure.fasta  # ★ CURRENT — 539 seqs (534 GH7 + 5 GH12 outgroup)
│   ├── GH7_GH12_relaxed.fasta       # Relaxed alignment version
│   ├── GH7_GH12_trimmed.fasta       # TrimAl-trimmed version
│   ├── final_matrix.fasta           # Previous 534-taxon GH7-only alignment
│   ├── sequences/                   # Curated sequence sets
│   │   ├── flagellate_only.fasta        (17 flagellate sequences)
│   │   ├── flagellate_GH7_protein.fasta
│   │   ├── bacterial_GH7_protein.fasta
│   │   ├── flagellate_CDS_all.fasta
│   │   ├── reference_CDS_BAB64565.fasta
│   │   ├── GH7_FINAL_525.fasta          (525-taxon publication set)
│   │   └── GH7_FINAL_PUB.fasta          (publication set)
│   ├── alignments/                  # MSA files
│   │   ├── GH7_alignment_mafft.fasta
│   │   ├── GH7_alignment_trimmed.fasta
│   │   ├── flagellate_codon_alignment.phy    (for selection analyses)
│   │   ├── flagellate_codon_alignment_clean.phy
│   │   └── GH7_linsi_curated_final.fasta     (final curated LINSi alignment)
│   ├── metadata/                    # Annotation tables
│   │   ├── GH7_metadata_full.xlsx
│   │   ├── GH7_metadata_enriched_v4.xlsx
│   │   ├── GH7_Tables_Summary.xlsx
│   │   └── GH7_signal_peptide.xlsx
│   ├── structure/                   # 3D structure files
│   │   ├── fold_2026_04_22_17_16_model_0.cif   (AlphaFold model)
│   │   └── receptor.pdb
│   └── raw/                         # Raw/intermediate sequences
├── scripts/
│   ├── plot_tree_clean.R            # ★ Current tree visualization (ggtree, species labels)
│   ├── compare_trees.R              # Side-by-side tree comparison (ggtree + patchwork)
│   ├── compare_trees_ape.R          # Side-by-side comparison (phytools fan trees)
│   ├── extract_subtree.R            # Extract subtree by accession/pattern
│   ├── download_cds.py              # Fetch CDS sequences from NCBI
│   ├── get_taxonomy.py              # Taxonomy annotation
│   ├── highlight_residues.py        # Structure visualization (positive sites)
│   ├── clean_with_whitelist.py      # Sequence whitelist filtering
│   ├── find_mrca.py                 # MRCA detection in tree
│   ├── pal2nal.pl                   # Protein-guided codon alignment
│   ├── fetch_species.sh
│   ├── iterative_clean.sh
│   └── visualization/               # Figure generation scripts (R), gitignored
├── packages/                        # Released bioinformatics packages
│   ├── treecolorize/                # R package — phylogenetic tree coloring
│   ├── phylofetch/                  # Python package — NCBI CDS/taxonomy fetch
│   └── msaclean/                    # Python package — MSA gap filtering
├── figures/                         # Output figures, gitignored
│   ├── publication/                 # Final publication figures (Fig1–Fig6)
│   └── supplementary/               # Supplementary figures
├── results/
│   ├── GH7_final_534.contree        # Bootstrap consensus tree (534-taxon GH7 run)
│   ├── GH7_final_534.treefile       # ML best tree
│   ├── GH7_Final_Rooted_Pure.contree # Rooted consensus tree (with GH12 outgroup)
│   ├── GH7_Final_Rooted_Pure.iqtree  # IQ-TREE report for rooted run
│   ├── BUSTED_534_subset.json       # BUSTED positive selection results
│   ├── BUSTED_analysis_report.txt
│   ├── positive_sites_amino_acids.tsv
│   ├── GH7_diamond_edges.tsv
│   ├── GH7_network.graphml
│   └── GH7_positive_sites.html
├── tools/
│   └── seqkit
└── logs/
    └── IQTREE_final.log
```

## Key Commands

### Build phylogenetic tree (with GH12 outgroup)
```bash
iqtree -s data/GH7_GH12_aligned_Pure.fasta -m LG+F+I+G4 -bb 1000 -T AUTO -o <GH12_accession>
```

### Visualize tree (R)
```r
source("scripts/plot_tree_clean.R")   # ggtree colored by taxonomic group
```

### Codon alignment (pal2nal)
```bash
perl scripts/pal2nal.pl data/sequences/flagellate_GH7_protein.fasta \
     data/sequences/flagellate_CDS_all.fasta -output fasta > data/alignments/codon_aln.fasta
```

### Positive selection (HyPhy BUSTED)
Input: `data/alignments/flagellate_codon_alignment_clean.phy` + pruned tree from `results/`

### Use packages
```r
# Install treecolorize
devtools::install_local("packages/treecolorize")
library(treecolorize)
treecolorize_help()
```
```bash
# Install phylofetch / msaclean
pip install -e packages/phylofetch
pip install -e packages/msaclean
phylofetch --help
msaclean --help
```

## Project Timeline

### Week of 2026-02-23 — Environment setup
- Downloaded and installed `seqkit` for sequence manipulation

### Week of 2026-04-13 — Raw data assembly
- Collected raw GH7 protein sequences from NCBI; ran CD-HIT clustering (`GH7_nr95`)
- First IQ-TREE run on `gh7_aligned.fasta` (legacy pipeline)
- Early figure prototypes in `legacy/sequence_processing/`

### Week of 2026-04-17 — First IQ-TREE production run
- Committed raw sequences (`data/raw/GH7_all_raw.fasta`, `GH7_cdhit_clustered.fasta`)
- Produced initial tree figures: `GH7_tree.pdf`, `GH7_rect_final.pdf`, `GH7_tree_large.pdf`
- `IQTREE_final.log` records the first full production run

### Week of 2026-04-20 — Tree visualization iteration
- Switched to ggtree + ape dual-backend coloring
- Added species name display (stripping accession prefix); colored by taxonomic group
- Figures: `GH7_tree_ape_colored.pdf`, `GH7_tree_final_species_names.pdf`, `GH7_tree_full_species.pdf`

### Week of 2026-04-21–23 — Network analysis & publication figures
- DIAMOND all-vs-all → `GH7_diamond_edges.tsv` → network in `GH7_network.rds` / `GH7_network.png`
- Supplementary figures generated: Fig1 (phylum pie), Fig2 (signal peptide), Fig3 (domain), Fig4 (heatmap), Fig5 (workflow)
- Publication figures finalized: Fig1 phylum bar, Fig4–6 (signal peptide, domain architecture, heatmap)
- AlphaFold structure predicted (`data/structure/fold_2026_04_22_17_16_model_0.cif`)

### Week of 2026-04-27 — Positive sites & domain analysis
- Mapped BUSTED positive selection sites onto AlphaFold structure (`positive_sites_structure.png`)
- Domain barcode figures: `GH7_domain_barcode.pdf/png`
- LINSi re-alignment: `GH7_tree_LINSi_colored.pdf`, `GH7_tree_LINSi_species_colored.pdf`
- First tree comparison figures: `tree_comparison_species.pdf`

### Week of 2026-05-04 — MSA curation & whitelist heatmaps
- Iterative MSA cleaning with whitelisted flagellate sequences
- Heatmap visualization of gap patterns: `GH7_MSA_heatmap.pdf`, `GH7_MSA_curated_clean.pdf`
- `GH7_whitelist_final_heatmap.pdf`, `Matrix_Optimization_Summary.pdf`
- Near-final MSA: `GH7_Final_MSA_538.pdf`, `GH7_MSA_heatmap_curated_v3.pdf`

### Week of 2026-05-11 — CDS subset & tree comparison scripts
- Extracted CDS subset for codon-aware analyses (`data/raw/cds_subset_*.fasta`)
- Developed multi-iteration tree comparison visualization scripts in `scripts/visualization/`
- Final tree polishing: `GH7_tree_final_pretty.pdf`, `GH7_11tree_complete.pdf`, `GH7_Complete_Tree.pdf`
- `plot_tree_with_species_names.R` → cleaned version becomes `scripts/plot_tree_clean.R`

### Week of 2026-05-15 — Final dataset & package release *(current)*
- **Dataset finalized**: 534 GH7 cellulase sequences + 5 GH12 sequences as outgroup = **539 taxa**
  - Main file: `data/GH7_GH12_aligned_Pure.fasta` (LINSi-aligned, manually curated)
  - Rooted IQ-TREE results: `results/GH7_Final_Rooted_Pure.contree`
- **Three packages released** (`packages/`):
  - `treecolorize` (R) — keyword-based tip classification, ape + ggtree backends, subtree extraction, tree comparison
  - `phylofetch` (Python) — batch NCBI CDS/taxonomy download with Entrez API
  - `msaclean` (Python) — gap-fraction filtering with whitelist support and iterative MAFFT re-alignment
- `scripts/plot_tree_clean.R` finalized — ggtree with species-name labels and group color palette

## Related Directories

- **Manuscript**: `../manuscript/` — paper drafts (v2, v3, v4)
- **Archive**: `../legacy/GH7_origin_pipeline/` — full pipeline history with all intermediate files
