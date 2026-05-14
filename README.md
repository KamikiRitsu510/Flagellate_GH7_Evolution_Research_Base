# Flagellate GH7 Evolution Research Base

Study of horizontal gene transfer and adaptive evolution of GH7 cellulases in flagellate gut symbionts of termites.

## Directory Structure

```
Flagellate_GH7_Evolution_Research_Base/
├── data/
│   ├── final_matrix.fasta          # Final alignment (534 sequences) — use for tree building
│   ├── sequences/                  # Curated sequence sets
│   │   ├── flagellate_only.fasta       (17 flagellate sequences)
│   │   ├── flagellate_GH7_protein.fasta
│   │   ├── bacterial_GH7_protein.fasta
│   │   ├── flagellate_CDS_all.fasta
│   │   ├── reference_CDS_BAB64565.fasta
│   │   ├── GH7_FINAL_525.fasta         (final 525-taxon set)
│   │   └── GH7_FINAL_PUB.fasta         (publication set)
│   ├── alignments/                 # MSA files
│   │   ├── GH7_alignment_mafft.fasta
│   │   ├── GH7_alignment_trimmed.fasta
│   │   ├── flagellate_codon_alignment.phy   (for selection analyses)
│   │   └── GH7_linsi_curated_final.fasta    (final curated LINSi alignment)
│   ├── metadata/                   # Annotation tables
│   │   ├── GH7_metadata_full.xlsx
│   │   ├── GH7_metadata_enriched_v4.xlsx
│   │   ├── GH7_Tables_Summary.xlsx
│   │   └── GH7_signal_peptide.xlsx
│   ├── structure/                  # 3D structure files
│   │   ├── fold_2026_04_22_17_16_model_0.cif   (AlphaFold model)
│   │   └── receptor.pdb
│   └── raw/                        # Raw/intermediate sequences
├── scripts/
│   ├── download_cds.py             # Fetch CDS sequences from NCBI
│   ├── get_taxonomy.py             # Taxonomy annotation
│   ├── highlight_residues.py       # Structure visualization (positive sites)
│   ├── clean_with_whitelist.py     # Sequence whitelist filtering
│   ├── find_mrca.py                # MRCA detection in tree
│   ├── pal2nal.pl                  # Protein-guided codon alignment
│   ├── fetch_species.sh
│   ├── iterative_clean.sh
│   └── visualization/              # Figure generation scripts (R)
│       ├── fig1.R
│       ├── figures_for_paper.R
│       ├── figures_with_legends.R
│       ├── generate_figures.R
│       └── [other plot_tree_*.R, final_tree_*.R, compare_*.R]
├── figures/
│   ├── publication/                # Final publication figures (Fig1–Fig6)
│   ├── supplementary/              # Supplementary figures
│   └── [working figures: trees, heatmaps, MSA, network, domain barcodes]
├── results/
│   ├── GH7_final_534.contree       # Bootstrap consensus tree (534 taxa)
│   ├── GH7_final_534.treefile      # ML best tree
│   ├── BUSTED_534_subset.json      # BUSTED positive selection results
│   ├── BUSTED_analysis_report.txt
│   ├── positive_sites_amino_acids.tsv
│   ├── signalp_results.csv
│   ├── GH7_diamond_edges.tsv
│   ├── GH7_network.rds
│   └── GH7_positive_sites.html
├── tools/
│   └── seqkit
└── logs/
    └── IQTREE_final.log
```

## Key Commands

### Build phylogenetic tree
```bash
iqtree -s data/final_matrix.fasta -m LG+F+I+G4 -bb 1000 -T AUTO
```

### Codon alignment (pal2nal)
```bash
perl scripts/pal2nal.pl data/sequences/flagellate_GH7_protein.fasta \
     data/sequences/flagellate_CDS_all.fasta -output fasta > data/alignments/codon_aln.fasta
```

### Positive selection (HyPhy BUSTED)
Input: `data/alignments/flagellate_codon_alignment_clean.phy` + pruned tree from `results/`

## Dataset Summary

| Dataset | Sequences | File |
|---------|-----------|------|
| Full alignment | 534 | data/final_matrix.fasta |
| Flagellates only | 17 | data/sequences/flagellate_only.fasta |
| Publication set | 525 | data/sequences/GH7_FINAL_525.fasta |

## Related Directories

- **Manuscript**: `../manuscript/` — paper drafts (v2, v3, v4)
- **Archive**: `../legacy/GH7_origin_pipeline/` — full pipeline history with all intermediate files
