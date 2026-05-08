# Flagellate GH7 Evolution Research Base

Final data for GH7 cellulase evolution study.

## Files
- `data/final_matrix.fasta`: 534 sequences alignment
- `results/`: IQ-TREE output files (.contree, .treefile, .iqtree)

## Build tree
```
iqtree -s data/final_matrix.fasta -m LG+F+I+G4 -bb 1000 -T AUTO
```
