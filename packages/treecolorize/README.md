# treecolorize

Keyword-based coloring of phylogenetic trees, with **ape** and **ggtree** backends.

## The problem it solves

When you have a phylogeny containing sequences from many kingdoms (fungi, bacteria, flagellates, …), coloring the tips and branches by group requires three steps that every paper rewrites from scratch:

1. Parse tip labels to extract species names  
2. Classify species by regex keyword rules  
3. Color branches by the *majority group among descendants*

`treecolorize` packages these into a clean, configurable API.

## Installation

```r
# Install from GitHub (once the repo is public)
# install.packages("remotes")
remotes::install_github("KamikiRitsu510/treecolorize")
```

You only need **ape** for the default backend.  
For the ggtree backend, also install:

```r
BiocManager::install("ggtree")
install.packages("ggplot2")
```

## Quick start

```r
library(ape)
library(treecolorize)

tree <- read.tree("GH7_final_534.contree")

# One-liner with sensible defaults
treecolorize(tree)

# Custom keyword rules
kmap <- list(
  "holomastigotoides|trichonympha|pseudotrichonympha" = "Flagellate",
  "chaetomium|aspergillus|trichoderma|fungi"          = "Fungi",
  "cellulomonas|bacillus|clostridium"                 = "Bacteria",
  "daphnia|pulex"                                     = "Crustacea"
)

treecolorize(tree,
             keyword_map = kmap,
             backend     = "ggtree",
             layout      = "fan",
             output_file = "tree_colored.pdf",
             width = 14, height = 14)
```

## Step-by-step API

```r
# 1. Classify tips
groups <- classify_tips(tree$tip.label, kmap)
table(groups)
#> Bacteria  Flagellate  Fungi  Other
#>        5          17    492     20

# 2. Build color palette
pal <- auto_palette(groups,
                    user_colors = c(Flagellate = "#CC0000"))

# 3a. Plot with ape (lightweight, no Bioconductor)
plot_tree_ape(tree, groups, palette = pal,
              type          = "phylogram",
              branch_method = "majority",   # color branches by descendant majority
              tip_cex       = 0.4,
              title         = "GH7 cellulase evolution")

# 3b. Plot with ggtree (ggplot2 ecosystem)
p <- plot_tree_ggtree(tree, groups, palette = pal,
                      layout         = "rectangular",
                      tip_label_size = 1.3)
ggplot2::ggsave("tree.pdf", p, width = 14, height = 28, limitsize = FALSE)

# 4. Get edge colors for custom ape plotting
ecols <- branch_colors_majority(tree, groups, pal)
plot(tree, edge.color = ecols, edge.width = 1.2, show.tip.label = FALSE)
```

## Key design decisions

| Feature | Detail |
|---------|--------|
| **Keyword order** | First match wins; put specific patterns before general ones |
| **Branch coloring** | Plurality vote among direct children (not all descendants) |
| **Mixed edges** | Edges where no group has ≥ 50 % are colored `#CCCCCC` by default |
| **Tip label parsing** | First whitespace-delimited token is stripped (accession removed) |
| **Color palette** | Colorblind-friendly (Wong 2011) base set; auto-extends with `rainbow()` for > 12 groups |

## Functions

| Function | Description |
|----------|-------------|
| `treecolorize()` | One-step wrapper |
| `classify_tips()` | Keyword → group classification |
| `strip_accession()` | Remove accession prefix from labels |
| `default_keyword_map()` | Built-in map for common kingdoms |
| `auto_palette()` | Assign colors to groups |
| `branch_colors_majority()` | Edge colors by descendant majority vote |
| `branch_colors_tip_only()` | Edge colors for terminal edges only |
| `plot_tree_ape()` | ape backend |
| `plot_tree_ggtree()` | ggtree backend |
| `save_tree_ape()` | Save ape plot to PDF |

## Citation

If you use treecolorize in a publication, please cite:

> KamikiRitsu (2026). treecolorize: Keyword-Based Coloring of Phylogenetic Trees.
> R package version 0.1.0. https://github.com/KamikiRitsu510/treecolorize
