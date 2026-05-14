#' Plot a colored phylogenetic tree using ggtree
#'
#' Produces a ggplot2-compatible tree plot where tip points and labels are
#' colored by group. Requires the **ggtree** and **ggplot2** packages
#' (Bioconductor).
#'
#' @param tree An object of class `"phylo"` or any class supported by
#'   [treeio::read.tree()].
#' @param tip_groups Named character vector of group assignments for each
#'   tip. See [classify_tips()].
#' @param palette Named character vector of hex colors. If `NULL`,
#'   [auto_palette()] is called automatically.
#' @param layout ggtree layout. One of `"rectangular"`, `"slanted"`,
#'   `"fan"`, `"circular"`, `"radial"`, `"equal_angle"`, `"daylight"`.
#'   Default: `"rectangular"`.
#' @param show_tip_labels Logical. Default: `TRUE`.
#' @param tip_label_fn Function applied to tip labels before display.
#'   Default: [strip_accession()].
#' @param tip_label_size Font size for tip labels. Default: `1.5`.
#' @param tip_point_size Size of tip point symbols. Default: `0.8`.
#' @param label_offset Horizontal offset of tip labels from tip points.
#'   Default: `0.5`.
#' @param legend_position ggplot2 legend position. Default: `"right"`.
#' @param title Optional plot title string.
#' @param ... Additional layers or theme modifications added to the ggplot
#'   object (e.g. `theme(plot.background = element_rect(fill = "white"))`).
#' @return A `ggplot` object. Use [ggplot2::ggsave()] to save it.
#' @examples
#' \dontrun{
#' library(ape); library(ggtree)
#' tree   <- read.tree("results/GH7_final_534.contree")
#' groups <- classify_tips(tree$tip.label, default_keyword_map())
#' p <- plot_tree_ggtree(tree, groups, layout = "fan", title = "GH7 phylogeny")
#' ggplot2::ggsave("tree.pdf", p, width = 14, height = 14)
#' }
#' @export
plot_tree_ggtree <- function(tree,
                             tip_groups,
                             palette          = NULL,
                             layout           = "rectangular",
                             show_tip_labels  = TRUE,
                             tip_label_fn     = strip_accession,
                             tip_label_size   = 1.5,
                             tip_point_size   = 0.8,
                             label_offset     = 0.5,
                             legend_position  = "right",
                             title            = NULL,
                             ...) {
  # Check Bioconductor dependencies
  if (!requireNamespace("ggtree",  quietly = TRUE)) {
    stop("Package 'ggtree' required. Install with: BiocManager::install('ggtree')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }

  if (is.null(palette)) palette <- auto_palette(tip_groups)

  # Build annotation data frame
  tip_labels_display <- tip_label_fn(names(tip_groups))
  df <- data.frame(
    label   = names(tip_groups),
    group   = unname(tip_groups),
    display = tip_labels_display,
    stringsAsFactors = FALSE
  )

  # Build the plot
  p <- ggtree::ggtree(tree, layout = layout,
                      branch.length = "branch.length") %<+% df

  p <- p + ggtree::geom_tippoint(
    ggplot2::aes(color = group),
    size = tip_point_size
  )

  if (show_tip_labels) {
    p <- p + ggtree::geom_tiplab(
      ggplot2::aes(label = display, color = group),
      size    = tip_label_size,
      offset  = label_offset,
      align   = TRUE,
      linesize = 0.15
    )
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = palette,
      name   = "Group",
      breaks = intersect(names(palette), unique(tip_groups))
    ) +
    ggplot2::theme(
      legend.position = legend_position,
      plot.margin     = ggplot2::unit(c(0.5, 3, 0.5, 0.5), "cm")
    )

  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title)
  }

  # Allow caller to add extra layers
  extras <- list(...)
  for (layer in extras) {
    p <- p + layer
  }

  p
}
