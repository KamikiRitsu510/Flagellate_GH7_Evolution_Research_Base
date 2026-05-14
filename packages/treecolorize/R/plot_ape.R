#' Plot a colored phylogenetic tree using ape
#'
#' Draws a phylogenetic tree with tip labels and branches colored by group,
#' using [ape::plot.phylo()] as the rendering backend. This is the
#' lightweight option that requires only **ape** (no Bioconductor).
#'
#' @param tree An object of class `"phylo"`.
#' @param tip_groups Named character vector of group assignments for each
#'   tip. See [classify_tips()].
#' @param palette Named character vector of hex colors. If `NULL`,
#'   [auto_palette()] is called automatically.
#' @param type Tree layout type passed to [ape::plot.phylo()].
#'   One of `"phylogram"`, `"cladogram"`, `"fan"`, `"unrooted"`,
#'   `"radial"`. Default: `"phylogram"`.
#' @param branch_method How to color internal branches. `"majority"`
#'   (default) uses [branch_colors_majority()]; `"tip_only"` colors only
#'   terminal edges; `"uniform"` uses a single `branch_color`.
#' @param branch_color Hex color for branches when `branch_method =
#'   "uniform"`. Default: `"#333333"`.
#' @param show_tip_labels Logical. Show tip labels? Default: `TRUE`.
#' @param tip_label_fn Function applied to tip labels before display.
#'   Default: [strip_accession()], so labels show species names only.
#'   Pass `identity` to show raw tip labels.
#' @param tip_cex Character expansion for tip labels. Default: `0.5`.
#' @param tip_col_alpha Alpha transparency for tip label colors (0–1).
#'   Default: `1.0`.
#' @param edge_width Line width for branches. Default: `0.8`.
#' @param legend Logical. Add a group legend? Default: `TRUE`.
#' @param legend_pos Legend position. One of `"topright"`, `"topleft"`,
#'   `"bottomright"`, `"bottomleft"`. Default: `"topright"`.
#' @param title Optional title string.
#' @param mixed_color Color for mixed-group internal edges. Default:
#'   `"#CCCCCC"`.
#' @param ... Additional arguments passed to [ape::plot.phylo()].
#' @return Invisibly returns a list with `tip_groups`, `palette`, and
#'   `edge_colors` so downstream code can reuse them.
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- read.tree("results/GH7_final_534.contree")
#' kmap <- list(
#'   "holomastigotoides|trichonympha|pseudotrichonympha" = "Flagellate",
#'   "chaetomium|aspergillus|fungi"                      = "Fungi",
#'   "cellulomonas|bacillus"                             = "Bacteria"
#' )
#' groups <- classify_tips(tree$tip.label, kmap)
#' plot_tree_ape(tree, groups)
#' }
#' @export
plot_tree_ape <- function(tree,
                          tip_groups,
                          palette          = NULL,
                          type             = "phylogram",
                          branch_method    = c("majority", "tip_only", "uniform"),
                          branch_color     = "#333333",
                          show_tip_labels  = TRUE,
                          tip_label_fn     = strip_accession,
                          tip_cex          = 0.5,
                          tip_col_alpha    = 1.0,
                          edge_width       = 0.8,
                          legend           = TRUE,
                          legend_pos       = "topright",
                          title            = NULL,
                          mixed_color      = "#CCCCCC",
                          ...) {
  branch_method <- match.arg(branch_method)

  # Build palette if not supplied
  if (is.null(palette)) palette <- auto_palette(tip_groups)

  # Ensure all tips are covered
  tip_groups_ordered <- tip_groups[tree$tip.label]
  if (anyNA(tip_groups_ordered)) {
    stop("Some tips have no group assignment. Check `tip_groups` names.")
  }

  # Tip label display and color
  tip_labels <- if (show_tip_labels) tip_label_fn(tree$tip.label) else rep("", ape::Ntip(tree))
  tip_cols   <- grDevices::adjustcolor(
    palette[tip_groups_ordered],
    alpha.f = tip_col_alpha
  )

  # Branch colors
  edge_cols <- switch(branch_method,
    majority  = branch_colors_majority(tree, tip_groups, palette,
                                        mixed_color = mixed_color),
    tip_only  = branch_colors_tip_only(tree, tip_groups, palette,
                                        mixed_color = mixed_color),
    uniform   = rep(branch_color, nrow(tree$edge))
  )

  # Draw
  ape::plot.phylo(
    tree,
    type            = type,
    show.tip.label  = show_tip_labels,
    tip.color       = tip_cols,
    edge.color      = edge_cols,
    edge.width      = edge_width,
    cex             = tip_cex,
    label.offset    = 0.5,
    ...
  )

  if (!is.null(title)) title(main = title)

  if (legend) {
    present_groups <- intersect(names(palette), unique(tip_groups_ordered))
    legend(
      legend_pos,
      legend  = present_groups,
      col     = palette[present_groups],
      lwd     = 3,
      cex     = 0.7,
      bty     = "n",
      title   = "Group"
    )
  }

  invisible(list(tip_groups = tip_groups_ordered,
                 palette    = palette,
                 edge_colors= edge_cols))
}


#' Save a colored ape tree to a PDF file
#'
#' A convenience wrapper around [plot_tree_ape()] that opens a PDF device,
#' draws the tree, and closes the device.
#'
#' @inheritParams plot_tree_ape
#' @param file Output PDF file path.
#' @param width PDF width in inches. Default: `12`.
#' @param height PDF height in inches. Default: `28`.
#' @param mar Plot margins (bottom, left, top, right). Default: `c(1,1,2,8)`.
#' @return Invisibly returns the path to the saved file.
#' @export
save_tree_ape <- function(tree, tip_groups, file,
                          palette = NULL,
                          width = 12, height = 28,
                          mar   = c(1, 1, 2, 8),
                          ...) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off())
  par(mar = mar)
  result <- plot_tree_ape(tree, tip_groups, palette = palette, ...)
  invisible(file)
}
