#' Compare two phylogenetic trees side by side (ggtree backend)
#'
#' Draws two trees in a side-by-side panel using ggtree and patchwork.
#' If the two trees have different tip sets, only shared tips are kept
#' before plotting (with a warning).
#'
#' @param old_tree A `"phylo"` object — the reference / older tree.
#' @param new_tree A `"phylo"` object — the updated / newer tree.
#' @param old_label Title for the left panel. Default: `"Old tree"`.
#' @param new_label Title for the right panel. Default: `"New tree"`.
#' @param layout ggtree layout passed to both panels. Default: `"circular"`.
#' @param tip_label_size Font size for tip labels. Default: `1.5`.
#' @param output_file If not `NULL`, saves to this path via [ggplot2::ggsave()].
#' @param width Plot width in inches when saving. Default: `16`.
#' @param height Plot height in inches when saving. Default: `8`.
#' @return A patchwork / ggplot object (invisibly when `output_file` is set).
#' @examples
#' \dontrun{
#' library(ape)
#' old <- read.tree("results/11_GH7_tree_consensus.contree")
#' new <- read.tree("results/GH7_final_534.contree")
#' compare_trees(old, new,
#'               old_label = "GH7 v1 (consensus)",
#'               new_label = "GH7 v2 (534 seq)",
#'               output_file = "trees_comparison.pdf")
#' }
#' @export
compare_trees <- function(old_tree,
                          new_tree,
                          old_label      = "Old tree",
                          new_label      = "New tree",
                          layout         = "circular",
                          tip_label_size = 1.5,
                          output_file    = NULL,
                          width          = 16,
                          height         = 8) {
  for (pkg in c("ggtree", "ggplot2", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' required. Install with: %s",
                   pkg,
                   if (pkg == "ggtree") "BiocManager::install('ggtree')"
                   else sprintf("install.packages('%s')", pkg)))
    }
  }

  old_tree <- .reconcile_tips(old_tree, new_tree)$t1
  new_tree <- .reconcile_tips(old_tree, new_tree)$t2

  p1 <- ggtree::ggtree(old_tree, layout = layout,
                        branch.length = "branch.length") +
    ggtree::geom_tiplab(size = tip_label_size) +
    ggplot2::labs(title = old_label) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p2 <- ggtree::ggtree(new_tree, layout = layout,
                        branch.length = "branch.length") +
    ggtree::geom_tiplab(size = tip_label_size) +
    ggplot2::labs(title = new_label) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  combined <- p1 + p2 +
    patchwork::plot_annotation(title = paste(old_label, "vs", new_label))

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, combined,
                    width = width, height = height, limitsize = FALSE)
    message("Saved: ", output_file)
    return(invisible(combined))
  }
  combined
}


#' Compare two phylogenetic trees side by side (ape / phytools backend)
#'
#' Draws two fan trees in a side-by-side PDF panel using base graphics
#' (ape + phytools). Does not require Bioconductor packages.
#'
#' @param old_tree A `"phylo"` object.
#' @param new_tree A `"phylo"` object.
#' @param old_label Title for the left panel.
#' @param new_label Title for the right panel.
#' @param output_file PDF path. If `NULL`, draws to the current device.
#' @param width PDF width in inches. Default: `12`.
#' @param height PDF height in inches. Default: `8`.
#' @param fsize Font size for tip labels (passed to [phytools::plotTree()]).
#'   Default: `0.6`.
#' @param type Tree type: `"fan"` (default) or `"phylogram"`.
#' @return Invisibly `NULL`. Called for side effects.
#' @examples
#' \dontrun{
#' library(ape)
#' old <- read.tree("results/11_GH7_tree_consensus.contree")
#' new <- read.tree("results/GH7_final_534.contree")
#' compare_trees_ape(old, new,
#'                   old_label = "GH7 v1", new_label = "GH7 v2",
#'                   output_file = "trees_comparison_ape.pdf")
#' }
#' @export
compare_trees_ape <- function(old_tree,
                               new_tree,
                               old_label   = "Old tree",
                               new_label   = "New tree",
                               output_file = NULL,
                               width       = 12,
                               height      = 8,
                               fsize       = 0.6,
                               type        = "fan") {
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' required: install.packages('phytools')")
  }

  reconciled <- .reconcile_tips(old_tree, new_tree)
  old_tree   <- reconciled$t1
  new_tree   <- reconciled$t2

  if (!is.null(output_file)) grDevices::pdf(output_file, width = width, height = height)

  graphics::par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
  phytools::plotTree(old_tree, type = type, fsize = fsize, lwd = 1,
                     ftype = "i", offset = 0.5)
  graphics::title(old_label, cex.main = 1.2)

  phytools::plotTree(new_tree, type = type, fsize = fsize, lwd = 1,
                     ftype = "i", offset = 0.5)
  graphics::title(new_label, cex.main = 1.2)

  if (!is.null(output_file)) {
    grDevices::dev.off()
    message("Saved: ", output_file)
  }
  invisible(NULL)
}


# ── Internal helper ────────────────────────────────────────────────────────────

.reconcile_tips <- function(t1, t2) {
  common <- intersect(t1$tip.label, t2$tip.label)
  if (length(common) < max(length(t1$tip.label), length(t2$tip.label))) {
    n_diff <- length(t1$tip.label) + length(t2$tip.label) - 2 * length(common)
    warning(sprintf(
      "%d tip(s) not shared between trees; using the %d common tips.",
      n_diff, length(common)
    ))
    t1 <- ape::keep.tip(t1, common)
    t2 <- ape::keep.tip(t2, common)
  }
  list(t1 = t1, t2 = t2)
}
