#' One-step colored tree plot
#'
#' A high-level wrapper that runs [classify_tips()], [auto_palette()], and
#' your choice of [plot_tree_ape()] or [plot_tree_ggtree()] in one call.
#' Suitable for quick exploration; for publication-quality output, call the
#' lower-level functions individually for full control.
#'
#' @param tree An object of class `"phylo"`.
#' @param keyword_map Named list mapping regex patterns to group names.
#'   If `NULL`, [default_keyword_map()] is used.
#' @param backend Plotting backend: `"ape"` (default, requires only **ape**)
#'   or `"ggtree"` (requires Bioconductor **ggtree** + **ggplot2**).
#' @param palette Optional named color vector. If `NULL`, colors are
#'   assigned automatically via [auto_palette()].
#' @param output_file If not `NULL`, saves the plot to this path.
#'   - For `backend = "ape"`: must be a `.pdf` path.
#'   - For `backend = "ggtree"`: any format accepted by [ggplot2::ggsave()]
#'     (e.g. `.pdf`, `.png`, `.svg`).
#' @param width Plot width (inches) when saving. Default: `12`.
#' @param height Plot height (inches) when saving. Default: `28`.
#' @param ... Passed to [plot_tree_ape()] or [plot_tree_ggtree()].
#' @return
#'   - `backend = "ape"`: Invisibly returns a list with `tip_groups`,
#'     `palette`, and `edge_colors`.
#'   - `backend = "ggtree"`: Returns the `ggplot` object.
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- read.tree("GH7_final_534.contree")
#'
#' # Quick ape plot
#' treecolorize(tree)
#'
#' # Custom keyword map, ggtree backend, save to PDF
#' kmap <- list(
#'   "pseudotrichonympha|holomastigotoides|trichonympha" = "Flagellate",
#'   "chaetomium|aspergillus" = "Fungi",
#'   "cellulomonas|bacillus"  = "Bacteria"
#' )
#' treecolorize(tree,
#'              keyword_map  = kmap,
#'              backend      = "ggtree",
#'              output_file  = "tree_colored.pdf",
#'              layout       = "fan",
#'              title        = "GH7 cellulase phylogeny")
#' }
#' @export
treecolorize <- function(tree,
                         keyword_map  = NULL,
                         backend      = c("ape", "ggtree"),
                         palette      = NULL,
                         output_file  = NULL,
                         width        = 12,
                         height       = 28,
                         ...) {
  backend <- match.arg(backend)

  if (is.null(keyword_map)) keyword_map <- default_keyword_map()

  tip_groups <- classify_tips(tree$tip.label, keyword_map)

  cat(sprintf("Classification summary (%d tips):\n", length(tip_groups)))
  print(table(tip_groups))
  cat("\n")

  if (is.null(palette)) palette <- auto_palette(tip_groups)

  if (backend == "ape") {
    if (!is.null(output_file)) {
      result <- save_tree_ape(tree, tip_groups,
                               file   = output_file,
                               palette = palette,
                               width  = width,
                               height = height,
                               ...)
      cat(sprintf("Tree saved to: %s\n", output_file))
      invisible(result)
    } else {
      plot_tree_ape(tree, tip_groups, palette = palette, ...)
    }
  } else {
    p <- plot_tree_ggtree(tree, tip_groups, palette = palette, ...)
    if (!is.null(output_file)) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required: install.packages('ggplot2')")
      }
      ggplot2::ggsave(output_file, p, width = width, height = height,
                      limitsize = FALSE)
      cat(sprintf("Tree saved to: %s\n", output_file))
    } else {
      print(p)
    }
    invisible(p)
  }
}
