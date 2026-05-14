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
#' Print a quick-reference help summary for treecolorize
#'
#' Prints a concise function reference to the console — useful when you
#' forget a parameter name and don't want to open the full documentation.
#' For full docs, use `?treecolorize`, `?classify_tips`, etc.
#'
#' @return Invisibly returns `NULL`. Called for its side effect.
#' @examples
#' treecolorize_help()
#' @export
treecolorize_help <- function() {
  msg <- c(
    "",
    "\033[1m── treecolorize 快速参考 ──────────────────────────────────────────────\033[0m",
    "",
    "\033[1m主函数（一站式）:\033[0m",
    "  treecolorize(tree, keyword_map = NULL, backend = 'ape',",
    "               palette = NULL, output_file = NULL, width = 12, height = 28, ...)",
    "  → 自动完成：分类 → 配色 → 绘图（→ 可选保存）",
    "",
    "\033[1m分步函数:\033[0m",
    "  classify_tips(tips, keyword_map, other_label = 'Other', strip_prefix = TRUE)",
    "    → 按关键词映射表给 tip 分组，返回命名字符向量",
    "",
    "  auto_palette(groups, user_colors = NULL, other_color = '#AAAAAA')",
    "    → 自动配色（Wong 2011 色盲友好调色板），返回命名颜色列表",
    "",
    "  plot_tree_ape(tree, tip_groups, palette = NULL,",
    "                type = 'phylogram', branch_method = 'majority',",
    "                tip_cex = 0.5, edge_width = 0.8, legend = TRUE, ...)",
    "    → 用 ape 出图（直接显示，不返回对象）",
    "",
    "  save_tree_ape(tree, tip_groups, file, palette = NULL,",
    "                width = 12, height = 28, ...)",
    "    → 保存为 PDF/PNG",
    "",
    "  plot_tree_ggtree(tree, tip_groups, palette = NULL,",
    "                   layout = 'rectangular', tip_label_size = 1.5, ...)",
    "    → 用 ggtree 绘图，返回 ggplot 对象（可继续 + 图层）",
    "",
    "\033[1mbranch_method 选项:\033[0m",
    "  'majority'  后序遍历，内部节点按子节点多数颜色着色（默认）",
    "  'tip_only'  仅 tip 着色，内部枝条统一 mixed_color",
    "  'uniform'   所有枝条同色（branch_color 参数指定）",
    "",
    "\033[1m辅助函数:\033[0m",
    "  strip_accession(labels)      去掉 tip label 的 accession 前缀",
    "  default_keyword_map()        内置常见界/门关键词映射表",
    "  branch_colors_majority(...)  手动计算内部节点颜色向量",
    "",
    "\033[1m完整文档:\033[0m  ?treecolorize  |  ?classify_tips  |  ?plot_tree_ape",
    "\033[1m示例教程:\033[0m  vignette('treecolorize-example')  或查看 EXAMPLE.md",
    "\033[1m────────────────────────────────────────────────────────────────────\033[0m",
    ""
  )
  cat(paste(msg, collapse = "\n"))
  invisible(NULL)
}


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
