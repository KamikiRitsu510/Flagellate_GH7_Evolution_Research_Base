#' Extract a subtree by keeping tips that match given patterns
#'
#' A convenience wrapper around [ape::keep.tip()] / [ape::drop.tip()] that
#' lets you specify which tips to keep by a character vector of regex patterns,
#' exact accession strings, or a combination of both.
#'
#' This function was generalised from the GH7 flagellate subtree extraction
#' workflow (`extract_subtree.R`), where specific accession prefixes and
#' genus names were used to pull flagellate + outgroup sequences from a
#' 534-tip tree.
#'
#' @param tree A `"phylo"` object.
#' @param patterns Character vector of regex patterns. Any tip whose label
#'   matches at least one pattern (case-insensitive by default) is kept.
#' @param exact Character vector of exact tip labels to keep in addition
#'   to pattern matches (e.g. specific accession numbers like
#'   `"XP_006965674.1"`). Labels not present in the tree are silently ignored.
#' @param ignore_case Logical. If `TRUE` (default), pattern matching is
#'   case-insensitive.
#' @param output_file If not `NULL`, writes the subtree to this file in
#'   Newick format via [ape::write.tree()].
#' @param verbose If `TRUE` (default), prints the original and retained
#'   tip counts.
#' @return A `"phylo"` object containing only the matched tips.
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- read.tree("results/GH7_final_534.contree")
#'
#' # Keep flagellate sequences + specific outgroup accessions
#' subtree <- extract_subtree(
#'   tree,
#'   patterns = c("AAY833", "BAB645", "BAC075",   # accession prefixes
#'                "Daphnia", "pulex"),              # outgroup by genus
#'   exact    = c("XP_006965674.1",                # specific accessions
#'                "XP_003054277.1",
#'                "AGT15838.1"),
#'   output_file = "results/subtree_flagellate.nwk"
#' )
#' cat("Retained:", length(subtree$tip.label), "tips\n")
#' }
#' @export
extract_subtree <- function(tree,
                            patterns    = character(0),
                            exact       = character(0),
                            ignore_case = TRUE,
                            output_file = NULL,
                            verbose     = TRUE) {
  if (!inherits(tree, "phylo")) stop("`tree` must be a 'phylo' object.")

  tips <- tree$tip.label

  # Tips matched by regex patterns
  pattern_hits <- if (length(patterns) > 0) {
    combined_pattern <- paste(patterns, collapse = "|")
    tips[grepl(combined_pattern, tips, ignore.case = ignore_case)]
  } else {
    character(0)
  }

  # Tips matched by exact label
  exact_hits <- intersect(exact, tips)
  missing    <- setdiff(exact, tips)
  if (length(missing) > 0) {
    warning(sprintf("%d exact label(s) not found in tree: %s",
                    length(missing), paste(missing, collapse = ", ")))
  }

  keep <- unique(c(pattern_hits, exact_hits))

  if (length(keep) == 0) stop("No tips matched; check your patterns / exact labels.")

  subtree <- ape::keep.tip(tree, keep)

  if (verbose) {
    cat(sprintf("Original tips : %d\n", length(tips)))
    cat(sprintf("Retained tips : %d\n", length(subtree$tip.label)))
    cat("Retained labels:\n")
    cat(paste(" ", sort(subtree$tip.label), collapse = "\n"), "\n")
  }

  if (!is.null(output_file)) {
    ape::write.tree(subtree, file = output_file)
    if (verbose) cat(sprintf("Subtree saved : %s\n", output_file))
  }

  subtree
}
