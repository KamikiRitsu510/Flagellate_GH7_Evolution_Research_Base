#' Compute edge colors by descendant-majority vote
#'
#' For each edge in the tree, determines the color of the dominant group
#' among all tips descended from the child node. Edges that subtend a mix
#' of groups receive the color of the plurality group; ties are broken by
#' the order of `palette`.
#'
#' This function is the core of the "branch coloring" approach: it lets
#' monophyletic clades stand out in a single color while mixed-origin
#' nodes fall back to their plurality group.
#'
#' @param tree An object of class `"phylo"` (from the **ape** package).
#' @param tip_groups Named character vector mapping each tip label to a
#'   group name. Names must match `tree$tip.label`. See [classify_tips()].
#' @param palette Named character vector of hex colors, one per group.
#'   See [auto_palette()].
#' @param mixed_color Hex color for edges where no group has a strict
#'   majority (> 50 %). Default: `"#CCCCCC"`. Set to `NULL` to use the
#'   plurality color even for mixed edges.
#' @param mixed_threshold Minimum proportion (0–1) for a group to "win"
#'   a node. Below this, `mixed_color` is used. Default: `0.5`.
#' @return A character vector of hex colors with length equal to
#'   `nrow(tree$edge)`, one color per edge.
#' @examples
#' library(ape)
#' tree   <- rtree(20)
#' groups <- setNames(sample(c("A","B","C"), 20, replace=TRUE), tree$tip.label)
#' pal    <- auto_palette(groups)
#' ecols  <- branch_colors_majority(tree, groups, pal)
#' plot(tree, edge.color = ecols, edge.width = 2)
#' @export
branch_colors_majority <- function(tree,
                                   tip_groups,
                                   palette,
                                   mixed_color      = "#CCCCCC",
                                   mixed_threshold  = 0.5) {
  n_tips  <- ape::Ntip(tree)
  n_nodes <- ape::Nnode(tree)
  n_all   <- n_tips + n_nodes

  # Validate tip_groups alignment
  matched <- tip_groups[tree$tip.label]
  if (anyNA(matched)) {
    missing <- tree$tip.label[is.na(matched)]
    stop("No group assignment for tips: ", paste(missing[1:min(5, length(missing))], collapse=", "))
  }

  # node_group[i] stores the *winning* group string at node i
  # (indices 1..n_tips = tips, (n_tips+1)..n_all = internal nodes)
  node_group <- character(n_all)
  node_group[seq_len(n_tips)] <- matched   # fill tips

  # node_frac[i] stores the winning group's proportion at node i
  node_frac <- numeric(n_all)
  node_frac[seq_len(n_tips)] <- 1.0       # tips are 100 % their own group

  # Process internal nodes in post-order (leaves → root)
  # tree$edge[i, ] = c(parent, child)
  # We need children before parents → iterate edges in reverse
  for (i in rev(seq_len(nrow(tree$edge)))) {
    parent <- tree$edge[i, 1]
    child  <- tree$edge[i, 2]

    if (child <= n_tips) {
      # Tip: already set
      next
    }

    # Collect all children of this internal node
    child_idx <- tree$edge[tree$edge[, 1] == child, 2]
    if (length(child_idx) == 0) next

    child_groups <- node_group[child_idx]
    # Weight each child equally (one vote per direct child)
    tab <- table(child_groups)
    best_group <- names(tab)[which.max(tab)]
    best_frac  <- max(tab) / sum(tab)

    node_group[child] <- best_group
    node_frac[child]  <- best_frac
  }

  # Assign edge colors based on child node's group
  edge_colors <- vapply(seq_len(nrow(tree$edge)), function(i) {
    child <- tree$edge[i, 2]
    grp   <- node_group[child]
    frac  <- node_frac[child]

    if (!is.null(mixed_color) && frac < mixed_threshold) {
      return(mixed_color)
    }
    col <- palette[grp]
    if (is.na(col)) "#CCCCCC" else col
  }, character(1))

  edge_colors
}


#' Get the group of each tip's *direct parent* edge
#'
#' A simpler alternative to [branch_colors_majority()]: each edge is
#' colored by the *tip's own group* rather than the majority among
#' descendants. Only meaningful for terminal edges; internal edges
#' receive `mixed_color`.
#'
#' @inheritParams branch_colors_majority
#' @return Character vector of hex colors, length = `nrow(tree$edge)`.
#' @export
branch_colors_tip_only <- function(tree, tip_groups, palette,
                                   mixed_color = "#CCCCCC") {
  n_tips <- ape::Ntip(tree)
  matched <- tip_groups[tree$tip.label]

  vapply(seq_len(nrow(tree$edge)), function(i) {
    child <- tree$edge[i, 2]
    if (child <= n_tips) {
      grp <- matched[child]
      col <- palette[grp]
      if (is.na(col)) mixed_color else col
    } else {
      mixed_color
    }
  }, character(1))
}
