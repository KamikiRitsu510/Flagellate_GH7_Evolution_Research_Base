#' Generate a color palette for phylogenetic groups
#'
#' Assigns a distinct color to each unique group label. Colors are chosen
#' from a curated set that is colorblind-friendly and print-safe, with
#' `"Other"` always rendered in light gray.
#'
#' @param groups Character vector of group labels (as returned by
#'   [classify_tips()]).
#' @param user_colors Optional named character vector of hex colors that
#'   *override* automatic assignments for specific groups. Any group not
#'   listed here gets an auto-assigned color.
#' @param other_color Hex color used for the `"Other"` group (and any group
#'   whose name matches `other_label`). Default: `"#AAAAAA"`.
#' @param other_label The label used for unclassified tips. Default:
#'   `"Other"`.
#' @return A named character vector of hex colors, one per unique group.
#' @examples
#' groups <- c("Flagellate", "Flagellate", "Fungi", "Bacteria", "Other")
#' auto_palette(groups)
#' @export
auto_palette <- function(groups,
                         user_colors = NULL,
                         other_color = "#AAAAAA",
                         other_label = "Other") {
  # Curated, colorblind-friendly palette (Wong 2011 + extensions)
  base_colors <- c(
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#4DAF4A",  # green
    "#984EA3",  # purple
    "#FF7F00",  # orange
    "#A65628",  # brown
    "#F781BF",  # pink
    "#1B9E77",  # teal
    "#D95F02",  # dark orange
    "#7570B3",  # lavender
    "#E7298A",  # magenta
    "#66A61E"   # olive green
  )

  unique_groups <- setdiff(unique(groups), other_label)
  n <- length(unique_groups)

  if (n > length(base_colors)) {
    # Fall back to rainbow if we have more groups than base colors
    extra <- grDevices::rainbow(n - length(base_colors), s = 0.8, v = 0.85)
    base_colors <- c(base_colors, extra)
  }

  palette <- setNames(base_colors[seq_len(n)], unique_groups)
  palette[other_label] <- other_color

  # Apply user overrides
  if (!is.null(user_colors)) {
    for (nm in names(user_colors)) {
      palette[nm] <- user_colors[nm]
    }
  }

  palette
}
