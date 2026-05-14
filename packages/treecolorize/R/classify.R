#' Strip accession prefix from a tip label
#'
#' Removes the first whitespace-delimited token from a tip label,
#' returning the remaining string as the "species name".
#'
#' @param labels Character vector of tip labels.
#' @param sep Separator between accession and species name. Default is any
#'   whitespace (`" "`). Use `"|"` for pipe-delimited NCBI-style labels.
#' @return Character vector of the same length with prefixes removed.
#'   Labels that contain no separator are returned unchanged.
#' @examples
#' strip_accession(c("AAY83390.3 Pseudotrichonympha grassii",
#'                   "XP_001212905.1 Chaetomium globosum"))
#' # [1] "Pseudotrichonympha grassii" "Chaetomium globosum"
#' @export
strip_accession <- function(labels, sep = " ") {
  sapply(labels, function(lbl) {
    parts <- strsplit(lbl, sep, fixed = TRUE)[[1]]
    if (length(parts) > 1) paste(parts[-1], collapse = sep) else lbl
  }, USE.NAMES = FALSE)
}


#' Classify phylogenetic tip labels by keyword matching
#'
#' Maps each tip label to a group name by matching against a user-supplied
#' list of regular expressions. Matching is case-insensitive and stops at
#' the first hit, so order the list from most specific to least specific.
#'
#' @param tips Character vector of tip labels (e.g. `tree$tip.label`).
#' @param keyword_map A named list where each *name* is a regular expression
#'   and each *value* is the group label to assign when that regex matches.
#'   Example:
#'   \code{list("pseudotrichonympha|trichonympha" = "Flagellate",
#'              "bacillus|clostridium" = "Bacteria")}
#' @param other_label Label assigned to tips that match no pattern.
#'   Default: `"Other"`.
#' @param strip_prefix If `TRUE` (default), the first whitespace-delimited
#'   token is stripped before matching, so that NCBI-style labels like
#'   `"XP_001212905.1 Chaetomium globosum"` are classified on the species
#'   name rather than the accession.
#' @param sep Separator used by [strip_accession()] when `strip_prefix = TRUE`.
#' @return A named character vector of group labels, with names equal to
#'   the original `tips`.
#' @seealso [strip_accession()], [auto_palette()]
#' @examples
#' tips <- c(
#'   "BAB64565.3 Holomastigotoides mirabile",
#'   "XP_006965674.1 Chaetomium globosum",
#'   "WP_165634028.1 Cellulomonas sp."
#' )
#' kmap <- list(
#'   "holomastigotoides|trichonympha|pseudotrichonympha" = "Flagellate",
#'   "chaetomium|aspergillus|fungi"                      = "Fungi",
#'   "cellulomonas|bacillus|clostridium"                 = "Bacteria"
#' )
#' classify_tips(tips, kmap)
#' @export
classify_tips <- function(tips,
                          keyword_map,
                          other_label  = "Other",
                          strip_prefix = TRUE,
                          sep          = " ") {
  if (!is.list(keyword_map) || is.null(names(keyword_map))) {
    stop("`keyword_map` must be a named list.")
  }

  targets <- if (strip_prefix) strip_accession(tips, sep = sep) else tips

  groups <- vapply(targets, function(txt) {
    txt_lower <- tolower(txt)
    for (pattern in names(keyword_map)) {
      if (grepl(pattern, txt_lower, perl = TRUE)) {
        return(keyword_map[[pattern]])
      }
    }
    other_label
  }, character(1))

  names(groups) <- tips
  groups
}


#' Default keyword map for common phylogenetic groups
#'
#' A ready-to-use keyword map covering the major kingdoms and several
#' ecologically common lineages. Pass it directly to [classify_tips()] or
#' use it as a starting point to build your own.
#'
#' @return A named list suitable for use as `keyword_map` in [classify_tips()].
#' @examples
#' kmap <- default_keyword_map()
#' classify_tips(tree$tip.label, kmap)
#' @export
default_keyword_map <- function() {
  list(
    # Flagellate gut symbionts
    paste(c("pseudotrichonympha", "holomastigotoides", "trichonympha",
            "spirotrichonympha", "oxymonad", "flagellate"), collapse = "|") = "Flagellate",
    # Bacteria
    paste(c("bacillus", "clostridium", "cellulomonas", "streptomyces",
            "saccharothrix", "bacterium", "proteobacteria",
            "actinobacteria"), collapse = "|") = "Bacteria",
    # Fungi / Oomycetes
    paste(c("chaetomium", "phytophthora", "pyrenophora", "aspergillus",
            "trichoderma", "neurospora", "penicillium",
            "ascomycota", "basidiomycota", "fungi", "fungus",
            "mucor"), collapse = "|") = "Fungi",
    # Crustacea
    "daphnia|pulex|water flea|crustacea" = "Crustacea",
    # Plants
    "arabidopsis|populus|oryza|zea mays|plant|viridiplantae" = "Plant",
    # Insects
    "drosophila|apis|bombyx|insecta" = "Insect"
  )
}
