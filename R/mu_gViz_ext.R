#' gviz_bwScheme
#' 
#' A Gviz scheme with black labels and transparent boxes for labels
#' @return (invisibly) the scheme structure
#' @export
gviz_bwScheme <- function() {
  require(Gviz)
  scheme <- getScheme()
  scheme$GdObject$fontface.title <- 1
  scheme$GdObject$background.title <- "transparent"
  scheme$GdObject$col.border.title <- "black"
  scheme$GdObject$fontcolor.title <- "black"
  scheme$GdObject$col.title <- "black"
  scheme$GdObject$fontcolor <- "black"
  scheme$IdeogramTrack$fontcolor <- "black"
  scheme$GenomeAxisTrack$fontcolor <- "black"
  scheme$GenomeAxisTrack$col <- "black"
  scheme$GenomeAxisTrack$labelPos <- "above"
  scheme$AnnotationTrack$col <- NULL
  scheme$AnnotationTrack$fontface.group <- 1
  scheme$AnnotationTrack$fontcolor.group <- "black"
  scheme$GeneRegionTrack$col <- NULL
  scheme$GeneRegionTrack$col.line <- "black"
  addScheme(scheme, "bwScheme")
  options(Gviz.scheme="bwScheme")
  invisible(scheme)
}
