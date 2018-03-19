################################################################################
# Utilities and wrappers for plotting heatmaps
################################################################################
#' plotCorPhm
#' 
#' Plot a correlation matrix as heatmap. Wraps around \code{pheatmap}. Note that no clustering will be performed if not supplied with an appropriate clustering dendrogram
#' @param cc        a correlation matrix (as returned by \code{cor()})
#' @param clustDend a clustering dendrogram to be used. Set to \code{NULL} to disable clustering dendrogram
#' @param sampleAnnot a data.frame containing sample information to color by (corrsponds to \code{annotation_row} and \code{annotation_col} parameters of \code{pheatmap()})
#' @param color     see \code{?pheatmap} for details
#' @param breaks    see \code{?pheatmap} for details. In this wrapper, the default value corresponds to splitting color across the full range of correlation coefficients [-1,1]
#' @param border_color     see \code{?pheatmap} for details
#' @param ...       parameters passed on to \code{pheatmap()}
#' @return invisibly the result of a call to \code{pheatmap()}
#' @author Fabian Mueller
#' @export 
plotCorPhm <- function(
		cc,
		clustDend=NULL,
		sampleAnnot=NA,
		color=colorRampPalette(rev(brewer.pal(n=11, name = "RdBu")))(100),
		breaks=seq(-1-1e-6, 1+1e-6, length.out=length(color)),
		border_color=NA,
		...){
	require(pheatmap)
	# color.base <- c("#EDF8B1","#41B6C4","#081D58")
	require(RColorBrewer)
	color.base <- rev(brewer.pal(11,"RdBu")[c(1,6,11)])
	color.panel <- colorpanel(100,color.base[1],color.base[2],color.base[3])
	clustr <- FALSE
	if (!is.null(clustDend)){
		if (!is.element("dendrogram", class(clustDend))){
			stop("clustDend must be a dendrogram")
		}
		clustr <- as.hclust(clustDend)
	}

	pheatmap(cc, color=color, breaks=breaks, border_color=border_color, annotation_row=sampleAnnot, annotation_col=sampleAnnot, cluster_rows=clustr, cluster_cols=clustr, ...)
}
