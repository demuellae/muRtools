################################################################################
# Utilities for clustering and correlation
################################################################################
#' dist.correlation
#' 
#' Compute a distance matrix based on 1-corelation
#' @param x     a matrix on which the distances should be computed
#' @param ...   parameters passed on to \code{cor()}
#' @return a distance matrix
#' @author Fabian Mueller
#' @export 
dist.correlation <- function(x, ...) {
	res <- as.dist(1 - cor(t(x), use="pairwise.complete.obs", ...))
	res[is.na(res)] <- 1
	return(res)
}
#' getClusteringDendrogram
#' 
#' Get a clustering dendrogram using hierarchical clustering (wrapper)
#' @param X     A matrix for which the sample clustering dendrogram should be computed. Samples correspond to columns and features correspond to rows.
#' @param samplesOrdered   character vector specifying the preferred order of samples
#' @param distMethod distance metric to be used for clusteing. must be either "cor" or a valid distance method for \code{dist()}
#' @param linkMethod linkage method (see \code{hclust} for details)
#' @param corMethod  method for computing correlation coefficients. Only relevant if \code{distMethod=="cor"}.
#' @return clustering dendrogram (\code{dendrogram} object)
#' @author Fabian Mueller
#' @export 
getClusteringDendrogram <- function(X, samplesOrdered=colnames(X), distMethod="cor", linkMethod="ward.D", corMethod="pearson"){
	dd <- NULL
	if (distMethod=="cor"){
		dd <- dist.correlation(t(X), method=corMethod)
	} else {
		dd <- dist(t(X), method=distMethod)
	}
	noDist <- is.na(dd)
	if (sum(noDist) > 0) {
		logger.warning(c("Some distances could not be computed (",sum(noDist),") --> replacing with max distance"))
		dd[noDist] <- max(dd, na.rm=TRUE) #hack: set to max distance
	}
	clustRes <- hclust(dd, method=linkMethod)
	clustDend <- as.dendrogram(clustRes)

	#order the samples in the dendrogram as much as possible
	desired.order.wts <- 1:length(samplesOrdered)
	names(desired.order.wts) <- samplesOrdered 
	wwts <- desired.order.wts[colnames(X)]
	clustDend <- reorder(clustDend, wwts, agglo.FUN=mean)
	return(clustDend)
}
