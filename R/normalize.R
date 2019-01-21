################################################################################
# data normalization
################################################################################
#' normalizeRank
#' 
#' Performs rank normalization on the columns of a matrix
#' @param X   A matrix which should be normalized
#' @param out output type. Either \code{"rank"} or \code{"percentile"}
#' @param ties.method method for breaking ties (see \code{?colRanks} for details)
#' @return a matrix containing the normalized avlues
#' @author Fabian Mueller
#' @export 
normalizeRank <- function(X, out="percentile", ties.method="average"){
	require(matrixStats)
	rankM <- colRanks(X, ties.method=ties.method, preserveShape=TRUE)
	dimnames(rankM) <- dimnames(X)
	if (out=="rank"){
		return(rankM)
	} else if (out=="percentile"){
		percM <- t(t(rankM)/colMaxs(rankM, na.rm=TRUE))
		return(percM)
	} else {
		stop(paste0("Invalid 'out' parameter (normalizeRank): ", out))
	}
}
