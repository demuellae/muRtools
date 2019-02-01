################################################################################
# data normalization
################################################################################
#' normalizeRank
#' 
#' Performs rank normalization on the columns of a matrix
#' @param X   A matrix which should be normalized
#' @param out output type. Either \code{"rank"} or \code{"percentile"}
#' @param ties.method method for breaking ties (see \code{?colRanks} for details)
#' @return a matrix containing the normalized values
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
		dimnames(percM) <- dimnames(percM)
		return(percM)
	} else {
		stop(paste0("Invalid 'out' parameter (normalizeRank): ", out))
	}
}

#' normalizePercentile
#' 
#' Performs percentile normalization on the columns of a matrix,
#' i.e. each element in a column will be the percentile it lies in in its column
#' @return a matrix containing the normalized values
#' @author Fabian Mueller
#' @export 
normalizePercentile <- function(X){
	# scale scores to their percentiles in columns
	res <- apply(X, 2, FUN=function(x){
		ecdf(x)(x)
	})
	dimnames(res) <- dimnames(X)
	return(res)
}
