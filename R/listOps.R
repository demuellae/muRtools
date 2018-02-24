#' indicesInList
#' 
#' Find the occurrences of items in a list of vectors
#' @param x vector of items to be found in the list
#' @param l list of vectors in which \code{x} should be found
#' @return a list containing an element for each item in \code{x} that contains the indices of its occurrence in \code{l}
#' @export 
#' @examples
#' l <- list(1:3, 4:5, 5:9)
#' x <- c(2,3,5,666,8,5)
#' indicesInList(x, l)
indicesInList <- function(x, l){
	# adapted from http://stackoverflow.com/questions/11002391/fast-way-of-getting-index-of-match-in-list
	g <- rep(seq_along(l), sapply(l, length))
	ll <- unlist(l)
	lu <- unique(ll)
	lf <- factor(ll, levels=lu)
	gg <- split(g, lf)
	gg[match(x, lu)]
}

#' combinationList
#' 
#' get a list of all combinations of vectors. Basically a wrapper around \code{\link{expand.grid}}
#' @param ...  vectors of elements. Ideally named
#' @return a list containing all combinations of elements in the input. Each element contains a unique combination
#' @export 
#' @examples
#' combinationList(a=letters[1:5], A=LETTERS[1:3], i=1:4) 
combinationList <- function(...){
	dd <- expand.grid(..., stringsAsFactors=FALSE)
	attr(dd, "out.attrs") <- NULL
	res <- lapply(1:nrow(dd), FUN=function(i){as.list(dd[i,])})
	return(res)
}
