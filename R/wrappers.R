################################################################################
# Wrappers, convenience functions and shortcuts for existing functions
################################################################################
#' readTab
#' 
#' Wrapper around \code{read.table} to read tab-separated tables by default
#' @param fn	 filename to read
#' @param sep	 see \code{?read.table}
#' @param header see \code{?read.table}
#' @param stringsAsFactors see \code{?read.table}
#' @param quote  see \code{?read.table}
#' @param comment.char see \code{?read.table}
#' @param ...    passed to \code{read.table}
#' @return the result of \code{read.table}
#' @author Fabian Mueller
#' @export
readTab <- function(fn, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", comment.char="",  na.strings="", ...){
	read.table(fn, header=header, sep=sep, comment.char=comment.char,  stringsAsFactors=stringsAsFactors, quote=quote, na.strings=na.strings, ...)
}

#' writeTab
#' 
#' Wrapper around \code{write.table} to write tab-separated tables by default
#' @param x      table to write to file
#' @param fn	 filename to write to
#' @param sep	 see \code{?write.table}
#' @param row.names	 see \code{?write.table}
#' @param col.names	 see \code{?write.table}
#' @param quote  see \code{?write.table}
#' @param stringsAsFactors see \code{?write.table}
#' @param ...    passed to \code{write.table}
#' @return the result of \code{write.table}
#' @author Fabian Mueller
#' @export
writeTab <- function(x, fn, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, ...){
	write.table(x, fn, sep=sep, row.names=row.names, col.names=col.names, quote=quote, ...)
}

#' aggregateDf
#' 
#' Wrapper around \code{aggregate} to merge rows of a data frame based on a specified grouping.
#' Merging is done by either returning the unique value for each group or if multiple different
#' values exists converting them to character and concatenating by ";"
#' @param df     \code{data.frame} to aggregate
#' @param groupBy vector whose unique values indicates the grouping of the rows in the data frame
#' @return aggregated \code{data.frame}
#' @author Fabian Mueller
#' @export
aggregateDf <- function(df, groupBy){
	if (length(groupBy) != nrow(df)) stop("invalid value for 'groupBy'")
	res <- data.frame(
		lapply(
			aggregate(df, list(groupBy), FUN=function(x){
				if (length(unique(x))==1) return(x[1])
				return(paste(sort(as.character(x)), collapse=";"))
			}, simplify=FALSE),
			FUN=function(x){
				hasFactor <- any(sapply(x, is.factor))
				hasChar   <- any(sapply(x, is.character))
				if (hasFactor && hasChar) x <- lapply(x, as.character)
				return(unlist(x))
			}
		),
		stringsAsFactors=FALSE
	)
	grpNames <- as.character(res[,1]) # first column "Group.1"
	res <- res[,-1] #remove the first column "Group.1"
	rownames(res) <- grpNames
	return(res)
}
