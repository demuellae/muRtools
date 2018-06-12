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
#' @param ...    passed to \code{read.table}
#' @return the result of \code{read.table}
#' @author Fabian Mueller
#' @export
readTab <- function(fn, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", ...){
	read.table(fn, header=header, sep=sep, comment.char="",  stringsAsFactors=stringsAsFactors, quote=quote, ...)
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
writeTab <- function(fn, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, ...){
	write.table(fn, sep=sep, row.names=row.names, col.names=col.names, quote=quote, ...)
}
