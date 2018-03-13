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
#' @param ...    passed to \code{read.table}
#' @return the result of \code{read.table}
#' @author Fabian Mueller
#' @export
readTab <- function(fn, sep="\t", header=TRUE, stringsAsFactors=FALSE, ...){
	read.table(fn, header=header, sep=sep, comment.char="",  stringsAsFactors=stringsAsFactors, ...)
}
