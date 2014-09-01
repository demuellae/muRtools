HISTONE.MOD.PATTERN <- "(h[[:digit:]]+)([k])([[:digit:]]+)(me[[:digit:]]?|ac|ub)"

#' Methods for recognizing histone modifications from strings
#'
#' @param s a string
#' @param s1 a string
#' @param s2 a string
#' \describe{
#'   \item{\code{containsHistoneModStr}}{
#'   	Does a string contain the pattern for histone modifications?
#'   }
#'   \item{\code{getHistoneFromHistoneModStr}}{
#'   	Retrive the histone from a string containing a histone modification pattern
#'   }
#'   \item{\code{getAaTypeFromHistoneModStr}}{
#'   	Retrive the amino acid type from a string containing a histone modification pattern
#'   }
#'   \item{\code{getAaPosFromHistoneModStr}}{
#'   	Retrive the amino acid position from a string containing a histone modification pattern
#'   }
#'   \item{\code{getModFromHistoneModStr}}{
#'   	Retrive the modification from a string containing a histone modification pattern
#'   }
#'   \item{\code{normalizeHistoneModStr}}{
#'   	normalize the histone modification if contained in a string
#'   }
#'   \item{\code{matchHistoneModStr}}{
#'   	Do two strings contain the same histone modifications
#'   }
#' }
#'	
#' @export
#' @rdname HistoneModStr
#' @aliases containsHistoneModStr,getHistoneFromHistoneModStr,getAaPosFromHistoneModStr,getModFromHistoneModStr,matchHistoneModStr
#' @examples 
#' s1 <- "H3K4me3"
#' s2 <- "h3k04ME3"
#' s3 <- "blubb5A27me3"
#' s4 <- "h3k27ac"
#' containsHistoneModStr(s1)
#' containsHistoneModStr(s3)
#' containsHistoneModStr(s4)
#' getHistoneFromHistoneModStr(s1)
#' getAaTypeFromHistoneModStr(s1)
#' getAaPosFromHistoneModStr(s1)
#' getModFromHistoneModStr(s1)
#' normalizeHistoneModStr(s1)
#' normalizeHistoneModStr(c(s1,s2,s3,s4))
#' matchHistoneModStr(s1,s2)
#' matchHistoneModStr(s1,s4)
containsHistoneModStr <- function(s){
	grepl(HISTONE.MOD.PATTERN,s,ignore.case=TRUE)
}
#' @rdname HistoneModStr
#' @export
getHistoneFromHistoneModStr <- function(s){
	ifelse(containsHistoneModStr(s),
		gsub(paste0(".*",HISTONE.MOD.PATTERN,".*"),"\\1",s,ignore.case=TRUE),
		NA
	)
}
#' @rdname HistoneModStr
#' @export
getAaTypeFromHistoneModStr <- function(s){
	ifelse(containsHistoneModStr(s),
		gsub(paste0(".*",HISTONE.MOD.PATTERN,".*"),"\\2",s,ignore.case=TRUE),
		NA
	)
}
#' @rdname HistoneModStr
#' @export
getAaPosFromHistoneModStr <- function(s){
	ifelse(containsHistoneModStr(s),
		as.integer(gsub(paste0(".*",HISTONE.MOD.PATTERN,".*"),"\\3",s,ignore.case=TRUE)),
		NA
	)
}
#' @rdname HistoneModStr
#' @export
getModFromHistoneModStr <- function(s){
	ifelse(containsHistoneModStr(s),
		gsub(paste0(".*",HISTONE.MOD.PATTERN,".*"),"\\4",s,ignore.case=TRUE),
		NA
	)
}
#' @rdname HistoneModStr
#' @export
normalizeHistoneModStr <- function(s){
	ifelse(containsHistoneModStr(s),
		tolower(paste0(getHistoneFromHistoneModStr(s),getAaTypeFromHistoneModStr(s),getAaPosFromHistoneModStr(s),getModFromHistoneModStr(s))),
		s
	)
}
#' @rdname HistoneModStr
#' @export
matchHistoneModStr <- function(s1,s2){
	if (!all(containsHistoneModStr(c(s1,s2)))){
		return(NA)
	}
	res <- (tolower(getHistoneFromHistoneModStr(s1)) == tolower(getHistoneFromHistoneModStr(s2))) &
		   (getAaPosFromHistoneModStr(s1) == getAaPosFromHistoneModStr(s2)) &
		   (tolower(getModFromHistoneModStr(s1)) == tolower(getModFromHistoneModStr(s2)))
	return(res)
}
