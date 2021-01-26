#' underscore2camel
#' 
#' converts underscores to camel case in a (vector of) string(s)
#' @param a string or string vector
#' @return the converted string(s)
#' @export 
#' @examples 
#' underscore2camel("bla_blubb")
underscore2camel <- function(x){
	x <- gsub("([[:alnum:]])_([0-9a-zA-Z])","\\1\\U\\2",x,perl=TRUE)
	return(x)
}

#' camel2underscore
#' 
#' converts camel case to underscores in a (vector of) string(s)
#' @param a string or string vector
#' @return the converted string(s)
#' @export 
#' @examples 
#' underscore2camel("bla_blubb")
camel2underscore <- function(x){
	x <- gsub("([a-z])([A-Z])","\\1_\\L\\2",x,perl=TRUE)
	return(x)
}

#' normalize.str
#' 
#' normalize a string by removing special characters and
#' replacing whitespaces and dots and afterwards
#' remove all leading and trailing whitespaces and special characters
#' By default, underscore is the replacement character. Avoids consecutive underscores.
#' @param a string or string vector
#' @param resolve.camel Is the string in camelCase and should this be resolved to snake_case?
#' @param return.camel should camelCase be outputed rather than undescores?
#' @return the normalized string(s)
#' @export 
#' @examples 
#' normalize.str("_ (b)lA BLu[bb.blA\tblubb- bla)\n_")
normalize.str <- function(x,resolve.camel=FALSE,return.camel=FALSE){
	if (resolve.camel){
		x <- camel2underscore(x)
	}
	x <- tolower(x)
	x <- gsub("[\\.[:space:]\\-]","_",x)
	x <- gsub("_+","_",x)
	x <- gsub("[^[:alnum:]_]","",x)
	x <- gsub("^[^[:alnum:]]+","",x)
	x <- gsub("[^[:alnum:]]+$","",x)
	if (return.camel){
		x <- underscore2camel(x)
	}
	return(x)
}

#' strTrim
#'
#' trim a character vector to have a desired length by taking the beginning of the string
#' and the end of the string
#'
#' @param ss       character vector
#' @param len.out  target output length to trim to
#' @param len.pref length of the prefix to be used from the original string
#' @param len.suf  length of the suffix to be used from the original string
#' @param trim.str string to place in between prefix and suffix
#' @return character vector in which each element has length\code{<=len.out}
#'
#' @author Fabian Mueller
#' @export
strTrim <- function(ss, len.out=50, trim.str="...", len.pref=ceiling((len.out-nchar(trim.str))/2), len.suf=len.out-len.pref-nchar(trim.str)){
	if (len.pref + len.suf + nchar(trim.str) > len.out){
		stop("Invalid string lengths specified")
	}
	lens <- nchar(ss)
	res <- ifelse(
		lens > len.out,
		paste0(substr(ss, 1, len.pref), trim.str, substr(ss, lens-len.suf+1, lens)),
		ss
	)
	return(res)
}

#' getHashString
#' 
#' Get a hash string, i.e. a string unlikely to occur again
#' @param pattern   a prefix that will be used in the returned hash string
#' @param useDate   Should the current time and date be used in the hash string to make it even more unique
#' @return a character string unlikely to occur again
#' @author Fabian Mueller
#' @export 
#' @examples 
#' getHashString()
getHashString <- function(pattern="", useDate=TRUE){
	pat <- pattern
	if (useDate) {
		pat <- format(Sys.time(), "%Y%m%d_%H%M%S_")
		if (nchar(pattern) > 0){
			pat <- paste(pattern, pat, sep="_")
		}
	}
	res <- basename(tempfile(pattern=pat))
	return(res)
}
