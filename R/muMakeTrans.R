#' col.text.2.hex
#' 
#' convert a color from a text string to a hex string (only works with scalars, not on vectors)
#' @param ss color string to be transformed
#' @param alpha maximum alpha
#' @return a string containing the hex code of the color
#' @export 
#' @examples 
#' col.text.2.hex("dark blue")
col.text.2.hex <- function(ss,alpha=255){
	c2r <- col2rgb(ss)
	return(rgb(c2r[1,1],c2r[2,1],c2r[3,1],alpha=alpha,maxColorValue=255))
}
#' makeTrans
#' 
#' make a vector of colors transparent
#' @param ccc color vector to make transparent
#' @param isText are the colors textstrings (as opposed to hex value strings)?
#' @param transparancy.val string value for transparancy
#' @return a vector of colors with transparancy (hex strings)
#' @export 
#' @examples 
#' makeTrans(rainbow(6))
#' makeTrans(c("red","dark blue","coral"),isText=TRUE)
makeTrans <- function(ccc,isText=FALSE,transparancy.val="44"){
	if (isText) ccc <- sapply(ccc,FUN=function(cc){col.text.2.hex(cc)})
	return(sapply(ccc,FUN=function(cc){paste(substr(cc,1,7),transparancy.val,sep="")}))
}
