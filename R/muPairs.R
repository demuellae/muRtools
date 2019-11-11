#' colorize.value
#' 
#' transforms a value to a color based on a colorscale and a value scale
#' @param val value to be transformed
#' @param rng value range
#' @param colscheme.col.val Color palette to pick a color from
#' @return color of the value
#' @export 
#' @examples 
#' colorize.value(0.5,rng=c(-1,1), colscheme.col.val=c(colorpanel(100,"blue","white"),colorpanel(100,"white","red")))
colorize.value <- function(val, rng=c(-1,1), colscheme.col.val=c(gplots::colorpanel(100,"blue","white"),gplots::colorpanel(100,"white","red"))){
	return(colscheme.col.val[round((val - rng[1]) / (rng[2] -rng[1])  * (length(colscheme.col.val)-1),0)+1])
}
#' panel.cor.col
#' 
#' creates a correlation panel. The correlation values are colored and fontsize is proportional to amount of correlation
#' @param x x
#' @param y y
#' @param digits number of digits to be displayed
#' @param prefix prefix for correlation text
#' @param cex.cor cex for correlation text
#' @param ... more plotting parameters
#' @return Nothing particularly interesting
panel.cor.col <- function(x, y, digits=2, prefix="", cex.cor, ...){
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	idx <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
	cc <- cor(x[idx], y[idx])
	cor.col <- colorize.value(cc,rng=c(-1,1),...)
	rect(0, 0, 1, 1, density = NULL, angle = 45,col = cor.col)
	r <- abs(cc)
	txt <- format(c(cc, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex.cor * r)
}
#' panel.density
#' 
#' creates a scatterplot density panel.
#' @param x x
#' @param y y
#' @param ... more plotting parameters
#' @return Nothing particularly interesting
panel.density <- function(x, y, ...){
	colorVals = densCols(x=x,y=y)
	points(x,y, col=colorVals, pch=20, cex=1)
}
#' pairsDensCor
#' 
#' Pair plot with correlation values colored and resized on the upper right diagonal and density scatter plot on the lower left diagonal
#' @param tt table to be visualized
#' @param colscheme colorscheme to be used for correlation
#' @param ... more plotting parameters
#' @return Nothing particularly interesting
#' @export 
#' @examples 
#' dd <- USJudgeRatings
#' dd[,5] <- 1/dd[,5]
#' tt <- dd[,1:5]
#' pairsDensCor(tt)
pairsDensCor <- function(tt,colscheme=c(gplots::colorpanel(100,"blue","white"),gplots::colorpanel(100,"white","red")),...){
	pairs(tt,lower.panel=panel.density,upper.panel=panel.cor.col,colscheme.col.val=colscheme,...)
}
