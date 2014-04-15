#' scatter.twogroups
#' 
#' density scatterplot with highlighting the points of a certain group in a different color
#' @param x x coordinates
#' @param y y coordinates
#' @param is.g logical vector indicating whether a point is in the highlighted group
#' @param cols.all colorscheme for all points
#' @param cols.g colorscheme for the points in the highlighted group
#' @param x x coordinates
#' @param ... more plotting parameters
#' @return Nothing particularly interesting
#' @export
scatter.twogroups <- function(x,y,is.g,cols.all=blues9[-(1:3)],cols.g=c("coral","darkred"),...){
	colorVals.all = densCols(x=x,y=y,colramp = colorRampPalette(cols.all))
	if (any(is.g,na.rm=TRUE)) colorVals.g = densCols(x=x[is.g],y=y[is.g],colramp = colorRampPalette(cols.g))
	plot(x,y, col=colorVals.all, pch=20, cex=1,...)
	if (any(is.g,na.rm=TRUE)) points(x[is.g],y[is.g], col=colorVals.g, pch=20, cex=1,...)
}
