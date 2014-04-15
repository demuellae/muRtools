#' kde.plot
#' 
#' plot a single vector as kernel density estimation (kde.plot.simple), plot columns of a matrix as kernel density estimation (kde.plot.matrix)
#' @param x values
#' @param X matrix
#' @param initial is the plot initial, i.e. should the plot() function be called
#' @param col color(s)
#' @param legend add a color legend with colnames of the matrix
#' @param ... more graphical parameters
#' @return Nothing particularly interesting
#' @export
#' @aliases kde.plot
#' @aliases kde.plot.simple
#' @aliases kde.plot.matrix
#' @rdname kde.plot
#' @examples 
#' kde.plot.simple(c(1:10,rep(5,3)),initial=TRUE,col=c("#00640044"),main="Some plot")
#' kde.plot.simple(rep(8,2),initial=FALSE,col=c("#64640044"))
#' 
#' dd <- USJudgeRatings
#' dd[,5] <- 1/dd[,5]
#' tt <- dd[,1:5]
#' kde.plot.matrix(tt)
kde.plot.simple <- function(x,initial=FALSE,col=c("#00640044"),...){
	#Kernel Density Estimation
	dens=density(x,na.rm=TRUE)
	dens$x <- c(min(dens$x),dens$x,max(dens$x)) #add end points to avoid weird corners
	dens$y <- c(0,dens$y,0)
	if(initial) plot(dens,frame.plot=FALSE,cex.axis=1,cex.main=1,...)
	polygon(dens, col=col)
}
#' @rdname kde.plot
#' @export
kde.plot.matrix <- function(X,col=makeTrans(rainbow(ncol(X))),legend=TRUE,...){
	kde.plot.simple(X[,1],col=col[1],initial=TRUE,...)
	for (i in 2:ncol(X)){
		kde.plot.simple(X[,i],col=col[i],initial=FALSE,...)
	}
	if (legend){
		legend('topright',colnames(X),fill=col)
	}
}
