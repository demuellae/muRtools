#' ggplot2.heatmap
#' 
#' converts a matrix or dataframe into a ggplot2 object for subsequent plotting.
#' @param mm matrix or dataframe to be plottet as heatmap
#' @param add.text logical. should the cells be labelled with the values
#' @return the ggplot2 object (can be extended for plotting)
#' @export
#' @aliases ggplot2.heatmap
#' @examples 
#' ggplot2.heatmap(airquality[1:15,],add.text=TRUE) + scale_fill_gradient(low = "red",high = "steelblue"))
#ggplot2.heatmap <- function(mm,add.text=FALSE){
#	ddd <- data.frame(mm)
#	ddd$rnames <- rownames(mm)
#	dd <- melt(ddd)
#	dd$rnames <- factor(dd$rnames,levels=rev(rownames(mm)))
#	p <-   ggplot(dd, aes(variable,rnames)) + geom_tile(aes(fill = value)) + scale_x_discrete(name="") + scale_y_discrete(name="")
#	if(add.text){
#		p <- p + geom_text(aes(label=value))
#	}
#	return(p)
#}
ggplot2.heatmap <- function(mm,add.text=FALSE){
	require(reshape)
	ddd <- data.frame(mm)
	ddd$rnames <- rownames(mm)
	dd <- melt(ddd)
	dd$rnames <- factor(dd$rnames,levels=rownames(mm))
	p <-   ggplot(dd, aes(variable,rnames)) + geom_tile(aes(fill = value)) + 
		   scale_x_discrete(name="") + scale_y_discrete(limits = rev(levels(dd$rnames)),name="") + coord_fixed()
	if(add.text){
		p <- p + geom_text(aes(label=value))
	}
	return(p)
}

#' Custom Color Paletes
#' 
#' \describe{
#'   \item{\code{colpal.cb}}{
#'         color blind friendly color palettes (adapted from http://wiki.stdout.org/rcookbook/Graphs/Colors%20%28ggplot2%29/)
#'   }
#' }
#'	
#' @export
#' @rdname colpal
#' @aliases colpal,colpal.cb,colpal.colpal.bde
#' @examples 
#' plot.new()
#' gradient.rect(0,0,1,1,col=colpal.cb,nslices=length(colpal.cb),gradient="x",border=NA)
colpal.cb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000")
#' \describe{
#'   \item{\code{colpal.bde}}{
#'         Enhanced Color Brewer palette 'Dark2'
#'   }
#' }
#' @rdname colpal
#' @export
#' @examples
#' plot.new()
#' gradient.rect(0,0,1,1,col=colpal.bde,nslices=length(colpal.bde),gradient="x",border=NA)
colpal.bde <- c("#2166AC","#B2182B","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")
