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
	ddd <- data.frame(mm, check.names=FALSE)
	ddd$rnames <- rownames(mm)
	dd <- melt(ddd)
	dd$rnames <- factor(dd$rnames,levels=rownames(mm))
	p <-   ggplot(dd, aes(variable,rnames)) + geom_tile(aes(fill = value)) + 
		   scale_x_discrete(name="") + scale_y_discrete(limits = rev(levels(dd$rnames)),name="") + coord_fixed() +
       theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
	if(add.text){
		p <- p + geom_text(aes(label=value))
	}
	return(p)
}

#' getPointDensity
#' 
#' Get point density of points in 2 dimensions. Code from http://slowkow.com/notes/ggplot2-color-by-density/
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @return The density within each square.
getPointDensity <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
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
colpal.bde <- c("#2166AC","#B2182B","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#00441B","#40004B","#053061")
#' \describe{
#'   \item{\code{colpal.nature}}{
#'         Color palette inspired by Nature journal color scheme's.
#'   }
#' }
#' @rdname colpal
#' @export
#' @examples
#' plot.new()
#' gradient.rect(0,0,1,1,col=colpal.nature,nslices=length(colpal.nature),gradient="x",border=NA)
colpal.nature <- c("#003D7C", "#D50911", "#0086A8", "#008136", "#7C68A4", "#8E1A47", "#E67800", "#709F28", "#008FB4", "#84486A", "#B5797F", "#7489A8",  "#6C9396", "#7D9FB1", "#84486A", "#7C698B", "#88A2C3")
#' \describe{
#'   \item{\code{colpal.nature}}{
#'         Color palette inspired by Nature journal color scheme's.
#'   }
#' }
#' @rdname colpal
#' @export
#' @examples
#' plot.new()
#' gradient.rect(0,0,1,1,col=colpal.mu.cat,nslices=length(colpal.nature),gradient="x",border=NA)
colpal.mu.cat <- c("#e69f00", "#56b4e9", "#74c476", "#cc79a7", "#d55e00", "#0072b2", "#009e73", "#6a51a3", "#f0e442", "#999999", "#000000")

#' theme_nogrid
#' 
#' A ggplot2 theme based on theme_bw but with no grid lines and axis only on top and bottom
#' @param base_size base size
#' @param base_family base family
#' @return the theme structure
#' @export
#' @aliases theme_nogrid
#' @examples 
#' theme_set(theme_nogrid())
#' dframe <- data.frame(x=runif(100),y=runif(100))
#' ggplot(dframe,aes(x=x,y=y)) + geom_point()
theme_nogrid <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line         = element_line(colour = "black",lineend="round"),
      axis.text         = element_text(size = rel(0.8), colour = "black"),
      axis.ticks        = element_line(colour = "black"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      panel.background  = element_blank(),
      panel.border      = element_blank(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey80", colour = "grey50", size = 0.2)
    )
}

#' ggsave4doc
#' 
#' Wrapper around ggsave that has default values for parameters fitting for embedding
#' plots into my documents
#' @param fn      file name
#' @param plot    see \code{?ggsave}
#' @param width   see \code{?ggsave}
#' @param height  see \code{?ggsave}
#' @param units   see \code{?ggsave}
#' @param family  see \code{?ggsave}
#' @param useDingbats see \code{?ggsave}
#' @param ...     see \code{?ggsave}
#' @return result of \code{ggsave} command
#' @export
#' @aliases ggsave4doc
ggsave4doc <- function(fn, plot=last_plot(), width=146, height=146, units="mm", family="Palatino", useDingbats=FALSE,...){
  ggsave(fn, plot=plot, width=width, height=height, units=units, family=family, useDingbats=useDingbats,...)
}

#' ggtemp
#' 
#' Wrapper for quickly saving plot to temporary file
#' @param plot    see \code{?ggsave}
#' @param fn      file name
#' @param ...     see \code{?ggsave}
#' @return result of \code{ggsave} command
#' @export
#' @aliases ggtemp
ggtemp <- function(plot=last_plot(), fn=paste0("~/tmp_work/", getHashString("ggplot"), ".pdf"), ...){
  ggsave4doc(fn, plot=plot, ...)
}

#' pdftemp
#' 
#' Wrapper for quickly saving plot to temporary pdf file. terminate using \code{dev.off()}
#' @param fn      file name
#' @param ...     see \code{?pdf}
#' @return nothing of particular interest
#' @export
#' @aliases pdftemp
pdftemp <- function(fn=paste0("~/tmp_work/", getHashString("rplot"), ".pdf"), ...){
  pdf(fn, ...)
}
#' pngtemp
#' 
#' Wrapper for quickly saving plot to temporary png file. terminate using \code{dev.off()}
#' @param fn      file name
#' @param ...     see \code{?png}
#' @return nothing of particular interest
#' @export
#' @aliases pngtemp
pngtemp <- function(fn=paste0("~/tmp_work/", getHashString("rplot"), ".png"), width=1024, height=1024, ...){
  png(fn, width=width, height=height, ...)
}
