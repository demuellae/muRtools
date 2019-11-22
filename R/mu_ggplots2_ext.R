#' @include stringOps.R
NULL

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

#' ggplot2.distr
#' 
#' Distribution plot combining violin and boxplot using ggplot2
#' @param x vector of values whose distribution is to be plottet
#' @param fillColor color to be used to fill the violin
#' @return the ggplot2 object (can be extended for plotting)
#' @export
#' @examples 
#' x <- rnorm(1000)
#' ggplot2.distr(x)
#}
ggplot2.distr <- function(x, fillColor="#676D8D"){
  if (!is.data.frame(x)){
    x <- data.frame(value=x)
  }
  x[,"group"] <- ".all"
  valCol <- colnames(x)[1]
  pp <- ggplot(x) + aes_string(x="group", y=valCol) +
        geom_violin(adjust=1, fill=fillColor) +
        geom_boxplot(aes(fill=NULL), outlier.shape=NA, width=0.2) +
        coord_flip() +
        # theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
        theme(axis.title.y=element_blank()) + guides(fill=FALSE) 
  return(pp)
}

#' ggMessagePlot
#'
#' Creates a plot, using \pkg{ggplot2}, with a single text message.
#'
#' @param txt Text to be plotted.
#' @return The newly initialized \code{ggplot} instance.
#'
#' @examples
#' \donttest{
#' ggMessagePlot("Missing data")
#' }
#' @export
ggMessagePlot <- function(txt) {
  if (!(is.character(txt) && length(txt) == 1 && (!is.na(txt)))) {
    stop("invalid value for txt")
  }
  ggplot(data.frame(x = 1, y = 1, labeltext = txt), aes_string("x", "y", label = "labeltext")) +
    geom_text(color = "grey50") +
    theme(axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
      panel.background = element_blank(), plot.background = element_blank())
}

#' getPointDensity
#' 
#' Get point density of points in 2 dimensions. Code from http://slowkow.com/notes/ggplot2-color-by-density/
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @return The density within each square
#' @export
getPointDensity <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

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
theme_nogrid <- function(base_size = 8, base_family = "Helvetica") {
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
      strip.background  = element_rect(fill = "white", colour = "black", size = 0.2)
    )
}

#' ggAutoColorScale
#' 
#' Automatical color scales for values for ggplots
#' @param x vector of values
#' @param method method for scaling: \code{"color"} or \code{"fill"}
#' @param symmetric treat numeric values as symmetric. If there are values smaller and larger than 0, a diverging color scheme will be applied
#' @return the theme structure
#' @export
#' @examples 
#' dframe.num.pos <- data.frame(x=runif(100),y=runif(100))
#' ggplot(dframe.num.pos, aes(x=x,y=y, color=x)) + geom_point() + ggAutoColorScale(dframe.num.pos[,"x"])
#' dframe.num.sym <- data.frame(x=rnorm(100),y=rnorm(100))
#' ggplot(dframe.num.sym, aes(x=x,y=y, color=x)) + geom_point() + ggAutoColorScale(dframe.num.sym[,"x"])
#' dframe.num.sym.lab <- data.frame(x=rnorm(100),y=rnorm(100), lab=sample(c("A", "B", "C", "D"), 100, replace=TRUE))
#' ggplot(dframe.num.sym.lab, aes(x=x,y=y, color=lab)) + geom_point() + ggAutoColorScale(dframe.num.sym.lab[,"lab"])
ggAutoColorScale <- function(x, method="color", symmetric=TRUE){
  if (!is.element(method, c("color", "colour", "fill"))) error("invalid scale method")
  res <- NULL
  params <- list()
  params[["na.value"]] <- "#C0C0C0"
  if (is.numeric(x)){
    
    lims <- c(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
    params[["colors"]] <- colpal.cont(n=9, name="viridis")
    if (symmetric){
      ctrVal <- 0L
      if(any(x < ctrVal, na.rm=TRUE) && any(x > ctrVal, na.rm=TRUE)){
        xDiff <- x - ctrVal
        lims <- max(abs(xDiff), na.rm=TRUE)
        lims <- c(ctrVal-lims, ctrVal+lims)
        params[["colors"]] <- rev(colpal.cont(n=9, name="cb.RdYlBu"))
      } 
    }
    params[["limits"]] <- lims

    mname <- paste0("scale_", method, "_gradientn")
  } else {
    res[["scale_method"]] <- "manual"
    if (!is.factor(x)){
      x <- factor(x)
    }
    lvls <- levels(x)
    cs <- rep(colpal.mu.cat, length.out=length(lvls))
    params[["limits"]] <- lvls
    params[["values"]] <- cs
    mname <- paste0("scale_", method, "_manual")
  }
  res <- do.call(mname, params)
  return(res)
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
#' @param dimPreset presets for figure dimensions. Possible values are \code{NULL} (don't use a preset; default),
#'                \code{"slide_nuc_wide_full"} (full slide using Fabian's wide nucleosome template),
#'                \code{"slide_nuc_wide_half"} (half a slide using Fabian's wide nucleosome template)
#'                Overwrites \code{width}, \code{height} and \code{units}.
#' @param useDingbats see \code{?ggsave}
#' @param ...     see \code{?ggsave}
#' @return result of \code{ggsave} command
#' @export
#' @aliases ggsave4doc
ggsave4doc <- function(fn, plot=last_plot(), width=192, height=192, units="mm", family="Helvetica", dimPreset=NULL, useDingbats=FALSE,...){
  if (!is.null(dimPreset)){
    if (dimPreset == "slide_nuc_wide_full"){
      width <- 338.6
      height <- 158
      units <- "mm"
    } else if (dimPreset == "slide_nuc_wide_half"){
      width <- 338.6/2
      height <- 158
      units <- "mm"
    } else {
      warning(paste("Undefined preset:", dimPreset))
    }
  }
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
