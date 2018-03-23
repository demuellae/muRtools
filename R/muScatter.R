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

########################################################################################################################

#' densRanks
#'
#' Rank the points accordind to density of the region they fall in. Densities are computed
#' as Kernel Density estimates. The method and parameters are implemented in analogy to
#' \code{grDevices::densCols}
#' @param x x-coordinate
#' @param y y-coordinate
#' @param nbin number of bins
#' @param bandwidth bandwidth
#' @author Fabian Mueller (RnBeads)
densRanks <- function (x, y = NULL, nbin = 128, bandwidth) 
{
    xy <- xy.coords(x, y)
    select <- is.finite(xy$x) & is.finite(xy$y)
    x <- cbind(xy$x, xy$y)[select, ]
    map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
    mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
    xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
    ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
    dens <- map$fhat[cbind(xbin, ybin)]
    dens[is.na(dens)] <- 0

    res <- rep(NA_integer_, length(select))
    rrs <- rank(dens,ties.method="max")
    res[select] <- rrs
    res
}

########################################################################################################################

#' create.densityScatter
#' 
#' Creates a density scatterplot highlighting points in sparsely populated plot regions
#' as well as points marked as special in a seperate color
#' @param df2p \code{data.frame} to be plotted. Only the fist two columns are taken into account as
#' 			x and y coordinates respectively
#' @param is.special boolean vector of length equal to the number of rows in \code{df2p}. Specifies
#' 			which points should be highlighed seperately in a different color
#' @param dens.subsample if the number of points exceeds this number, subsample the number of points for the
#'			density estimation to that number. Any non-numeric value disables subsampling.
#' @param dens.special Flag indicating whether the points of the special population should be colored
#' 			according to their density
#' @param sparse.points Either percentage (\code{<=1,>=0}) or the absolute number 
#'          of points in the sparsely populated area that should be drawn seperately. A value of 0 means that these points
#'			will not be drawn.
#' @param dens.n passed on to \code{ggplot2::stat_density2d}: argument: \code{n}
#' @param add.text.cor flag indicating whether a text token with the correlation coefficient should be included in the lower
#'          right corner of the plot
#' @return \code{ggplot} object
#' @author Fabian Mueller (RnBeads)
#' @export
#' @examples
#' \donttest{
#' d <- data.frame(x=rnorm(1000),y=rnorm(1000))
#' s <- rep(FALSE,1000)
#' s[sample(1:length(s),100)] <- TRUE
#' create.densityScatter(d,s)
#' }
create.densityScatter <- function(df2p,is.special=NULL,dens.subsample=FALSE,dens.special=TRUE,
		sparse.points=0.01,dens.n=100,add.text.cor=FALSE){
	if (!(is.numeric(sparse.points) && sparse.points>=0)) {
		stop("Invalid parameter value: sparse.points")
	}
	if (!is.null(is.special)) is.special[is.na(is.special)] <- FALSE
	if (sum(is.special)<1){
		is.special <- NULL
	}
	if (!is.null(is.special)){
		df2p$is.special <- is.special
	}
	if (is.null(df2p) || nrow(df2p)<1){
		logger.warning(c("Could not create density scatterplot"))
		pp <- rnb.message.plot("Could not create plot")
		return(pp)
	}
	df2p <- na.omit(df2p)
	if (is.null(df2p) || nrow(df2p)<1){
		logger.warning(c("Could not create density scatterplot (NA omission removed all entries)"))
		pp <- rnb.message.plot("Could not create plot")
		return(pp)
	}
	df2p.sub <- df2p
	dens.ranks <- NULL
	tryCatch(
		dens.ranks <- densRanks(x=df2p[,1],y=df2p[,2]),
		error=function(ee){
			logger.warning(c("Could not assess density ranking:",ee$message))
		}
	)
	if (is.numeric(dens.subsample) && dens.subsample>0){
		ss <- as.integer(dens.subsample)
		if (nrow(df2p) > ss) {
			df2p.sub <- df2p[sample(nrow(df2p),ss),]
		}
	}

	#the standard bandwith function of MASS::kde2d is unstable when looking at
	#distributions with very low variance. Here's a more stable version
	stable.bandwidth.fun <- function(x,eps=1e-4){
	    r <- quantile(x, c(0.25, 0.75))
	    h <- (r[2] - r[1])/1.34
	    if (h==0) h <- eps
	    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
	}
	stable.h <- c(stable.bandwidth.fun(df2p.sub[,1]),stable.bandwidth.fun(df2p.sub[,2]))

	if (is.null(dens.ranks)){
  		pp <- rnb.message.plot("Could not assess density")
  	} else {		
  		pp <- ggplot(df2p.sub) + aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2]) + 
		  stat_density2d(geom="tile", fill=DENS.COLORS.LOW[1], aes(,alpha=..density..^0.25), contour=FALSE, n=dens.n, h=stable.h) +
		  scale_alpha(range = c(0.0, 1))
		if (sparse.points > 0){
			if (sparse.points <= 1){
				thres <- ceiling(nrow(df2p)*sparse.points)
			} else {
				thres <- sparse.points
			}
			df2p.loose <- df2p[dens.ranks<=thres,]#the sub data.frame in of the least dens points
			pp <- pp + geom_point(data=df2p.loose,aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2]),colour=DENS.COLORS.LOW[1],size=0.4)
		}
		if (!is.null(is.special)){
			df2p.special <- df2p[df2p$is.special,]
			colors.dmp <- DENS.COLORS.LOW[2]
			if (dens.special && nrow(df2p.special) > 1){
				tryCatch(
					colors.dmp   <- densCols(x=df2p.special[,1],y=df2p.special[,2],colramp = colorRampPalette(c(DENS.COLORS.LOW[2],DENS.COLORS.HIGH[2]))),
					error=function(ee){
						logger.warning(c("Could not assess density colors using densCols:",ee$message))
					}
				)
			}
			df2p.special$color <- colors.dmp

			pp <- pp + geom_point(data=df2p.special,aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2],colour="color"),size=1) + scale_color_identity()
		}
		if (add.text.cor) {
			cc <- cor(df2p[,1],df2p[,2],use="pairwise.complete.obs")
			txt.cor <- paste0('rho',paste0("==",round(cc,4)))
			pp <- pp + annotate("text", x=max(df2p[,1],na.rm=TRUE),y=min(df2p[,2],na.rm=TRUE),label=txt.cor,parse=TRUE,hjust=1,vjust=1,size=4)
		}
		pp <- pp + theme(legend.position="none")
  	}
	return(pp)
}
