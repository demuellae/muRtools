################################################################################
# dimension reduction
################################################################################
#' getDimRedCoords.pca
#' 
#' Get dimension reduction coordinates (PCA)
#' @param X   A matrix on which the dimension reduction is to be performed
#' @param components   principal component to be returned
#' @return a matrix containing two columns for the reduced dimensions and the same number of
#'         rows as \code{X}
#' @author Fabian Mueller
#' @export 
getDimRedCoords.pca <- function(X, components=c(1,2)){
	numNA <- colSums(is.na(X) | is.infinite(X))
	has.noNA <- numNA==0
	if (any(numNA>0)){
		logger.info(c("retained",sum(has.noNA),"of",ncol(X),"observations because the remaining ones contained NAs/Inf"))
		X <- X[,has.noNA]
	}
	pca <- prcomp(X, center = TRUE, scale. = FALSE)
	coords <- pca$x[,components]
	rownames(coords) <- rownames(X)
	percVar <- 100 *(pca$sdev)^2 / sum(pca$sdev^2)
	names(percVar) <- colnames(pca$x)
	attr(coords, "percVar") <- percVar[components]
	return(coords)
}
#' getDimRedCoords.mds
#' 
#' Get dimension reduction coordinates (Multidimensional Scaling)
#' @param X   A matrix on which the dimension reduction is to be performed
#' @param distMethod   distance metric to be employed
#' @return a matrix containing two columns for the reduced dimensions and the same number of
#'         rows as \code{X}
#' @author Fabian Mueller
#' @export 
getDimRedCoords.mds <- function(X, distMethod="euclidean"){
	dist.matrix <- dist(X, method=distMethod)
	if (any(dist.matrix == 0)) {
		## Get rid of zeros in the distance matrix by adding an epsilon
		dist.matrix <- dist.matrix + min(dist.matrix[dist.matrix != 0]) / 1000
	}
	if (any(is.na(dist.matrix))){
		warning("Distance matrix contained NAs --> overwrite with maximum")
		dist.matrix[is.na(dist.matrix)] <- max(dist.matrix, na.rm=TRUE)
	}
	coords <- isoMDS(dist.matrix, k = 2, maxit = 100, trace = FALSE)$points
	colnames(coords) <- paste0("MDS", 1:ncol(coords))
	rownames(coords) <- rownames(X)
	return(coords)
}
#' getDimRedCoords.tsne
#' 
#' Get dimension reduction coordinates (t-SNE)
#' @param X   A matrix on which the dimension reduction is to be performed
#' @param distMethod   distance metric to be employed
#' @param dims dimensions to return from the reduction
#' @return a matrix containing two columns for the reduced dimensions and the same number of
#'         rows as \code{X}
#' @author Fabian Mueller
#' @export 
getDimRedCoords.tsne <- function(X, distMethod="euclidean", dims=c(1,2)){
	require(tsne)
	dist.matrix <- dist(X, method=distMethod)
	if (any(is.na(dist.matrix))){
		warning("Distance matrix contained NAs --> overwrite with maximum")
		dist.matrix[is.na(dist.matrix)] <- max(dist.matrix, na.rm=TRUE)
	}
	k <- max(dims[1:2])
	coords <- tsne(dist.matrix, k=k)
	colnames(coords) <- paste0("tSNE", 1:ncol(coords))
	rownames(coords) <- rownames(X)
	return(coords)
}
#' getDimRedPlot
#' 
#' Generate a plot from dimension reduction coordinates
#' @param coords   dimension reduction coordinates
#' @param annot    annotation matrix with the same number of rows as \code{coord}
#' @param colorCol name or index in the annotation matrix (\code{annot}) that should be used for coloring the points
#'                 if \code{colorCol} not supplied but \code{annot} is supplied, it defaults to the first annotation column
#' @param shapeCol name or index in the annotation matrix (\code{annot}) that should be used for point shapes
#'                 if \code{shapeCol} not supplied but \code{annot} is supplied and has more than one column, 
#'                 it defaults to the second annotation column
#' @param colScheme color sheme to be used in coloring the points
#' @param addLabels should observation labels be added to each point
#' @param addDensity should Gaussian Kernel density estimation be performed and the contour lines plotted for each color group
#' @param annot.text optional text to be added in the lower right corner of the plot
#' @return a \code{ggplot2} object containing the dimension reduction plot
#' @author Fabian Mueller
#' @export 
getDimRedPlot <- function(coords, annot=NULL, colorCol=NULL, shapeCol=NULL, colScheme=NULL, addLabels=FALSE, addDensity=FALSE, annot.text=NULL){
	if (!is.null(annot)){
		if (nrow(annot)!=nrow(coords)){
			stop("Non-matching number of rows for dimension reduction coordinates and annotation")
		}
		if (ncol(annot) > 0){
			if (is.null(colorCol)){
				colorCol <- 1L
			}
		} else {
			stop("annot must have at least one column")
		}
		if (ncol(annot) > 1){
			if (is.null(shapeCol)){
				shapeCol <- 2L
			}
		}
	}
	xLab <- colnames(coords)[1]
	yLab <- colnames(coords)[2]

	df2p <- data.frame(coords)
	if (!is.null(colorCol)) {
		if (length(colorCol)!=1){
			stop("colorCol must be of length 1")
		}
		if (is.numeric(colorCol)){
			if (!is.null(colnames(annot))){
				colorCol <- colnames(annot)[colorCol]
			}
		}
		if (is.character(colorCol)){
			df2p[,colorCol] <- annot[,colorCol]
		} else if (is.numeric(colorCol)){
			df2p[, "color"] <- annot[,colorCol]
			colorCol <- "color"
		} else {
			stop("invalid value for colorCol")
		}
	}
	if (!is.null(shapeCol)) {
		if (length(shapeCol)!=1){
			stop("shapeCol must be of length 1")
		}
		if (is.numeric(shapeCol)){
			if (!is.null(colnames(annot))){
				shapeCol <- colnames(annot)[shapeCol]
			}
		}
		if (is.character(shapeCol)){
			df2p[,shapeCol] <- annot[,shapeCol]
		} else if (is.numeric(shapeCol)){
			df2p[, "shape"] <- annot[,shapeCol]
			shapeCol <- "shape"
		} else {
			stop("invalid value for shapeCol")
		}
	}

	if (!is.null(rownames(coords))) df2p$observation <- rownames(coords)
	pp <- ggplot(df2p, aes_string(x=xLab, y=yLab, color=colorCol))
	if (!is.null(colScheme)){
		pp <- pp + scale_color_manual(na.value = "#C0C0C0", values=colScheme)
	}
	if (addDensity){
		#remove 1-sample groups
		grpCounts <- table(df2p[,colorCol])
		grps.1sample <- names(grpCounts)[grpCounts<2]
		df2p2 <- df2p[!(df2p[,colorCol] %in% grps.1sample),]
		pp <- pp + stat_density2d(data=df2p2, alpha=0.3)
	}
	if (!is.null(shapeCol)){
		pp <- pp + geom_point(aes_string(shape=shapeCol), size=3)
	} else {
		pp <- pp + geom_point(size=3)
	}
	if (addLabels && is.element("observation", colnames(df2p))){
		pp <- pp + geom_text(aes(label=observation), size=2)
	}
	# pp <- pp + coord_fixed()
	if (!is.null(annot.text)){
		pp <- pp + annotate("text", x=max(df2p[,xLab],na.rm=TRUE),y=min(df2p[,yLab],na.rm=TRUE),label=annot.text,parse=FALSE,hjust=1,vjust=1,size=4)
	}
	#add percent variance explained to axis labels
	percVar <- attr(coords,"percVar")
	if (!is.null(percVar)){
		pp <- pp + xlab(paste0(pp$labels$x, " (", round(percVar[pp$labels$x], 2), "%)"))
		pp <- pp + ylab(paste0(pp$labels$y, " (", round(percVar[pp$labels$y], 2), "%)"))
	}
	return(pp)
}
#' getDimRedPlot
#' 
#' Generate a plot from a feature matrix
#' @param X         feature matrix containing one row for each observation and one column for each feature
#' @param dimRedFun function to do dimension reduction. E.g. \code{getDimRedCoords.pca}, \code{getDimRedCoords.mds}, \code{getDimRedCoords.tsne},
#' @param annot     annotation matrix with the same number of rows as \code{X}
#' @param colorCol name or index in the annotation matrix (\code{annot}) that should be used for coloring the points
#'                 if \code{colorCol} not supplied but \code{annot} is supplied, it defaults to the first annotation column
#' @param shapeCol name or index in the annotation matrix (\code{annot}) that should be used for point shapes
#'                 if \code{shapeCol} not supplied but \code{annot} is supplied and has more than one column, 
#'                 it defaults to the second annotation column
#' @param colScheme color sheme to be used in coloring the points
#' @param addLabels should observation labels be added to each point
#' @param addDensity should Gaussian Kernel density estimation be performed and the contour lines plotted for each color group
#' @param annot.text optional text to be added in the lower right corner of the plot
#' @param ...       arguments to be passed on to \code{dimRedFun}
#' @return a \code{ggplot2} object containing the dimension reduction plot
#' @author Fabian Mueller
#' @export 
plotDimRed <- function(X, dimRedFun=getDimRedCoords.pca,
		annot=NULL, colorCol=NULL, shapeCol=NULL, colScheme=NULL, addLabels=FALSE, addDensity=FALSE, annot.text=NULL,
		...){
	coords <- dimRedFun(X, ...)
	pp <- getDimRedPlot(coords, annot=annot, colorCol=colorCol, shapeCol=shapeCol, colScheme=colScheme, addLabels=addLabels, addDensity=addDensity, annot.text=annot.text)
	return(pp)
}
#' plotAllDimRed
#' 
#' Generate a plots with multiple methods and parameter settings from a feature matrix
#' @param X           feature matrix containing one row for each observation and one column for each feature
#' @param fn.prefix   file prefix to be used for the resulting plots. 
#'                    If \code{NULL}, no plot is actually created, but a list of resulting plot objects is returned
#' @param fn.suffix   file suffix to be used for the resulting plots
#' @param annot       annotation matrix with the same number of rows as \code{X}
#' @param distMethods distance methods for MDS and t-SNE
#' @param width       width of the resulting plot
#' @param height      height of the resulting plot
#' @param ...       arguments to be passed on to \code{getDimRedPlot}
#' @return (invisibly) a list of lists containing the created plots as ggplot objects and additional info for each plot
#' 
#' @details
#' Currently, PCA, MDS and t-SNE are employed by default with euclidean and manhattan distance metrics where applicable
#' @author Fabian Mueller
#' @export 
plotAllDimRed <- function(X, fn.prefix=NULL, fn.suffix="", annot=NULL, distMethods=c(euc="euclidean",man="manhattan"), width=10, height=10,
		 ...){
	
	res <- list()
	pp <- plotDimRed(X, dimRedFun=getDimRedCoords.pca, annot=annot, ...)
	res <- c(res, list(list(plot=pp, method="pca", dist=NA)))
	for (i in 1:length(distMethods)){
		distMeth <- distMethods[i]
		distMethName <- names(distMethods)[i]
		pp <- plotDimRed(X, dimRedFun=getDimRedCoords.mds, annot=annot, distMethod=distMethods[i], ...)
		res <- c(res, list(list(plot=pp, method="mds", dist=distMethName)))
		pp <- plotDimRed(X, dimRedFun=getDimRedCoords.tsne, annot=annot, distMethod=distMethods[i], ...)
		res <- c(res, list(list(plot=pp, method="tsne", dist=distMethName)))
	}
	if (!is.null(fn.prefix)) {
		suff <- fn.suffix
		if (nchar(fn.suffix)>0) suff <- paste0("_",fn.suffix)
		for (ple in res){
			fName <- paste0(fn.prefix,"_", ple$method)
			if (!is.na(ple$dist)) fName <- paste0(fName, "_", ple$dist)
			fName <- paste0(fName, suff, ".pdf")
			ggsave(fName, ple$plot, width=width, height=height)
		}
	}
	invisible(res)
}
