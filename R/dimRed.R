################################################################################
# dimension reduction
################################################################################
#' getDimRedCoords.pca
#' 
#' Get dimension reduction coordinates (PCA)
#' @param X   A matrix on which the dimension reduction is to be performed
#' @param components   principal component to be returned
#' @param method Method/package to be used for computing principal components. Currently \code{prcomp} and \code{irlba} are supported.
#' @param ... other arguments passed on to the PCA method
#' @return a matrix containing two columns for the reduced dimensions and the same number of
#'         rows as \code{X}
#' @author Fabian Mueller
#' @export 
getDimRedCoords.pca <- function(X, components=c(1,2), method="prcomp", ...){
	csFun <- colSums
	tFun <- t
	if (is.character(attr(class(X), "package")) && attr(class(X), "package")=="Matrix"){
		csFun <- Matrix::colSums
		tFun <- Matrix::t
	}
	numNA <- csFun(is.na(X) | is.infinite(X))
	has.noNA <- numNA==0
	if (any(numNA>0)){
		logger.info(c("retained",sum(has.noNA),"of",ncol(X),"features because the remaining ones contained NAs/Inf"))
		X <- X[,has.noNA]
	}

	if (method=="prcomp"){
		pca <- prcomp(X, center = TRUE, scale. = FALSE, ...)
		coords <- pca$x[,components, drop=FALSE]
		rownames(coords) <- rownames(X)
		percVar <- 100 *(pca$sdev)^2 / sum(pca$sdev^2)
		names(percVar) <- colnames(pca$x)
		attr(coords, "percVar") <- percVar[components]
		attr(coords,"PCAclass") <-"PCcoord"
		attr(coords,"PCAmethod") <- "prcomp"
	} else if (method=="irlba"){
		require(irlba)
		nComps <- max(components)
		pca <- prcomp_irlba(X, n=nComps, center = TRUE, scale. = FALSE, ...)
		coords <- pca$x[,components, drop=FALSE]
		rownames(coords) <- rownames(X)
		percVar <- 100 *(pca$sdev)^2 / sum(pca$sdev^2)
		names(percVar) <- colnames(pca$x)
		attr(coords, "percVar") <- percVar[components]
		attr(coords,"PCAclass") <-"PCcoord"
		attr(coords,"PCAmethod") <- "irldba"
	} else if (method=="irlba_svd") {
		require(irlba)
		nComps <- max(components)
		svdRes <- irlba::irlba(tFun(X), nComps, nComps, maxit=1000, ...)
		svdDiag <- matrix(0, nrow=length(svdRes$d), ncol=length(svdRes$d))
		diag(svdDiag) <- svdRes$d
		coords <- t(svdDiag %*% t(svdRes$v))
		rownames(coords) <- rownames(X)
		colnames(coords) <- paste0('PC', 1:ncol(coords))

		vv <- svdRes$d^2/nrow(X)
		attr(coords, "percVar") <- 100 * vv/sum(vv)
		attr(coords, "PCAclass") <-"PCcoord"
		attr(coords, "PCAmethod") <- "irldba_svd"
		attr(coords, "SVD_D") <- svdDiag
		attr(coords, "SVD_U") <- svdRes$u
		attr(coords, "SVD_V") <- svdRes$v
		# projection of a new matrix (M) into PC space: coord_M <- M %*% svdRes$u
	} else {
		stop(paste("Unknown PCA method:", method))
	}
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
	require(MASS)
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
#' getDimRedCoords.umap
#' 
#' Get dimension reduction coordinates using the UMAP method
#' @param X   A matrix on which the dimension reduction is to be performed
#' @param distMethod   distance metric to be employed
#' @param dims dimensions to return from the reduction
#' @param ... parameters passed on to \code{uwot::umap()}
#' @return a matrix containing two columns for the reduced dimensions and the same number of
#'         rows as \code{X}
#' @author Fabian Mueller
#' @export 
getDimRedCoords.umap <- function(X, distMethod="euclidean", dims=c(1,2), ...){
	require(uwot)
	dotArgs <- list(...)
	numNA <- colSums(is.na(X) | is.infinite(X)) # colSums(is.na(X) | is.infinite(X))
	has.noNA <- numNA==0
	if (any(numNA>0)){
		logger.info(c("retained",sum(has.noNA),"of",ncol(X),"features because the remaining ones contained NAs/Inf"))
		X <- X[,has.noNA]
	}
	k <- max(dims)
	callArgL <- list(
		X,
		n_components=k,
		metric=distMethod,
		ret_model=TRUE
	)
	callArgL <- c(callArgL, list(n_neighbors=15))
	if (is.element("n_neighbors", names(dotArgs))){
		callArgL[["n_neighbors"]] <- dotArgs[["n_neighbors"]]
		dotArgs[["n_neighbors"]] <- NULL
	}
	if (callArgL[["n_neighbors"]] > nrow(X)){
		nn <- as.integer(nrow(X)/2)
		logger.warning(c("UMAP: number of neighbors can't be > N_samples. --> reducing to", nn, "neighbors"))
		callArgL[["n_neighbors"]] <- nn
	}
	callArgL <- c(callArgL, dotArgs)

	uRes <- do.call("umap", callArgL)
	coords <- uRes$embedding[,dims]
	colnames(coords) <- paste0("UMAP", 1:ncol(coords))
	rownames(coords) <- rownames(X)
	attr(coords,"umapRes") <- uRes
	return(coords)
}
# getDimRedCoords.umap <- function(X, distMethod="euclidean", dims=c(1,2)){
# 	require(reticulate)
# 	require(umap)
# 	if (!py_module_available("umap")){
# 		stop("could not load Python module 'umap-learn'")
# 		# if the module cannot be loaded this can be because not the correct python environment is loaded (check py_config())
# 		# also: umap requires the correct version of the C++ compiler. So, be sure to load GCC beforehand
# 	}
# 	numNA <- colSums(is.na(X) | is.infinite(X)) # colSums(is.na(X) | is.infinite(X))
# 	has.noNA <- numNA==0
# 	if (any(numNA>0)){
# 		logger.info(c("retained",sum(has.noNA),"of",ncol(X),"features because the remaining ones contained NAs/Inf"))
# 		X <- X[,has.noNA]
# 	}
# 	k <- max(dims[1:2])
# 	umapRes <- umap(X, method="umap-learn", n_components=k, metric=distMethod)
# 	coords <- umapRes$layout
# 	colnames(coords) <- paste0("UMAP", 1:ncol(coords))
# 	rownames(coords) <- rownames(X)
# 	return(coords)
# }

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
#' @param colScheme color sheme to be used in coloring the points. can be a character vector with the supplied colors. Alternatively, if it is a one-element character vector \code{"[auto]"} the color scheme will be selected automatically using \code{muRtools::ggAutoColorScale}. If \code{NULL}, ggplots default color scheme will be used.
#' @param ptSize   size of the points in the scatterplot
#' @param addLabels should observation labels be added to each point
#' @param addDensity should Gaussian Kernel density estimation be performed and the contour lines plotted for each color group
#' @param addVoronoi should a Voronoi tessalation grid (based on \code{colorCol}) be added to the plot
#' @param annot.text optional text to be added in the lower right corner of the plot
#' @param orderCol name or index in the annotation matrix (\code{annot}) that should be used for ordering the points. If not \code{NULL} Points will be ordered increasingly by their value, i.e. higher-valued points are plottet over lower-valued points
#' @param facetCols name (string) of columns to be used for faceting the resulting plot. Each facet will contain all the points not in the facet as grey points.
#' @return a \code{ggplot2} object containing the dimension reduction plot
#' @author Fabian Mueller
#' @export 
#' 
#' @examples
#' \donttest{
#' df <- data.frame(
#' 	x = c(rnorm(20, mean=0, sd=0.2), rnorm(10, mean=1, sd=0.4), rnorm(15, mean=1, sd=0.2)),
#' 	y = c(rnorm(20, mean=0, sd=0.2), rnorm(10, mean=1, sd=0.4), rnorm(15, mean=0.5, sd=0.3)),
#' 	group = rep(c("group1", "group2", "group3"), times=c(20,10,15)),
#' 	stringsAsFactors=FALSE
#' )
#' getDimRedPlot(df[,c("x", "y")], annot=df[,c("group"), drop=FALSE], colorCol="group")
#' getDimRedPlot(df[,c("x", "y")], annot=df[,c("group"), drop=FALSE], colorCol="group", addDensity=TRUE)
#' }
getDimRedPlot <- function(coords, annot=NULL, colorCol=NULL, shapeCol=NULL, colScheme="[auto]", ptSize=3, addLabels=FALSE, addDensity=FALSE, addVoronoi=FALSE, annot.text=NULL, orderCol=NULL, facetCols=NULL){
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
	colorNumeric <- FALSE
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
		if (is.numeric(df2p[,colorCol])){
			colorNumeric <- TRUE
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
		} else if (is.logical(shapeCol) && !shapeCol) {
			shapeCol <- NULL
		} else {
			stop("invalid value for shapeCol")
		}
		if (is.numeric(df2p[,shapeCol])){
			warning("Currently only non-numeric columns are supported for dimRed shapes. --> converting to factor")
			df2p[,shapeCol] <- factor(df2p[,shapeCol])
		}
	}

	if (!is.null(orderCol)) {
		if (length(orderCol)!=1){
			stop("colorCol must be of length 1")
		}
		if (is.numeric(orderCol)){
			if (!is.null(colnames(annot))){
				orderCol <- colnames(annot)[orderCol]
			}
		}
		if (is.character(orderCol) || is.numeric(orderCol)){
			oo <- order(df2p[,orderCol], decreasing=FALSE, na.last=FALSE)
		} else {
			stop("invalid value for orderCol")
		}
		df2p <- df2p[oo,]
	}

	if (!is.null(facetCols)){
		if (length(facetCols) > 2){
			stop("Too many facet columns")
		}
		for (cn in facetCols){
			df2p[,cn] <- annot[,cn]
		}
	}

	if (!is.null(rownames(coords))) df2p$observation <- rownames(coords)
	pp <- ggplot(df2p, aes_string(x=xLab, y=yLab, color=colorCol))
	if (!is.null(colScheme)){
		if (is.character(colScheme) && length(colScheme)==1 && colScheme=="[auto]"){
			pp <- pp + ggAutoColorScale(df2p[,colorCol])
		} else {
			if (colorNumeric){
				pp <- pp + scale_color_gradientn(colours=colScheme, na.value = "#C0C0C0")
			} else {
				pp <- pp + scale_color_manual(na.value = "#C0C0C0", values=colScheme)
			}
		}
	}
	if (addDensity){
		if (colorNumeric){
			warning("Currently only non-numeric columns are supported for dimRed coloring in combination with density contours. --> skipping density contours")
		} else {
			#remove 1-sample groups
			grpCounts <- table(df2p[,colorCol])
			grps.1sample <- names(grpCounts)[grpCounts<2]
			df2p2 <- df2p[!(df2p[,colorCol] %in% grps.1sample),]
			pp <- pp + stat_density2d(data=df2p2, alpha=0.3)
		}
	}
	if (addVoronoi){
		if (colorNumeric){
			warning("Currently only non-numeric columns are supported for dimRed coloring in combination with density contours. --> skipping density contours")
		} else {
			require(deldir)
			center.x <- tapply(df2p[,xLab], df2p[,colorCol], FUN=function(x){mean(x, na.rm=TRUE)})
			center.y <- tapply(df2p[,yLab], df2p[,colorCol], FUN=function(x){mean(x, na.rm=TRUE)})

			voronoi <- deldir(center.x, center.y,
				rw=c(min(df2p[,xLab], na.rm=TRUE), max(df2p[,xLab], na.rm=TRUE), min(df2p[,yLab], na.rm=TRUE), max(df2p[,yLab], na.rm=TRUE))
			)
			pp <- pp + geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), size=1, data=voronoi$dirsgs, linetype=1, color="#C0C0C0")
		}
	}
	if (!is.null(facetCols)){
		pp <- pp + geom_point(data=df2p[,!(colnames(df2p) %in% facetCols)], color="#C0C0C0", size=ptSize, shape=16)
	}
	if (!is.null(shapeCol)){
		pp <- pp + geom_point(aes_string(shape=shapeCol), size=ptSize) + scale_shape_manual(values=c(19,15,17,4,3,18,8,1,0,2,6))
	} else {
		pp <- pp + geom_point(size=ptSize, shape=16)
	}
	if (addLabels && is.element("observation", colnames(df2p))){
		pp <- pp + geom_text(aes(label=observation), size=2)
	}
	# pp <- pp + coord_fixed()
	if (!is.null(annot.text)){
		pp <- pp + annotate("text", x=max(df2p[,xLab],na.rm=TRUE),y=min(df2p[,yLab],na.rm=TRUE),label=annot.text,parse=FALSE,hjust=1,vjust=1,size=4)
	}
	if (!is.null(facetCols)){
		if (length(facetCols)==1){
			nr <- ceiling(sqrt(length(unique(df2p[,facetCols]))))
			pp <- pp + facet_wrap(as.formula(paste0("~", facetCols)), nrow=nr)
		} else {
			pp <- pp + facet_grid(as.formula(paste(facetCols, collapse="~")))
		}	
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
#' @param ptSize   size of the points in the scatterplot
#' @param addLabels should observation labels be added to each point
#' @param addDensity should Gaussian Kernel density estimation be performed and the contour lines plotted for each color group
#' @param annot.text optional text to be added in the lower right corner of the plot
#' @param ...       arguments to be passed on to \code{dimRedFun}
#' @return a \code{ggplot2} object containing the dimension reduction plot
#' @author Fabian Mueller
#' @export 
plotDimRed <- function(X, dimRedFun=getDimRedCoords.pca,
		annot=NULL, colorCol=NULL, shapeCol=NULL, colScheme=NULL, ptSize=3, addLabels=FALSE, addDensity=FALSE, annot.text=NULL,
		...){
	coords <- dimRedFun(X, ...)
	pp <- getDimRedPlot(coords, annot=annot, colorCol=colorCol, shapeCol=shapeCol, colScheme=colScheme, ptSize=ptSize, addLabels=addLabels, addDensity=addDensity, annot.text=annot.text)
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

################################################################################
# Associations
################################################################################
#' getAssocTestRes.pca
#' 
#' Test associations of annotations with principal components (PCA)
#' @param X       A matrix on which the dimension reduction is to be performed. Alternatively, it can be a matrix of PC coordinates computed by \code{getDimRedCoords.pca}.
#' @param ph      annotation table for the datapoints. The columns of this table will be used to test associations with the PCs.
#'                Should be a \code{matrix} or \code{data.frame}.
#' @param nComp   number of PCs to be considered
#' @param nPerm   number of permutation tests to be conducted if an annotation in \code{ph} is numeric (i.e. a correlation permutation test is performed)
#' @return A nested list of tested associations one element for each column in \code{ph} (1st level), each PC (2nd level).
#'         Each element is again a list with the name of the test being used (\code{test}), the test statistic (\code{statistic}) and p-value (\code{pvalue})
#' @author Fabian Mueller
#' @export 
getAssocTestRes.pca <- function(X, ph, nComp=10, nPerm=1000){
	if (!is.matrix(X)) logger.error("Invalid X object. Expected matrix")
	if (!is.matrix(ph) && !is.data.frame(ph)) logger.error("Invalid ph object. Expected matrix or data frame")
	if (nrow(ph)!=nrow(X)) logger.error("Invalid ph object. Expected matrix or data frame WITH MATCHING DIMENSIONS")
	nPoints <- nrow(X)

	nComp <- min(nComp, nPoints)
	# check if X is already a PCcoord object
	if (is.null(attr(X,"PCAclass")) || attr(X,"PCAclass")!="PCcoord"){
		coords <- getDimRedCoords.pca(X, components=1:nComp)
	} else {
		coords <- X
		nComp <- ncol(coords)
	}

	#factorize matrix
	for (j in which(sapply(ph, FUN=function(x){is.character(x) || is.logical(x)}))){
		ph[,j] <- factor(ph[,j])
	}
	permMat <- NULL
	if (nPerm > 0 && any(!sapply(ph, is.factor))) {
		permMat <- mapply(sample, rep(nPoints, times=nPerm))
		# logger.info(c("Created", ncol(permMat), "sample permutations for correlation permutation testing"))
	}
	assocL <- lapply(1:ncol(ph), FUN=function(i){
		rr <- lapply(1:nComp, FUN=function(j){
			testAssoc(ph[,i], coords[,j], permMat)
		})
		names(rr) <- colnames(coords)
		return(rr)
	})
	names(assocL) <- colnames(ph)
	# pVals <- sapply(assocL, FUN=function(x){sapply(x, FUN=function(y){y[["pvalue"]]})})
	# testNames <- sapply(assocL, FUN=function(x){sapply(x, FUN=function(y){y[["test"]]})})
	return(assocL)
}

################################################################################
# Grid search UMAP parameters
################################################################################
#' umapParamGridReport
#' 
#' Generate a report with plots of UMAP dimension reduction plots for parameter combinations
#' @param X         feature matrix containing one row for each observation and one column for each feature
#' @param outDir    output directory
#' @param metric    parameters passed on to \code{uwot::umap()}
#' @param min_dist  parameters passed on to \code{uwot::umap()}
#' @param n_neighbors parameters passed on to \code{uwot::umap()}
#' @param ...       parameters passed on to \code{getDimRedPlot}
#' @return a \code{muReportR} report (HTML) showing dimension reduction plots for the grid search
#' 
#' @author Fabian Mueller
#' @export
umapParamGridReport <- function(X, outDir, metric=c("euclidean"), min_dist=c(0.01, 0.05, seq(0.1, 0.9, by=0.1)), n_neighbors=c(5,15,25,50), ...){
	require(muReportR)
	require(ggplot2)
	require(muLogR)
	umap_params <- expand.grid(metric, min_dist, n_neighbors, stringsAsFactors=FALSE)
	colnames(umap_params) <- c("metric", "min_dist", "n_neighbors")
	coordL <- lapply(1:nrow(umap_params), FUN=function(i){
		logger.status(c("Computing UMAP dimensions (", i, "of", nrow(umap_params), ")"))
		getDimRedCoords.umap(X, distMethod=umap_params[i,"metric"], dims=c(1,2), min_dist=umap_params[i,"min_dist"], n_neighbors=umap_params[i,"n_neighbors"])
	})
	rr <- muReportR::createReport(file.path(outDir, "umapParams.html"), "UMAP parameter grid search", page.title = "UMAP parameters", init.configuration=TRUE, theme="stanford")
	txt <- "Parameter grid search for UMAP parameters"
	rr <- muReportR::addReportSection(rr, "Dimension reduction plots", txt, level=1L, collapsed=FALSE)
	
	umap_params_strings <- umap_params
	for (cn in colnames(umap_params)){
		umap_params_strings[[cn]] <- normalize.str(as.character(umap_params_strings[[cn]]), return.camel=TRUE)
	}

	plotL <- lapply(seq_along(coordL), FUN=function(i){
		logger.status(c("Plotting reduced dimension scatterplot (", i, "of", nrow(umap_params), ")"))
		pp <-  getDimRedPlot(coordL[[i]], ...) + coord_fixed()
		figFn <- paste("umap", umap_params_strings[i,"metric"], umap_params_strings[i,"min_dist"], umap_params_strings[i,"n_neighbors"], sep="_")
		repPlot <- muReportR::createReportGgPlot(pp, figFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
		repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
		return(repPlot)
	})
	figSettings.metric <- unique(umap_params[,"metric"])
	names(figSettings.metric) <- unique(umap_params_strings[,"metric"])
	figSettings.min_dist <- unique(umap_params[,"min_dist"])
	names(figSettings.min_dist) <- unique(umap_params_strings[,"min_dist"])
	figSettings.n_neighbors <- unique(umap_params[,"n_neighbors"])
	names(figSettings.n_neighbors) <- unique(umap_params_strings[,"n_neighbors"])
	figSettings <- list(
		"metric" = figSettings.metric,
		"min_dist" = figSettings.min_dist,
		"n_neighbors" = figSettings.n_neighbors
	)
	rr <- muReportR::addReportFigure(rr, "UMAP dimension reduction", plotL, figSettings)

	muReportR::off(rr)
	invisible(rr)
}
# umapParamGridReport(dre$pcaCoord, "reports", metric=c("euclidean"), min_dist=c(0.01, 0.51), n_neighbors=c(5,15,30), annot=data.frame(cluster=dre$clustAss), ptSize=0.25)
# umapParamGridReport(dre$pcaCoord, "reports", annot=data.frame(cluster=dre$clustAss), ptSize=0.25)
