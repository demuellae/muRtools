################################################################################
# Statistical testing
################################################################################

#' testAssoc
#' 
#' Tests for association between two vectors. Based on \code{RnBeads:::test.traits}
#'
#' @param x           Sample values for the first trait. This must be a vector of type \code{factor}, \code{integer} or
#'                    \code{numeric}.
#' @param y           Sample values for the second trait. This must be a vector of type \code{factor}, \code{integer} or
#'                    \code{numeric}.
#' @param permMat     Matrix of sample permutations (indices of x that will be used in permutation tests) in case none of the traits is a \code{factor}, and thus
#'                    permutation-based p-value from correlations is computed. If this parameter is \code{NULL} and
#'                    both \code{x} and \code{y} are sequences of numbers, no p-value is calculated.
#' @return            List of four elements:
#'                    \describe{
#'                      \item{error}{Error, if any, that prevented this function from computing a p-value for trait
#'                           association.}
#'                      \item{test}{Type of test performed. This is one of \code{"Fisher"}, \code{"Wilcoxon"},
#'                           \code{"Kruskal-Wallis"}, \code{"Correlation"} or \code{NA}. The last value indicates that
#'                           the traits cannot be tested for association.}
#'                      \item{correlation}{Value of the pearson correlation coefficient between \code{x} and \code{y},
#'                           or \code{NA} if any of them is \code{factor}.}
#'                      \item{pvalue}{Calculated p-value, or \code{NA} if the traits cannot be tested for association.}
#'                    }
#' @export
testAssoc <- function(x, y, permMat=NULL) {
	result <- list(
			"error" = as.character(NA),
			"test" = as.character(NA),
			"statistic" = as.double(NA),
			"pvalue" = as.double(NA))

	if (length(x) != length(y)) logger.error(c("x and y must have matching lengths"))
	if (class(x) == "Date") { x <- as.integer(x) }
	if (class(y) == "Date") { y <- as.integer(y) }
	if (is.character(x) || is.logical(x)) { x <- factor(x) }
	if (is.character(y) || is.logical(y)) { y <- factor(y) }
	disc <- is.na(x) | is.na(y)
	if (sum(disc) > 0){
		logger.warning(c("Discarding", sum(disc), "observations with NA (testAssoc)"))
		inds <- which(!disc)
		x <- x[inds]
		y <- y[inds]
	}
	
	if (length(x) < 2) {
		result[["error"]] <- "not enough shared values between x and y"
		return(result)
	}
	if (is.factor(x)) {
		x <- as.factor(as.character(x))
		if (nlevels(x) < 2) {
			## Not enough categories in y
			result[["error"]] <- "not enough categories in x"
			return(result)
		}
	}
	if (is.factor(y)) {
		y <- as.factor(as.character(y))
		if (nlevels(y) < 2) {
			## Not enough categories in y
			result[["error"]] <- "not enough categories in y"
			return(result)
		}
	}

	## Perform a test or compute correlation
	get.stat <- function(expr, statName) { tryCatch(suppressWarnings(expr[[statName]]), error = function(er) { as.double(NA) }) }
	if (is.factor(x)) {
		if (is.factor(y)) {
			simulate <- (nlevels(x) > 2 || nlevels(y) > 2)
			testRes <- fisher.test(x, y, conf.int = FALSE, simulate.p.value = simulate, B = 50000)
			result[["test"]] <- "Fisher"
			result[["pvalue"]] <- get.stat(testRes, "p.value")
			result[["statistic"]] <- get.stat(testRes, "estimate")
		} else if (nlevels(x) == 2) {
			result[["test"]] <- "Wilcoxon"
			values <- tapply(y, x, identity)
			testRes <- wilcox.test(values[[1]], values[[2]], alternative = "two.sided")
			result[["pvalue"]] <- get.stat(testRes, "p.value")
			result[["statistic"]] <- get.stat(testRes, "statistic")
		} else {
			testRes <- kruskal.test(y, x)
			result[["test"]] <- "Kruskal-Wallis"
			result[["pvalue"]] <- get.stat(testRes, "p.value")
			result[["statistic"]] <- get.stat(testRes, "statistic")
		}
	} else if (is.factor(y)) {
		if (nlevels(y) == 2) {
			result[["test"]] <- "Wilcoxon"
			values <- tapply(x, y, identity)
			testRes <- wilcox.test(values[[1]], values[[2]], alternative = "two.sided")
			result[["pvalue"]] <- get.stat(testRes, "p.value")
			result[["statistic"]] <- get.stat(testRes, "statistic")
		} else {
			testRes <- kruskal.test(x, y)
			result[["test"]] <- "Kruskal-Wallis"
			result[["pvalue"]] <- get.stat(testRes, "p.value")
			result[["statistic"]] <- get.stat(testRes, "statistic")
		}
	} else {
		result[["test"]] <- "Pearson correlation"
		if (is.null(permMat)) {
			result[["statistic"]] <- cor(x, y)
		} else {
			N <- length(x)
			values <- apply(permMat, 2, function(i) { cor(x[i[i <= N]], y) })
			result[["statistic"]] <- values[1]
			values <- abs(values)
			result[["pvalue"]] <- mean(values[1] <= values)
		}
	}

	if (is.na(result[["pvalue"]])) {
		result[["error"]] <- "test failed"
	}
	return(result)
}

#' plotFisherTest
#' 
#' Conduct a Fisher's exact test and plot the results as a heatmap
#' @param x    factor object or one that can be coerced to one. Alternative a 2x2 contingency matrix
#' @param y    factor object or one that can be coerced to one
#' @param name.x optional character string specifying the name for the first grouping
#' @param name.y optional character string specifying the name for the second grouping
#' @param ...  arguments passed on to \code{fisher.test}
#' @return an \code{S3} object containing the test result object as returned by \code{fisher.test} and a \code{ggplot} object
#' @author Fabian Mueller
#' @export
plotFisherTest <- function(x, y=NULL, name.x=NULL, name.y=NULL, ...){
	testRes <- fisher.test(x=x, y=y, ...)

	if (is.matrix(x) && nrow(x)==2 && ncol(x)==2){
		ct <- as.table(x)
		if (is.null(name.x)) name.x <- names(dimnames(conMat))[1]
		if (is.null(name.y)) name.y <- names(dimnames(conMat))[2]
	} else {
		ct <- table(x, y)
		if (is.null(name.x)) name.x <- deparse(substitute(x))
		if (is.null(name.y)) name.y <- deparse(substitute(y))
	}

	# pp <- ggplot2.heatmap(ct) + xlab(name.x) + ylab(name.y)
	df2p <- data.frame(ct, check.names=FALSE)
	colnames(df2p) <- c("var1", "var2", "count")

	pp <- ggplot(df2p) +  aes(var2, var1) + geom_tile(aes(fill=count)) + scale_y_discrete(limits = rev(levels(df2p[,2]))) + coord_fixed() +
		  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + geom_text(aes(label=count)) + xlab(name.y) + ylab(name.x) +
		  labs(subtitle=paste0("OR=", round(testRes$estimate, 4), " (Alternative: '", testRes$alternative,"'; p-value=", round(testRes$p.value, 4), ")"))
	res <- list(
		testRes=testRes,
		plot=pp
	)
	
	class(res) <- "FisherTestPlot"
	return(res)
}

#' rowTtest
#'
#' performs a two-sided Welch's t-test (unequal variances, equal or unequal sample sizes) on each row of a matrix X with the indices inds.1 vs indices idx2 as group assignments.
#' @author Fabian Mueller
#' @param X Matrix on which the test is performed for every row
#' @param idx1 column indices of group 1 members
#' @param idx2 column indices of group 2 members
#' @param na.rm Should NAs be removed (logical)
#' @param alternative Testing alternative. Must be one of "two.sided" (default),"less","greater" or "all".
#' 		  in case of "all" a data frome with corresping alternative variables is returned. 
#' 		  Otherwise the result is a vector.
#' @return vector (or data.frame if alternative=="all") of p-values resulting from the Welch's t-test
#' @export
#' @note Requires \code{matrixStats} package
rowTtest <- function(X, idx1, idx2=-idx1, na.rm=FALSE, alternative="two.sided"){
	require(matrixStats)
	if (!(alternative %in% c("two.sided","less","greater","all"))) {
		stop("invalid value for testing alternative")
	}
	X.1 <- X[,idx1]
	X.2 <- X[,idx2]
	# sparse matrices not fully supported yet: convert to dense matrix
	if (!is.matrix(X.1)) X.1 <- as.matrix(X.1)
	if (!is.matrix(X.2)) X.2 <- as.matrix(X.2)
		
	if (na.rm){
		n.1 <- rowSums(!is.na(X.1), FALSE)
		n.2 <- rowSums(!is.na(X.2), FALSE)
	} else {
		n.1 <- length(idx1)
		n.2 <- length(idx2)
	}

	rm.1 <- rowMeans(X.1, na.rm = na.rm)
	rm.2 <- rowMeans(X.2, na.rm = na.rm)
	rv.1 <- rowVars(X.1, na.rm = na.rm)
	rv.2 <- rowVars(X.2, na.rm = na.rm)
	rq.1 <- rv.1/n.1
	rq.2 <- rv.2/n.2
	t.stat <- (rm.1 - rm.2)/sqrt(rq.1 + rq.2)
	rdf <- (rq.1 + rq.2)^2/(rq.1^2/(n.1-1) + rq.2^2/(n.2-1)) #degrees of freedom
	rp <- rep(NA,nrow(X))
	if (alternative == "two.sided" || alternative == "all") {
		rp.2s <- 2*pt(-abs(t.stat),rdf)
	}
	if (alternative == "less" || alternative == "all") {
		rp.l <- pt(t.stat,rdf)
	}
	if (alternative == "greater" || alternative == "all") {
		rp.g <- pt(t.stat,rdf,lower.tail=FALSE)
	}
	if (alternative == "two.sided") rp <- rp.2s
	if (alternative == "greater")   rp <- rp.g
	if (alternative == "less")      rp <- rp.l
	if (alternative == "all")		rp <- data.frame(less=rp.l,greater=rp.g,two.sided=rp.2s)

	res <- data.frame(
	  pval = rp, 
	  mean1 = rm.1, 
	  mean2 = rm.2, 
	  var1 = rv.1,
	  var2 = rv.2,
	  n1 = n.1,
	  n2 = n.2
	)
	return(res)
}
