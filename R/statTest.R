################################################################################
# Statistical testing
################################################################################

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

	## Focus on common values
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
