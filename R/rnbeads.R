################################################################################
# Extended RnBeads functionality
################################################################################

#' rnb.xml2options.backwardsCompatible
#' 
#' Extends \code{RnBeads::rnb.xml2options} to drop legacy options, so that old XML files can be read
#' @param fname   file name of the XML file to be read
#' @param ...     arguments passed on to
#' @return result of \code{RnBeads::rnb.xml2options} with dropped legacy options
#' @author Fabian Mueller
#' @noRd
#' 
rnb.xml2options.backwardsCompatible <- function(fname, ...){
	legacyNames <- c(
		"differential.enrichment"
	)
	
	xml.tree <- tryCatch(XML::xmlRoot(XML::xmlTreeParse(fname)), error = function(e) { return(e) })
	if (inherits(xml.tree, "error")) {
		stop(xml.tree$message)
	}

	anyRemoved <- FALSE
	for (nn in legacyNames){
		if (!is.null(xml.tree[[nn]])){
			logger.warning(c("Removed legacy option:", nn))
			xml.tree[[nn]] <- NULL
			anyRemoved <- TRUE
		}
	}
	if (anyRemoved){
		newXmlFname <- tempfile(fileext=".xml")
		saveXML(xml.tree, file=newXmlFname)
		return(rnb.xml2options(newXmlFname, ...))
	}
	return(rnb.xml2options(fname, ...))	
}

#' loadRnBeadsAnalysis
#' 
#' Loads RnBeads analysis results (RnBSet object, differential methylation) and optionally
#' runs the corresponding preanalysis script and sets options
#' @param input   input directory
#' @param type    analysis type. Currently only "cluster_run" is supported
#' @param preprocessed determines the type of RnBSet object to be loaded. If \code{TRUE} (default)
#'                the preprossed/filtered object will be loaded. Otherwise the raw, imported object
#' @param setOptions should the analysis options be set from the options settings of the RnBeads run?
#' @param preAnalysis should the corresponding preanalysis script be called
#' @param diffMeth should the differential analysis result be loaded in addition?
#' @return A \code{list} containing the \code{RnBSet} object (list item \code{rnbs}) and optionally
#'         the differential methylaiton result (list item \code{diffmeth})
#' @author Fabian Mueller
#' @export
loadRnBeadsAnalysis <- function(input, type="cluster_run", preprocessed=TRUE, setOptions=TRUE, preAnalysis=TRUE, diffMeth=TRUE){
	rnbs <- NULL
	dm <- NULL
	if (type=="cluster_run"){
		baseDir <- input
		path.rnbs <- "preprocessing_RnBSet"
		if (!preprocessed) path.rnbs <- "import_RnBSet"
		path.rnbs <- file.path(baseDir, "cluster_run", path.rnbs)
		path.dm <- file.path(baseDir, "cluster_run", "differential_rnbDiffMeth")
		path.opts <- file.path(baseDir, "cluster_run", "options.xml")
		print(path.opts)
	} else {
		logger.error("invalid type of analysis")
	}
	logger.start("Loading RnBeads Analysis")
		if (!is.null(path.opts)){
			rnb.settings <- rnb.xml2options.backwardsCompatible(path.opts, return.full.structure=TRUE)
			path.preana <- rnb.settings[["preanalysis.script"]]
			if (setOptions){
				logger.start("Setting RnBeads options")
					if (length(rnb.settings$options) != 0) {
						do.call(rnb.options, rnb.settings$options)
					}
				logger.completed()
			}
		}
		if (!is.null(path.preana) && preAnalysis){
			logger.start("Running preanalysis script")
				source(path.preana)
			logger.completed()
		}
		logger.start("Loading RnBSet")
			rnbs <- load.rnb.set(path.rnbs)
		logger.completed()
		if (!is.null(path.dm) && diffMeth){
			logger.start("Loading differential methylation")
				dm <- load.rnb.diffmeth(path.dm)
			logger.completed()
		}
	logger.completed()
	res <- list(
		rnbset=rnbs,
		diffmeth=dm
	)
	return(res)
}
