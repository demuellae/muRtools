################################################################################
# Utilities for configuring analysis projects
################################################################################

################################################################################
# Parsers for config elements
################################################################################
cfg_parse_string <- function(x){
	if (!is.character(x)) stop("Config element is not a string")
	return(x)
}

cfg_parse_directory <- function(x){
	if (!is.character(x)) stop("Config element is not a string (directory)")
	return(x)
}

cfg_parse_directory_existing <- function(x){
	x <- cfg_parse_directory(x)
	if (!dir.exists(x)) stop("directory must exist")
	return(x)
}

cfg_parse_directory_non_existing <- function(x){
	x <- cfg_parse_directory(x)
	if (dir.exists(x)) stop("directory must be non-existing")
	return(x)
}

cfg_parse_directory_warn_existing <- function(x){
	x <- cfg_parse_directory(x)
	if (dir.exists(x)) logger.warning(c("Directory already exists:", x))
	return(x)
}

cfg_parse_directory_create <- function(x){
	x <- cfg_parse_directory_warn_existing(x)
	if (!dir.exists(x)){
		logger.info(c("Creating directory:", x))
		dir.create(x)
	}
	return(x)
}

cfg_parse_file <- function(x){
	if (!is.character(x)) stop("Config element is not a string (file)")
	return(x)
}

cfg_parse_file_existing <- function(x){
	x <- cfg_parse_file(x)
	if (!file.exists(x)) stop("file must exist")
	return(x)
}

cfg_parse_file_non_existing <- function(x){
	x <- cfg_parse_file(x)
	if (file.exists(x)) stop("file must be non-existing")
	return(x)
}

cfg_parse_file_warn_existing <- function(x){
	x <- cfg_parse_file(x)
	if (file.exists(x)) logger.warning(c("File already exists:", x))
	return(x)
}

cfg_parse_json_colormap <- function(x){
	fn <- cfg_parse_file_existing(x)
	cmaps.parsed <- fromJSON(fn)
	cmaps <- lapply(cmaps.parsed, FUN=function(x){
		mm <- unlist(x)
		names(mm) <- names(x)
		return(mm)
	})
	names(cmaps) <- names(cmaps.parsed)
	return(cmaps)
}

cfg_parse_commandr_slurm <- function(x){
	require(muPipeR)
	commandRobjs <- lapply(x, FUN=function(y){
		req <- as.character(y[["REQ"]])
		names(req) <- names(y[["REQ"]])
		print(req)
		CommandRslurm(logDir=y[["LOG_DIR"]], req=req)
	})
	names(commandRobjs) <- names(x)
	return(commandRobjs)
}
################################################################################
# Main functions
################################################################################
#' getConfig
#' 
#' Get analysis configuration object from a JSON file. Automatically runs parsers
#' for config elements specified in the '.MU_ANA_CONFIG' section of the JSON file
#' @param cfgFn   JSON file name
#' @param anaName  name of the analysis to be conducted
#' @param addDirs  should the analysis directories be created and added to the object
#' @return an S3 object containing analysis config elements
#' @author Fabian Mueller
#' @export
getConfig <- function(cfgFn, anaName, addDirs=TRUE){
	require(jsonlite)
	cfgJson <- fromJSON(cfgFn)

	cfgElems <- cfgJson[[".MU_ANA_CONFIG"]][["CONFIG_ELEMS"]]
	if (is.null(cfgElems)) logger.error("Invalid mu analysis config file")
	cfgElems.reqCols <- c("NAME", "TYPE", "OPTIONAL")
	if (!is.data.frame(cfgElems) || !all(cfgElems.reqCols %in% colnames(cfgElems))) logger.error("invalid CONFIG_ELEMS data frame")
	
	baseDir <- cfgJson[[cfgJson[[".MU_ANA_CONFIG"]][["ANALYSIS_DIR_BASE"]]]]
	if (is.null(baseDir)) logger.warning("Invalid base directory specified in JSON (ANALYSIS_DIR_BASE)")
	version <- cfgJson[[cfgJson[[".MU_ANA_CONFIG"]][["VERSION"]]]]
	if (is.null(version)) logger.warning("Invalid analysis version specified in JSON (VERSION)")
	config <- list(
		.anaName=anaName,
		.anaBaseDir=baseDir,
		.anaVersion=version,
		.fromJson=cfgJson
	)

	req.config.there <- cfgElems[,"OPTIONAL"] || cfgElems[, "NAME"] %in% names(cfgJson)
	if (!all(req.config.there)){
		logger.error(c("Missing from config file:", paste(cfgElems[!req.config.there,"NAME"], collapse=",")))
	}


	for (i in 1:nrow(cfgElems)){
		eRes <- NULL
		eRes <- tryCatch(
			do.call(paste0("cfg_parse_", cfgElems[i, "TYPE"]), list(cfgJson[[cfgElems[i, "NAME"]]])),
			error = function(ee) {
				logger.error(c("Could not parse config element: ", cfgElems[i, "NAME"], "; Reason:", ee$message))
				NULL
			}
		)
		config <- c(config, list(eRes))
		names(config)[length(config)] <- cfgElems[i, "NAME"]
	}
	if (addDirs){
		config[[".anaDir"]] <- cfg_parse_directory_create(file.path(baseDir, paste(anaName, version, sep="_")))
	}

	#source utility function if a corresponding file exists
	scrptDir <- cfgJson[[cfgJson[[".MU_ANA_CONFIG"]][["SCRIPT_DIR_BASE"]]]]
	if (!is.null(scrptDir)){
		utilsFn <- file.path(scrptDir, "utils", "utils.R")
		if (file.exists(utilsFn)){
			logger.info(c("Sourcing utils.R from", utilsFn))
			source(utilsFn)
		}
	}

	class(config) <- "muConfig"
	return(config)
}

#' getRelatedAnaDirFromConfig
#' 
#' Get analysis configuration object from a JSON file. Automatically runs parsers
#' for config elements specified in the '.MU_ANA_CONFIG' section of the JSON file
#' @param cfg   configuration object as returned by \code{getConfig} 
#' @param anaName  name of the analysis
#' @param anaVersion  version of the analysis. If \code{NULL}, it will look for the most recent one
#' @return path of the related analysis directory
#' @author Fabian Mueller
#' @export
getRelatedAnaDirFromConfig <- function(cfg, anaName, anaVersion=cfg[[".anaVersion"]]){
	if (!is.element(".anaBaseDir", names(config))) logger.error("To get related analysis directories, cfg must contain '.anaBaseDir'")
	if (is.null(anaVersion)){
		logger.info("getting most recent analysis directory according to version number")
		dds <- list.dirs(cfg[[".anaBaseDir"]], recursive=FALSE)
		dds <- grep(paste0(anaName, "_v[0-9]+$"), dds, value=TRUE)
		if (length(dds) < 1) logger.error("no related directory found")
		vers <- gsub(".*_(v[0-9]+)$", "\\1", dds)
		anaVersion <- vers[vers==max(vers)][1]
	}
	res <- file.path(cfg[[".anaBaseDir"]], paste(anaName, anaVersion, sep="_"))
	logger.info(c("Related analysis directory:", res))
	if (!dir.exists(res)){
		logger.error(c("Related analysis directory does not exist"))
	}
	return(res)
}
