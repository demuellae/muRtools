################################################################################
# Utilities for configuring analysis projects
################################################################################

################################################################################
# Parsers for config elements
################################################################################
cfg_parse_free <- function(x){
	return(x)
}

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
	if (!all(dir.exists(x))) stop("directories must exist")
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
	if (!all(file.exists(x))) stop("file(s) must exist")
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
		user <- ""
		if (is.element("USER", names(y))) user <- y[["USER"]][1]
		# print(req)
		CommandRslurm(logDir=y[["LOG_DIR"]], req=req, user=user)
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
			source(utilsFn, chdir=TRUE)
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

#' getScriptDependencyGraph
#' 
#' Get a script dependency graph by parsing analysis R scripts
#' @param scriptDirs   directories containing R scripts
#' @return \code{igraph} object containing dependencies between scripts
#' @author Fabian Mueller
#' @noRd
getScriptDependencyGraph <- function(scriptDirs="."){
	# scriptDirs <- "~/muResearchCode/BrainDev/src/"
	require(igraph)
	scrptFns <- list.files(scriptDirs, pattern="*.R$", full.names=TRUE, recursive=TRUE)
	sdf <- data.frame(
		script=basename(scrptFns),
		path=scrptFns,
		anaName=as.character(NA),
		stringsAsFactors=FALSE
	)
	sdf[is.na(sdf[,"anaName"]),"anaName"] <- sdf[is.na(sdf[,"anaName"]), "script"]

	importL <- rep(list(c()), nrow(sdf))
	for (i in 1:nrow(sdf)){
		lls <- trimws(readLines(sdf[i,"path"]))
		lls <- lls[!grepl("^#", lls)] # remove comments
		anaNameLines <- grep("^ANA_NAME ?<-", lls, value=TRUE)
		if (length(anaNameLines) > 0){
			anaNameLines <- anaNameLines[length(anaNameLines)]
			sdf[i, "anaName"] <- gsub("^.+\"(.+)\"$", "\\1", anaNameLines)
		}
		importLines <- grep('getRelatedAnaDirFromConfig\\(', lls, value=TRUE)
		if (length(importLines) > 0){
			importL[[i]] <- sort(unique(gsub("^.+getRelatedAnaDirFromConfig\\((config ?, ?)['\"]([^'\"]+)['\"][,\\)].*$", "\\2", importLines)))
		}
	}

	vvs <- unique(sort(c(sdf[,"anaName"], unlist(importL))))
	adjM <- matrix(FALSE, nrow=length(vvs), ncol=length(vvs))
	rownames(adjM) <- colnames(adjM) <- vvs
	for (i in 1:nrow(sdf)){
		if (length(importL[[i]]) > 0) adjM[importL[[i]], sdf[i,"anaName"]] <- TRUE
	}

	gg <- graph_from_adjacency_matrix(adjM, mode="directed", diag=FALSE)
	V(gg)$hasScript <- V(gg)$name %in% sdf[,"anaName"]
	V(gg)[V(gg)$hasScript]$scriptNames <- sapply(V(gg)[V(gg)$hasScript]$name, FUN=function(vn){
		paste(sdf[sdf[,"anaName"]==vn,"script"], collapse=";")
	})

	V(gg)$color <- ifelse(V(gg)$hasScript, "#a6cee3", "#bdbdbd")
	V(gg)$label.color <- "black"
	E(gg)$color <- "black"

	# plot(gg, edge.arrow.size=0.75)
	# plot(make_ego_graph(gg, order=99999, nodes=c("integrative_sc_regulation"), mode="in")[[1]], edge.arrow.size=0.75)
	return(gg)
}
