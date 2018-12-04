########################################################################################################################
## lolaUtils.R
## created: 2017-07-06
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for LOLA enrichment analysis
########################################################################################################################


################################################################################
# LOLA enrichment analyses
################################################################################
#' downloadLolaDbs
#' 
#' Downloading prepared LOLA DBs from server
#' @param dest    destination directory
#' @param dbs     vector of names of LOLA DBs to be downloaded. Currently 'LOLACore' and 'LOLAExt'
#'                are supported
#' @details
#' Requires a stable internet connection. Could take a while depending on the size of the database and the internet connection
#' @return a list containing vectors of directory names for each available genome assembly
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' }
downloadLolaDbs <- function(dest, dbs=c("LOLACore")){
	res <- list()
	for (dbName in dbs){
		dbDownloadLink <- NULL
		if (dbName == "LOLACore"){
			# dbDownloadLink <- "http://regiondb.s3.amazonaws.com/LOLACoreCaches_170206.tgz"
			dbDownloadLink <- "http://cloud.databio.org.s3.amazonaws.com/regiondb/LOLACoreCaches_170206.tgz"
		} else if (dbName == "LOLAExt"){
			# dbDownloadLink <- "http://regiondb.s3.amazonaws.com/LOLAExtCaches_170206.tgz"
			dbDownloadLink <- "http://cloud.databio.org.s3.amazonaws.com/regiondb/LOLAExtCaches_170206.tgz"
		}
		if (is.null(dbDownloadLink)){
			logger.error(c("Invalid DB name:", dbName))
		}

		dbDir <- file.path(dest, dbName)
		if (dir.exists(dbDir)){
			logger.warning(c("Directory", dbDir, "already existing --> skipping download"))
		} else {
			tmpFn <- tempfile(fileext=".tgz")
			logger.start(c("Downloading from", dbDownloadLink))
			download.file(dbDownloadLink, destfile=tmpFn, mode = "wb")
			exDir <- tempfile()
			untar(tmpFn, exdir=exDir)
			baseDir <- grep(paste0(dbName, "$"), list.dirs(exDir, recursive=TRUE), value=TRUE)
			if (length(baseDir) < 1) logger.error(c("Downloaded archive does not contain a subdirectory", dbName))
			dir.create(dbDir)
			file.copy(baseDir, dest, recursive=TRUE)
			unlink(c(tmpFn, exDir), recursive=TRUE) # remove temp files
			logger.completed()
		}
		assemblies <- list.dirs(dbDir, full.names=FALSE, recursive=FALSE)
		for (aa in assemblies){
			res[[aa]] <- c(res[[aa]], file.path(dbDir, aa))
		}
	}
	return(res)
}

#' loadLolaDbs
#' 
#' Load LOLA databases from disk and merge them
#' @param lolaDbPaths  vector of names of LOLA DB paths to be loaded
#' 
#' @return LOLA DB list as returned by \code{LOLA::loadRegionDB}
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' }
loadLolaDbs <- function(lolaDbPaths){
	require("data.table") #explicitely load data.table to adress LOLA namespace issues
	require("LOLA")
	require("simpleCache") # TODO: include requirement in dependencies
	if (length(lolaDbPaths) < 1) logger.error(c("No LOLA DB paths specified"))
	logger.start("Loading LOLA DBs")
		lolaDb <- loadRegionDB(lolaDbPaths[1])
		if (length(lolaDbPaths)>1){
			for (i in 2:length(lolaDbPaths)){
				lolaDb <- mergeRegionDBs(lolaDb, loadRegionDB(lolaDbPaths[i]))
			}
		}
	logger.completed()
	return(lolaDb)
}

#' getRoadmapMdFromId
#'
#' Helper function to retrieve metadata from Roadmap sample IDs
#'
#' @param sampleIds	vector of sample ids in Roadmap fort (e.g. "E014")
#' @return a atble containing metadata for the supplied sample ids
#'
#' @author Fabian Mueller
#' @noRd
getRoadmapMdFromId <- function(sampleIds=NULL){
	remcAnnotFn <- system.file(file.path("extdata", "remc_metadata_2013.tsv"), package = "RnBeads")
	eidColname <- make.names("Epigenome ID (EID)")
	nameColname <- make.names("Epigenome Mnemonic")
	annot <- read.table(remcAnnotFn, header=TRUE, sep="\t", quote="", skip=0, comment.char="", stringsAsFactors=FALSE)
	annot <- annot[!is.na(annot[,eidColname]) & nchar(annot[,eidColname])>0, ]
	rownames(annot) <- annot[,eidColname]
	if (length(sampleIds)<1) sampleIds <- annot[,eidColname]

	snames <- annot[sampleIds,nameColname]
	snames[is.na(snames)] <- sampleIds[is.na(snames)]
	res <- data.frame(sampleId=sampleIds, sampleName=snames, stringsAsFactors=FALSE)
	return(res)
}

#' getCellTypesFromLolaDb
#'
#' retrieve or guess cell types from a LOLA DB object
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @return character vector with cell types
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getCellTypesFromLolaDb(lolaDb)
#' }
getCellTypesFromLolaDb <- function(lolaDb){
	tt <- data.frame(lolaDb$regionAnno)
	res <- tt$cellType

	# special treatment for "sheffield_dnase" from LOLACore
	isInCollection <- tt$collection=="sheffield_dnase"
	if (sum(isInCollection) > 0) res[isInCollection] <- tt$description[isInCollection]

	# special treatment for "roadmap_epigenomics" from LOLAExt
	isInCollection <- tt$collection=="roadmap_epigenomics"
	if (sum(isInCollection) > 0){
		sampleIds <- gsub("^(.+)-(.+)$", "\\1", tt$filename[isInCollection])
		res[isInCollection] <- paste0(getRoadmapMdFromId(sampleIds)[,"sampleName"], " (", sampleIds,")")
	}

	return(res)
}

#' getTargetFromLolaDb
#'
#' retrieve or guess the target from a LOLA DB object. Here, target typically
#' refers to antibodies for ChIP-seq experiments, but could also refer to other annotations
#' (e.g. motifs in TF motif databases, annotation according to UCSC features etc.)
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @return character vector with targets
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getTargetFromLolaDb(lolaDb)
#' }
getTargetFromLolaDb <- function(lolaDb){
	tt <- data.frame(lolaDb$regionAnno)
	res <- tt$antibody

	# special treatment for "ucsc_features" from LOLACore
	isInCollection <- tt$collection=="ucsc_features"
	if (sum(isInCollection) > 0) res[isInCollection] <- gsub("^UCSC ", "", tt$description[isInCollection])

	# special treatment for "sheffield_dnase" from LOLACore
	isInCollection <- tt$collection=="sheffield_dnase"
	if (sum(isInCollection) > 0) res[isInCollection] <- "DNase"

	# special treatment for "encode_segmentation" from LOLACore
	isInCollection <- tt$collection=="encode_segmentation"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub(" Segments$", "", tt$description[isInCollection])
	}

	# special treatment for "jaspar_motifs" from LOLAExt
	isInCollection <- tt$collection=="jaspar_motifs"
	if (sum(isInCollection) > 0) res[isInCollection] <- gsub("\\.bed$", "", tt$filename[isInCollection])

	# special treatment for "roadmap_epigenomics" from LOLAExt
	isInCollection <- tt$collection=="roadmap_epigenomics"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub("(\\.(hotspot\\.(fdr0\\.01|all)|macs2))?\\.(peaks\\.bed|narrowPeak)$", "", 
			gsub("^(.+)-(.+)$", "\\2", tt$filename[isInCollection])
		)
	}

	# special treatment for "encodeTFBSmm10" from LOLAExt which don't have a good description
	isInCollection <- tt$collection=="encodeTFBSmm10" & tt$description=="encodeTFBSmm10"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub("^wgEncodePsuTfbs", "", 
			gsub("RepPeaksRep[1-9]\\.broadPeak$", "", tt$filename[isInCollection])
		)
	}

	return(res)
}

#' getNamesFromLolaDb
#'
#' get human readable names from a LOLA DB object
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param addCollectionNames  attach the name of the collection to the name
#' @param addDbId  attach the index of the item in the LOLA DB object to the name
#' @return character vector with human readable names
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getNamesFromLolaDb(lolaDb)
#' }
getNamesFromLolaDb <- function(lolaDb, addCollectionNames=FALSE, addDbId=TRUE){
	tt <- data.frame(lolaDb$regionAnno)
	res <- paste0(tt$description)

	cellTypes <- getCellTypesFromLolaDb(lolaDb)
	targets <- getTargetFromLolaDb(lolaDb)

	# if there is a valid "target" (see getTargetFromLolaDb) inferred from the annotation, construct a name from it
	hasDesc <- !is.na(targets)
	if (sum(hasDesc) > 0) res[hasDesc] <- targets[hasDesc]

	# if both cellType and target annotation were inferred, construct a name from them
	hasDesc <- !is.na(cellTypes) & !is.na(targets)
	if (sum(hasDesc) > 0) res[hasDesc] <- paste0(targets[hasDesc], " - ", cellTypes[hasDesc])

	# add suffixes if requested
	if (addCollectionNames){
		res <- paste0(res, " [", tt$collection, "]")
	}
	if (addDbId){
		res <- paste0(res, " [db", 1:nrow(tt), "]")
	}
	return(res)
}

#' lolaPrepareDataFrameForPlot
#'
#' Helper function for preparing a data.frame for multiple LOLA plots
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}.
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut p-value cutoff to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \code{\link{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @noRd
lolaPrepareDataFrameForPlot <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, signifCol="qValue", includedCollections=c(), pvalCut=0.01, maxTerms=50, perUserSet=FALSE, groupByCollection=TRUE, orderDecreasing=NULL){
	#dedect by column name whether decreasing order needs to be used
	if (is.null(orderDecreasing)){
		oset.dec <- c("pValueLog", "qValueLog", "logOddsRatio", "oddsRatio", "log2OR")
		oset.inc <- c("maxRnk", "meanRnk", "qValue")
		if (!is.element(orderCol, c(oset.dec, oset.inc))){
			logger.error(c("Could not determine whether to use increasing or decreasing order for column:", orderCol))
		}
		orderDecreasing <- is.element(orderCol, oset.dec)
	}

	if (!is.element(signifCol, c("pValueLog", "qValue"))){
		logger.error(c("Invalid significance column name:", signifCol))
	}

	lolaRes <- as.data.frame(lolaRes)

	if (signifCol == "qValue"){
		lolaRes[["qValueLog"]] <- -log10(lolaRes[["qValue"]])
		signifCol <- "qValueLog"
	}

	#check if column name is in the LOLA results
	# take older versions of LOLA into account (oddsRatio was named logOddsRatio before)
	isOldLolaObj <-is.element("logOddsRatio", colnames(lolaRes))
	if (isOldLolaObj){
		logger.warning("Detected old version of LOLA. Renaming 'logOddsRatio' to 'oddsRatio'")
		colnames(lolaRes)[colnames(lolaRes)=="logOddsRatio"] <- "oddsRatio"
	}
	if (scoreCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		scoreCol <- "oddsRatio"
	}
	if (orderCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		orderCol <- "oddsRatio"
	}

	if (scoreCol=="log2OR" || orderCol=="log2OR") {
		if (isOldLolaObj){
			logger.error("I don't want to take the log of logs")
		}
		lolaRes[["log2OR"]] <- log2(lolaRes[["oddsRatio"]])
	}

	lolaRes$name <- getNamesFromLolaDb(lolaDb, addCollectionNames=!groupByCollection)[lolaRes$dbSet]
	lolaRes$name <- strTrim(lolaRes$name, len.out=50, trim.str="...", len.pref=30)
	lolaRes$target <- getTargetFromLolaDb(lolaDb)[lolaRes$dbSet]

	if (!is.element(scoreCol, colnames(lolaRes))){
		logger.error(c("Missing score column from LOLA result:", scoreCol))
	}
	if (!is.element(orderCol, colnames(lolaRes))){
		logger.error(c("Missing order column from LOLA result:", orderCol))
	}

	if (length(includedCollections) > 0){
		lolaRes <- lolaRes[lolaRes[["collection"]] %in% includedCollections,]
	}

	calledSignif <- lolaRes[[signifCol]] > -log10(pvalCut)
	lolaRes[["isSignif"]] <- calledSignif

	if (sum(calledSignif) < 1){
		# lolaRes.signif <- lolaRes
		logger.warning("No significant association found")
		return(NULL)
	} else {
		lolaRes.signif <- lolaRes[calledSignif,]
	}
	lolaRes.signif <- lolaRes.signif[order(lolaRes.signif[,orderCol], decreasing=orderDecreasing),]

	# data frame for plotting sorted according to significance/odds ratio
	df2p <- lolaRes.signif
	dbSets2keep <- unique(df2p[["dbSet"]])
	if (length(dbSets2keep) > maxTerms){
		if (perUserSet){
			logger.info(c("Reduced the number of region sets in the plot to", maxTerms, "per userSet"))
			uSetF <- lolaRes[["userSet"]]
			if (!is.factor(uSetF)) uSetF <- factor(uSetF)
			dbSets2keep <- c()
			for (ss in levels(uSetF)){
				dbSets.cur <- unique(df2p[uSetF==ss,"dbSet"])
				dbSets2keep <- union(dbSets2keep, dbSets.cur[1:min(maxTerms, length(dbSets.cur))])
			}
		} else {
			logger.info(c("Reduced the number of region sets in the plot to", maxTerms))
			dbSets2keep <- dbSets2keep[1:maxTerms]
		}
		df2p <- df2p[df2p[["dbSet"]] %in% dbSets2keep, ]
	}

	df2p[["term"]] <- factor(df2p[["name"]], levels=unique(df2p[["name"]])) #sort factor by occurence (data.frame is ordered)
	df2p[["target"]] <- factor(df2p[["target"]])

	#bound score at maximum (excluding +/-Inf)
	scoresFinite <- df2p[[scoreCol]][is.finite(df2p[[scoreCol]])]
	cutVal <- max(scoresFinite, na.rm=TRUE)
	df2p[[scoreCol]][df2p[[scoreCol]]>cutVal] <- cutVal
	cutVal <- min(scoresFinite, na.rm=TRUE)
	df2p[[scoreCol]][df2p[[scoreCol]]<cutVal] <- cutVal

	return(df2p)
}


#' lolaVolcanoPlot
#'
#' plot a volcano plot showing LOLA enrichment results: LOLA p-value against the log-odds score. Colored by rank
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}.
#' @param colorBy  annotation/column in the the LOLA DB that should be used for point coloring
#' @param colorpanel colors to be used for coloring the points
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # example taken from RnBeads
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' # compute differential methylation
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' # perform enrichment analysis
#' res <- performLolaEnrichment.diffMeth(rnb.set.example,dm,lolaDirs[["hg19"]])
#' # select the 500 most hypermethylated tiling regions in ESCs compared to iPSCs
#' # in the example dataset
#' lolaRes <- res$region[["hESC vs. hiPSC (based on Sample_Group)"]][["tiling"]]
#' lolaRes <- lolaRes[lolaRes$userSet=="rankCut_500_hyper",]
#' # plot
#' lolaVolcanoPlot(res$lolaDb, lolaRes, signifCol="qValue")
#' }
lolaVolcanoPlot <- function(lolaDb, lolaRes, includedCollections=c(), signifCol="qValue", colorBy="maxRnk", colorpanel=c()){
	if (length(unique(lolaRes[["userSet"]])) > 1){
		logger.warning("Multiple userSets contained in LOLA result object")
	}

	#prepare data.frame for plotting
	# adjust maxTerms, pvalCut to not filter anything
	df2p <- lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol="pValueLog", signifCol=signifCol, orderCol="maxRnk", includedCollections=includedCollections, pvalCut=1.1, maxTerms=Inf, perUserSet=FALSE, groupByCollection=TRUE, orderDecreasing=NULL)
	#reverse the order of points s.t. the plot looks nice
	df2p <- df2p[nrow(df2p):1,]

	if (signifCol == "qValue"){
		signifCol <- "qValueLog"
	}

	is.color.gradient <- FALSE
	is.color.discrete <- FALSE
	if (!is.null(colorBy) && nchar(colorBy) > 0 && !is.element(colorBy, colnames(df2p))){
		logger.warning(c("Invalid colorBy argument:", colorBy, "--> using black as color"))
	} else {
		is.color.gradient <- is.numeric(df2p[[colorBy]])
		is.color.discrete <- !is.color.gradient
	}

	oddsRatioCol <- "oddsRatio"
	if (!is.element(oddsRatioCol, colnames(df2p))){
		logger.error("Invalid LOLA result. Could not find a valid column containing odds ratios.")
	}
	pp <- ggplot(df2p) + aes_string(oddsRatioCol, signifCol, color=colorBy) + geom_point()
	cpanel <- colorpanel
	if (length(cpanel) < 1){
		if (is.color.gradient){
			#default: viridis 
			cpanel <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")
		}
		if (is.color.discrete){
			if (is.factor(df2p[[colorBy]])){
				cpanel <- rep(colpal.mu.cat, length.out=nlevels(df2p[[colorBy]]))
				names(cpanel) <- levels(df2p[[colorBy]])
			} else if (is.character(df2p[[colorBy]])){
				cpanel.names <- unique(df2p[[colorBy]])
				cpanel <- rep(colpal.mu.cat, length.out=length(cpanel.names))
				names(cpanel) <- cpanel.names
			} else {
				logger.error("invalid discrete coloring column type")
			}
		}
	}
	if (is.color.gradient){
		pp <- pp + scale_color_gradientn(colors=cpanel)
	}
	if (is.color.discrete){
		pp <- pp + scale_color_manual(na.value="#C0C0C0", values=cpanel)
	}
		
	return(pp)
}

#' lolaBarPlot
#'
#' plot a barplot of LOLA enrichment results
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut p-value cutoff to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \code{\link{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # example taken from RnBeads
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' # compute differential methylation
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' # perform enrichment analysis
#' res <- performLolaEnrichment.diffMeth(rnb.set.example,dm,lolaDirs[["hg19"]])
#' # select the 500 most hypermethylated tiling regions in ESCs compared to iPSCs
#' # in the example dataset
#' lolaRes <- res$region[["hESC vs. hiPSC (based on Sample_Group)"]][["tiling"]]
#' lolaRes <- lolaRes[lolaRes$userSet=="rankCut_500_hyper",]
#' # plot
#' lolaBarPlot(res$lolaDb, lolaRes, scoreCol="oddsRatio", orderCol="maxRnk", pvalCut=0.05)
#' }
lolaBarPlot <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, signifCol="qValue", includedCollections=c(), pvalCut=0.01, maxTerms=50, colorpanel=sample(rainbow(maxTerms,v=0.5)), groupByCollection=TRUE, orderDecreasing=NULL){
	if (length(unique(lolaRes[["userSet"]])) > 1){
		logger.warning("Multiple userSets contained in LOLA result object")
	}
	if (scoreCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		scoreCol <- "oddsRatio"
	}
	#prepare data.frame for plotting
	df2p <- lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol=scoreCol, orderCol=orderCol, signifCol=signifCol, includedCollections=includedCollections, pvalCut=pvalCut, maxTerms=maxTerms, perUserSet=FALSE, groupByCollection=groupByCollection, orderDecreasing=orderDecreasing)

	if (is.null(df2p)){
		return(ggMessagePlot("No significant association found"))
	}

	aesObj <- aes_string("term", scoreCol)
	if (length(colorpanel) > 0) aesObj <- aes_string("term", scoreCol, fill="target")

	pp <- ggplot(df2p) + aesObj + geom_col() + 
		  scale_x_discrete(name="") + 
		  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

	if (length(colorpanel) > 0){
		# check if the defined color panel already contains a mapping of targets to colors
		# otherwise create an arbitrary mapping
		if (!is.null(names(colorpanel)) && all(levels(df2p[["target"]]) %in% names(colorpanel))){
			cpanel <- colorpanel
		} else {
			cpanel <- rep(colorpanel, length.out=nlevels(df2p[["target"]]))
			names(cpanel) <- levels(df2p[["target"]])
		}
		pp <- pp + scale_fill_manual(na.value="#C0C0C0", values=cpanel, guide=FALSE)
	}

	if (groupByCollection){
		pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
	}
	return(pp)
}

#' lolaBoxPlotPerTarget
#'
#' plot a boxplot showing LOLA enrichment results per "target" group (see \code{\link{getTargetFromLolaDb}} for an explanation of
#' "target").
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut  p-value cutoff to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \code{\link{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @param scoreDecreasing flag indicating whether the value in \code{scoreCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # example taken from RnBeads
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' # compute differential methylation
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' # perform enrichment analysis
#' res <- performLolaEnrichment.diffMeth(rnb.set.example,dm,lolaDirs[["hg19"]])
#' # select the 500 most hypermethylated tiling regions in ESCs compared to iPSCs
#' # in the example dataset
#' lolaRes <- res$region[["hESC vs. hiPSC (based on Sample_Group)"]][["tiling"]]
#' lolaRes <- lolaRes[lolaRes$userSet=="rankCut_500_hyper",]
#' # plot
#' lolaBoxPlotPerTarget(res$lolaDb, lolaRes, scoreCol="oddsRatio", orderCol="maxRnk", pvalCut=0.05)
#' }
lolaBoxPlotPerTarget <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, signifCol="qValue", includedCollections=c(), pvalCut=0.01, maxTerms=50, colorpanel=c(), groupByCollection=TRUE, orderDecreasing=NULL, scoreDecreasing=NULL){
	if (length(unique(lolaRes[["userSet"]])) > 1){
		logger.warning("Multiple userSets contained in LOLA result object")
	}
	# detect by column name whether decreasing order needs to be used for the score column
	if (is.null(scoreDecreasing)){
		oset.dec <- c("pValueLog", "logOddsRatio", "oddsRatio")
		oset.inc <- c("maxRnk", "meanRnk")
		if (!is.element(scoreCol, c(oset.dec, oset.inc))){
			logger.error(c("Could not determine whether to use increasing or decreasing order for score column:", scoreCol))
		}
		scoreDecreasing <- is.element(scoreCol, oset.dec)
	}
	if (scoreCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		scoreCol <- "oddsRatio"
	}
	#prepare data.frame for plotting
	df2p <- lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol=scoreCol, orderCol=orderCol, signifCol=signifCol, includedCollections=includedCollections, pvalCut=pvalCut, maxTerms=Inf, perUserSet=FALSE, groupByCollection=groupByCollection, orderDecreasing=orderDecreasing)

	if (is.null(df2p)){
		return(ggMessagePlot("No significant association found"))
	}

	# maxTerms now applies to the targets, not to LOLA DB items
	if (nlevels(df2p$target) > maxTerms){
		# if there are too many targets, select the targets from the top of the ordered table
		targets2keep <- unique(as.character(df2p$target))[1:maxTerms]
		df2p <- df2p[df2p$target %in% targets2keep,]
		df2p[["target"]] <- factor(as.character(df2p[["target"]]))
		logger.info(c("Reduced the number of targets in the plot to", maxTerms))
	}

	# sort targets by median score
	targetScore <- tapply(df2p[[scoreCol]], df2p[["target"]], FUN=function(x){median(x, na.rm=TRUE)})
	df2p[["target"]] <- factor(as.character(df2p[["target"]]),
		levels=names(targetScore)[order(targetScore, decreasing=scoreDecreasing)]
	)

	aesObj <- aes_string("target", scoreCol)
	if (length(colorpanel) > 0) aesObj <- aes_string("target", scoreCol, fill="target")

	pp <- ggplot(df2p) + aesObj + geom_boxplot(outlier.shape=NA) +
		  geom_dotplot(aes(fill=NULL), binaxis="y", stackdir="center") +
		  scale_x_discrete(name="") + 
		  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

	if (length(colorpanel) > 0){
		# check if the defined color panel already contains a mapping of targets to colors
		# otherwise create an arbitrary mapping
		if (!is.null(names(colorpanel)) && all(levels(df2p[["target"]]) %in% names(colorpanel))){
			cpanel <- colorpanel
		} else {
			cpanel <- rep(colorpanel, length.out=nlevels(df2p[["target"]]))
			names(cpanel) <- levels(df2p[["target"]])
		}
		pp <- pp + scale_fill_manual(na.value="#C0C0C0", values=cpanel, guide=FALSE)
	} 
	if (groupByCollection){
		pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
	}
	return(pp)
}

#' lolaRegionSetHeatmap
#' Plot a heatmap in which the rows are different user sets and the color corresponds to an enrichment score
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \link{\code{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}
#' @param markSignif mark significant enrichments in the heatmap
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut p-value cutoff (negative log10) to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colorpanel for the heatmap gradient
#' @param groupByCollection facet the plot by collection
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @export
lolaRegionSetHeatmap <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, signifCol="qValue", markSignif=FALSE, includedCollections=c(), pvalCut=0.01, maxTerms=50, colorpanel=colpal.cont(9, "cb.OrRd"), groupByCollection=TRUE, orderDecreasing=NULL){
	#prepare data.frame for plotting
	df2p <- muRtools:::lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol=scoreCol, orderCol=orderCol, signifCol=signifCol, includedCollections=includedCollections, pvalCut=pvalCut, maxTerms=maxTerms, perUserSet=TRUE, groupByCollection=groupByCollection, orderDecreasing=orderDecreasing)

	if (is.null(df2p)){
		return(rnb.message.plot("No significant association found"))
	}

	if (!is.factor(df2p$userSet)) df2p$userSet <- factor(df2p$userSet)
	if (markSignif) df2p$signifTxt <- ifelse(df2p$isSignif, "*", NA) 

	pp <- ggplot(df2p) + aes(term, userSet) + geom_tile(aes_string(fill=scoreCol)) + 
		  scale_x_discrete(name="") + scale_y_discrete(limits=rev(levels(df2p$userSet)), name="") + 
		  scale_fill_gradientn(colors=colorpanel)
	if (markSignif){
		pp <- pp + geom_text(aes(label=signifTxt))
	}
	pp <- pp + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
	if (groupByCollection){
		pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
	}
	return(pp)
}
