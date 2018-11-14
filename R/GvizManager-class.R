################################################################################
# GvizManager
# a class for easier managing of Gviz browser plots
################################################################################
#' @include txdb.R
NULL
require(GenomicRanges)

#' GvizManager
#'
#' A class for storing information on pipeline jobs.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{assembly}}{
#'		Genome assembly as string
#'   }
#'   \item{\code{txdb}}{
#'		\code{TxDb} object for the gene model
#'   }
#'   \item{\code{geneAnnot}}{
#'		\code{GRanges} object with gene annotation (such as gene identifiers and symbols)
#'   }
#' }
#'
#' @section Methods:
#' \describe{
#'    \item{\code{\link{firstMethod,GvizManager-method}}}{
#'      Description of the first method
#'    }
#' }
#'
#' @name GvizManager-class
#' @rdname GvizManager-class
#' @author Fabian Mueller
#' @exportClass GvizManager
setClass("GvizManager",
	slots = list(
		assembly			= "character",
		geneModelName		= "character",
		geneAnnot			= "GRanges",
		iTrackL             = "list",
		geneTrackL          = "list"
	),
	package = "muRtools"
)
setMethod("initialize", "GvizManager",
	function(
		.Object,
		assembly,
		geneModel,
		txdb,
		geneAnnot
	) {
		require(Gviz)

		.Object@assembly			<- assembly
		.Object@geneModelName	    <- geneModel
		.Object@geneAnnot			<- geneAnnot

		chroms <- names(getSeqlengths4assembly(assembly, onlyMainChrs=TRUE, adjChrNames=TRUE))
		.Object@iTrackL <- lapply(chroms, function(cc){IdeogramTrack(genome=assembly, chromosome=cc)})
		names(.Object@iTrackL) <- chroms

		.Object@geneTrackL <- lapply(chroms, function(cc){
			rr <- GeneRegionTrack(txdb, chromosome=cc, name=geneModel, transcriptAnnotation="symbol")
			symbol(rr) <- elementMetadata(geneAnnot)[gene(rr), ".symbol"]
			return(rr)
		})
		names(.Object@geneTrackL) <- chroms

		.Object
	}
)

#' @param assembly		Parameter description
#' @noRd
#' @export
GvizManager <- function(assembly, geneModel=NULL){
	geneModelMap <- c(
		"hg38"="gencode.v27",
		"hg19"="gencode.v19",
		"mm10"="gencode.vM16",
		"mm9"="gencode.vM1"
	)
	if (is.null(geneModel)) geneModel <- geneModelMap[assembly]
	txdb <- getTxDb.gencode(geneModel)
	geneAnnot <- genes(txdb)
	elementMetadata(geneAnnot)[, ".id_ensembl_long"] <- elementMetadata(geneAnnot)[, "gene_id"]
	elementMetadata(geneAnnot)[, ".id_ensembl"] <- gsub("\\.[0-9]+$", "", elementMetadata(geneAnnot)[, "gene_id"])
	elementMetadata(geneAnnot)[, ".symbol"] <- getGeneAnnotMap(assembly, from="ENSEMBL", to="SYMBOL", multiMap="paste")[elementMetadata(geneAnnot)[, ".id_ensembl"]]
	obj <- new("GvizManager",
		assembly,
		geneModel,
		txdb,
		geneAnnot
	)
	return(obj)
}

################################################################################
# Getters
################################################################################
if (!isGeneric("getIdeogramTrack")) {
	setGeneric(
		"getIdeogramTrack",
		function(object, ...) standardGeneric("getIdeogramTrack"),
		signature=c("object")
	)
}
#' getIdeogramTrack-methods
#'
#' Get a \code{Gviz} ideogram track for a specified chromosome
#'
#' @param object	\code{\linkS4class{GvizManager}} object
#' @param chrom		name of the chromosome for which the track should be returned
#' @return \code{Gviz::IdeogramTrack} for the specified chromosome
#'
#' @rdname getIdeogramTrack-GvizManager-method
#' @docType methods
#' @aliases getIdeogramTrack
#' @aliases getIdeogramTrack,GvizManager-method
#' @author Fabian Mueller
#' @export
setMethod("getIdeogramTrack",
	signature(
		object="GvizManager"
	),
	function(
		object,
		chrom
	) {
		return(object@iTrackL[[chrom]])
	}
)

################################################################################

if (!isGeneric("getGeneTrack")) {
	setGeneric(
		"getGeneTrack",
		function(object, ...) standardGeneric("getGeneTrack"),
		signature=c("object")
	)
}
#' getGeneTrack-methods
#'
#' Get a \code{Gviz} gene annotation track for a specified chromosome
#'
#' @param object	\code{\linkS4class{GvizManager}} object
#' @param chrom		name of the chromosome for which the track should be returned
#' @return \code{Gviz::GeneRegionTrack} for the specified chromosome
#'
#' @rdname getGeneTrack-GvizManager-method
#' @docType methods
#' @aliases getGeneTrack
#' @aliases getGeneTrack,GvizManager-method
#' @author Fabian Mueller
#' @export
setMethod("getGeneTrack",
	signature(
		object="GvizManager"
	),
	function(
		object,
		chrom
	) {
		return(object@geneTrackL[[chrom]])
	}
)

################################################################################

if (!isGeneric("getGeneRegionBySymbol")) {
	setGeneric(
		"getGeneRegionBySymbol",
		function(object, ...) standardGeneric("getGeneRegionBySymbol"),
		signature=c("object")
	)
}
#' getGeneRegionBySymbol-methods
#'
#' Retrieve the coordinates of a gene given the gene symbol
#'
#' @param object	\code{\linkS4class{GvizManager}} object
#' @param symbol	character specifying the gene symbol for which coordinates should be retrieved
#' @param offsetUp  offset for retrieving a flanking region upstream of the gene. Can be an \code{integer} (e.g. \code{5L}) for the absolute number of bases
#'                  or a \code{double} that specifies the relative length of the gene
#' @param offsetDown offset for retrieving a flanking region downstream of the gene. Can be an \code{integer} (e.g. \code{5L}) for the absolute number of bases
#'                  or a \code{double} that specifies the relative length of the gene
#' @param exactMatch should the symbol match exactly or should it just be contained in the symbol annotation
#' @param ignore.case ignore the case when matching the symbol
#' @return a \code{GRanges} object containing the gene coordinates
#'
#' @rdname getGeneRegionBySymbol-GvizManager-method
#' @docType methods
#' @aliases getGeneRegionBySymbol
#' @aliases getGeneRegionBySymbol,GvizManager-method
#' @author Fabian Mueller
#' @export
setMethod("getGeneRegionBySymbol",
	signature(
		object="GvizManager"
	),
	function(
		object,
		symbol,
		offsetUp=0.1,
		offsetDown=0.1,
		exactMatch=TRUE,
		ignore.case=FALSE
	) {
		
		if (exactMatch){
			if (ignore.case) {
				idx <- which(toupper(elementMetadata(object@geneAnnot)[,".symbol"])==toupper(symbol))
			} else {
				idx <- which(elementMetadata(object@geneAnnot)[,".symbol"]==symbol)
			}
		} else {
			idx <- grep(symbol, elementMetadata(object@geneAnnot)[,".symbol"], ignore.case=ignore.case)
		}
		gr <- object@geneAnnot[idx]

		if (is.double(offsetUp) && offsetUp >= 0){
			offsetUp <- ceiling(width(gr)*offsetUp)
		} else if (!is.integer(offsetUp) || offsetUp < 0){
			stop("Invalid parameter: offsetUp")
		}
		if (is.double(offsetDown) && offsetDown >= 0){
			offsetDown <- ceiling(width(gr)*offsetDown)
		} else if (!is.integer(offsetDown) || offsetDown < 0){
			stop("Invalid parameter: offsetDown")
		}
		gr <- resize(gr, width=width(gr)+offsetDown, fix="start")
		gr <- resize(gr, width=width(gr)+offsetUp, fix="end")
		return(gr)
	}
)

