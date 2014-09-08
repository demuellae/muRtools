#' getRegionSet
#' 
#' retrieves custom region sets by id
#' @param id region set id
#' @param assembly genome assembly
#' @return a GRanges object containing the region set
#' @export
#' @aliases getRegionSet
#' @examples 
#' regs <- getRegionSet("deep_chip_ctrl_regions","hg19")
getRegionSet <- function(id,assembly){
	reg.fname <- system.file(file.path("extdata", paste(id,assembly,"bed",sep=".")), package = "muRtools")
	res <- NULL
	if (!file.exists(reg.fname)){
		stop(paste0("Error in retrieving region set. Could not open region file: ",reg.fname))
	}
	#parse the bed file to GRanges
	res <- bed2GRanges(reg.fname,assembly=assembly)
	return(res)
}


#' bed2GRanges
#' 
#' loads a bed file and converts it to a GRanges object. Note that bed files are 0-based, right-exclusive by definition
#' The output of this function will be 1-based, right inclusive as defined by GRanges.
#' @param id region set id
#' @param assembly genome assembly
#' @return a GRanges object containing the region set
#' @export
#' @aliases bed2GRanges
#' @examples
#' reg.fname <- system.file(file.path("extdata","deep_chip_ctrl_regions.hg19.bed"), package = "muRtools")
#' regs <- bed2GRanges(reg.fname)
bed2GRanges <- function(fname,assembly=NA){
	require(GenomicRanges)
	BED.COLUMNS <- c("chrom", "start", "end", "id", "score", "strand")
	if (assembly=="hg19"){
		require(BSgenome.Hsapiens.UCSC.hg19)
		sls <- seqlengths(Hsapiens)
	} else if (assembly=="GRCh38") {
		require(BSgenome.Hsapiens.NCBI.GRCh38)
	} else if (assembly=="mm9"){
		require(BSgenome.Mmusculus.UCSC.mm9)
		sls <- seqlengths(Mmusculus)
	} else if (assembly=="mm10"){
		require(BSgenome.Mmusculus.UCSC.mm10)
		sls <- seqlengths(Mmusculus)
	} else {
		warning(paste0("Unknown assembly: ",assembly,". Did not incorporate seqlengths"))
		assembly <- NA
	}
	tbl <- tryCatch(suppressWarnings(read.delim(fname, header = FALSE, quote = "", comment.char = "#",
				stringsAsFactors = FALSE, na.strings = "")), error = function(e) { e })
	if (inherits(tbl, "error")) {
		if (grepl("cannot open", tbl, fixed = TRUE)) {
			stop(paste0("Could not open bed file: ",fname))
		}
		stop("invalid file format")
	}
	if (ncol(tbl) < 3) {
		stop("invalid file format; expected at least 3 columns")
	}
	if (!(is.integer(tbl[[2]]) && is.integer(tbl[[3]]))) {
		stop("invalid file format; expected start and end positions in columns 2 and 3, respectively")
	}
	bed.col.inds <- 1:min(length(BED.COLUMNS),ncol(tbl))
	colnames(tbl)[bed.col.inds] <- BED.COLUMNS[bed.col.inds]

	sns <- tbl$chrom
	if (!is.na(assembly)) {
		sns <- factor(sns,levels=names(sls))
	}
	res <- GRanges(
		seqnames=sns,
		ranges=IRanges(tbl$start + 1,end=tbl$end) #0-based, right exclusive --> 1-based, right inclusive
	)
	if (is.element("strand",colnames(tbl))) {
		strand(res) <- tbl$strand
	}
	if (is.element("score",colnames(tbl))) {
		elementMetadata(res)[["score"]] <- tbl$score
	}
	if (is.element("id",colnames(tbl))) {
		elementMetadata(res)[["id"]] <- tbl$id
		names(res) <- tbl$id
	}
	if (ncol(tbl)>6) {
		emd <- tbl[,7:ncol(tbl),drop=FALSE]
		elementMetadata(res) <- data.frame(elementMetadata(res),emd)
		names(res) <- tbl$id
	}
	if (!is.na(assembly)) {
		seqlengths(res) <- sls
	}
	return(res)
}
