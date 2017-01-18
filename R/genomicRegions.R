################################################################################
# Extensions and utitility functions for GRanges
################################################################################
#' granges2bed
#' 
#' Save a GRanges object to a bed file
#' 
#' @param gr GRanges object
#' @param fn filename to save bed file to
#' @param sc score vector or column in elementMetadata of GRanges
#' @param addAnnotCols add the columns strored in elementMetadata of GRanges
#' @param colNames add column names
#' @param doSort sort the regions before writing the output
#' @return result of writing the table (see \code{write.table})
#' @export 
granges2bed <- function(gr, fn, score=NULL, addAnnotCols=FALSE, colNames=FALSE, doSort=TRUE){
	if (doSort){
		oo <- order(as.integer(seqnames(gr)),start(gr), end(gr), as.integer(strand(gr)))
		gr <- gr[oo]
	}
	nns <- rep(".", length(gr))
	if (!is.null(names(gr))) nns <- names(gr)
	sc <- rep(".", length(gr))
	if (!is.null(score)){
		sc <- score
	}
	tt <- data.frame(
		chrom=seqnames(gr),
		start=format(start(gr)-1, trim=TRUE, scientific=FALSE),
		end=format(end(gr), trim=TRUE, scientific=FALSE),
		name=nns,
		score=sc,
		strand=strand(gr),
		stringsAsFactors=FALSE
	)
	if (addAnnotCols){
		tt <- data.frame(tt, elementMetadata(gr))
	}
	write.table(tt, file=fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=colNames)
}

#' getSeqlengths4assembly
#'
#' retrieve chromosomes/contigs and their sequence lengths for known assemblies
#'
#' @param assembly    assembly
#' @return named vector of chromosomes/contigs and sequence lengths
getSeqlengths4assembly <- function(assembly){
	res <- c()
	if (is.element(assembly, c("hg19"))){
		res <- seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
	} else if (is.element(assembly, c("hg38", "GRCh38"))){
		res <- seqlengths(BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens)
	} else if (is.element(assembly, c("hg38_chr", "GRCh38_chr"))){
		res <- seqlengths(BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens)
		prependChr <- c(1:22, "X", "Y", "MT")
		repNames <- names(res) %in% prependChr
		names(res)[repNames] <- paste0("chr", names(res)[repNames])
		names(res)[names(res)=="chrMT"] <- "chrM"
	} else if (is.element(assembly, c("mm9"))){
		res <- seqlengths(BSgenome.Mmusculus.UCSC.mm9::Mmusculus)
	} else if (is.element(assembly, c("mm10"))){
		res <- seqlengths(BSgenome.Mmusculus.UCSC.mm10::Mmusculus)
	} else {
		stop(paste0("Unknown assembly:", assembly))
	}
	return(res)
}
#' matchStrand
#'
#' match commonly used strand names to \code{"+", "-", "*"}
#'
#' @param values    character vector or factor of strand names
#' @return Factor of genomic strand (with levels \code{"+", "-", "*"})
matchStrand <- function(values) {
	if (is.factor(values) && setequal(levels(values), c("+", "-", "*"))) {
		values[is.na(values)] <- "*"
	} else {
		values <- as.character(values)
		values[is.na(values)] <- "*"
		i.positive <- values %in% c("+", "1", "+1", "F", "f", "TOP")
		i.negative <- values %in% c("-", "-1", "R", "r", "BOT")
		values[i.positive] <- "+"
		values[i.negative] <- "-"
		values[!(i.positive | i.negative)] <- "*"
		values <- factor(values, levels = c("+", "-", "*"))
	}
	return(values)
}
#' df2granges
#'
#' Converts a \code{data.frame} that defines genomic regions to object of type \code{GRanges}.
#'
#' @param df         Table defining genomic regions.
#' @param ids        Region names (identifiers) as a \code{character} vector, or \code{NULL} if no names are present.
#' @param chrom.col  Column name or index that lists the chromosome names.
#' @param start.col  Column name or index that lists the start positions of the regions.
#' @param end.col    Column name or index that lists the end positions of the regions.
#' @param strand.col Column name or index that lists the strands on which the regions are located. Set this to
#'                   \code{NULL} if this region set is not strand-specific.
#' @param coord.format Coordinate format \code{"B1RI"} for 1-based right-inclusive (default), \code{"B0RE"} for
#' 					 0-based right-exclusive.
#' @param assembly   Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported
#'                   genomes.
#' @param doSort     Should the resulting table be sorted
#' @return \code{GRanges} object encapsulating of regions included in \code{df}.
#' 	       As GRanges, the coordinates will be 1-based right-inclusive.
#'         Columns other that the ones listed as parameters in this function are included as elementMetadata.
#'
#' @export
#' @examples
#' df <- data.frame(chrom=c(rep("chr5", 7), rep("chr21", 3)), start=1:10, end=seq(20, by=10, length.out=10), strand=rep(c("+","+", "-", "*"), length.out=10), letter=letters[1:10], score=rnorm(10))
#' df
#' df2granges(df, assembly="GRCh38_chr")
df2granges <- function(df, ids=rownames(df), chrom.col=1L, start.col=2L, end.col=3L, strand.col=NULL, coord.format="B1RI", assembly=NULL, doSort=FALSE) {
	if (is.character(chrom.col)) { chrom.col <- which(colnames(df) == chrom.col) }
	if (is.character(start.col)) { start.col <- which(colnames(df) == start.col) }
	if (is.character(end.col)) { end.col <- which(colnames(df) == end.col) }
	if (is.character(strand.col)) { strand.col <- which(colnames(df) == strand.col) }

	# get sequence lengths
	validChroms <- c()
	seqlengths <- c()
	if (!is.null(assembly)){
		seqlengths <- getSeqlengths4assembly(assembly)
		validChroms <- names(seqlengths)
	}

	# Process the chromosome names
	chroms <- as.character(df[, chrom.col])
	# chroms <- paste0("chr", sub("^chr", "", chroms))
	param.list <- list()
	param.list[["seqnames"]] <- chroms

	df <- df[!(is.na(df[, start.col]) | is.na(df[, end.col])), ]
	1
	# adjust format
	if (coord.format=="B1RI"){
		#1-based, right-inclusive --> nothing to do
		tmp <- NA
	} else if (coord.format=="B0RE"){
		#0-based, right-exclusive --> adjust start coordinate
		df[, start.col] <- df[, start.col] + 1L
	} else {
		stop("unknown coordinate format")
	}

	# Region coordinates
	param.list[["ranges"]] <- IRanges(start=df[, start.col], end=df[, end.col], names=ids)

	# Match strands
	if (!is.null(strand.col)) {
		param.list[["strand"]] <- matchStrand(df[, strand.col])
	}

	# other columns
	for (cname in colnames(df)[-c(chrom.col, start.col, end.col, strand.col)]) {
		param.list[[cname]] <- df[[cname]]
	}

	# valid chromosomes
	if (length(validChroms) > 0) {
		i.valid <- (chroms %in% validChroms)
		if (sum(i.valid) < 1) stop("No valid chromosome/contig name was matched. Consider using the *_chr version of the assembly fo deal with missing 'chr' prefix")
		param.list <- lapply(param.list, function(x) { x[i.valid] })
	}
	res <- do.call(GRanges, param.list)

	# Assembly info
	if (!is.null(assembly)) {
		seqlevels(res) <- validChroms
		seqlengths(res) <- seqlengths
		genome(res) <- assembly
	}

	# optional sorting
	if (doSort){
		res <- order(as.integer(seqnames(res)), start(res), end(res), as.integer(strand(res)))
		res <- res[oo]
	}
	return(res)
}

