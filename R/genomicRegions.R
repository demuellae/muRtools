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
#' @param addAnnotCols add the columns stored in elementMetadata of GRanges
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

#' granges2igv
#' 
#' Save a GRanges object to a IGV file
#' 
#' @param gr GRanges object
#' @param fn filename to save IGV file to
#' @param sc score vector or column in elementMetadata of GRanges
#' @param addAnnotCols add the columns stored in elementMetadata of GRanges
#' @param doSort sort the regions before writing the output
#' @param toTDF convert to TDF file. Requires that "igvtools" is executable from the current path
#' @return result of writing the table (see \code{write.table})
#' @export 
granges2igv <- function(gr, fn, addStrand=FALSE, addAnnotCols=TRUE, doSort=TRUE, toTDF=FALSE){
	if (doSort){
		oo <- order(as.integer(seqnames(gr)),start(gr), end(gr), as.integer(strand(gr)))
		gr <- gr[oo]
	}
	nns <- rep(".", length(gr))
	if (!is.null(names(gr))) nns <- names(gr)
	tt <- data.frame(
		Chromosome=seqnames(gr),
		Start=format(start(gr)-1, trim=TRUE, scientific=FALSE),
		End=format(end(gr), trim=TRUE, scientific=FALSE),
		Feature=nns,
		stringsAsFactors=FALSE
	)
	if (addStrand){
		tt <- data.frame(tt, Strand=strand(gr), stringsAsFactors=FALSE)
	}
	if (addAnnotCols){
		tt <- data.frame(tt, elementMetadata(gr), stringsAsFactors=FALSE)
	}
	res <- write.table(tt, file=fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	if (toTDF){
		assembly <- unique(genome(gr))
		if (length(assembly) != 1){
			stop("Could not convert to TDF: invalid value for genome(GRangesObject)")
		}
		convertRes <- system2("igvtools", c("toTDF", fn, paste0(fn, ".tdf"), assembly), stdout=TRUE)
	}
	invisible(res)
}

#' granges2bed.igv
#' 
#' Save a GRanges object to a bed file which can be displayed by IGV
#' 
#' @param gr        GRanges object
#' @param fn        filename to save bed file to
#' @param trackName track name to be displayed
#' @param scoreCol  the score column (in the GRanges elementMetadata) that is optionally used for coloring
#' @param na.rm     flag indicating whether items with NA score should be removed
#' @param nameCol   the name column (in the GRanges elementMetadata) that is used for labelling the items
#' @param col.cat   color panel for coloring categorical scores
#' @param col.cont  color panel for coloring numerical scores
#' @param col.na    color used for NA scores
#' @param col.range vector of length 2 indicating the range of scores for the color scales to be applied (continuous scores only)
#' @param na.rm     flag indicating whether items with NA score should be removed
#' @param doSort sort the regions before writing the output
#' @return invisibly, the resulting data frame containing the bed file columns
#' @export 
granges2bed.igv <- function(gr, fn, trackName=NULL, scoreCol=NULL, na.rm=FALSE, nameCol=NULL, col.cat=colpal.bde, col.cont=c("#EDF8B1","#41B6C4","#081D58"), col.na="#bdbdbd", col.range=NULL, doSort=TRUE){
	if (doSort){
		oo <- order(as.integer(seqnames(gr)),start(gr), end(gr), as.integer(strand(gr)))
		gr <- gr[oo]
	}
	naCol <- rep(".", length(gr))
	itemNames <- naCol
	if (!is.null(nameCol)){
		itemNames <- elementMetadata(gr)[, nameCol]
	}

	scores <- naCol
	itemColors <- naCol
	col.na.rgb <- as.integer(col2rgb(col.na)[,1])
	sc <- NULL
	if (!is.null(scoreCol)){
		sc <- elementMetadata(gr)[, scoreCol]
		isCategorical <- FALSE

		#score column is categorical
		if (is.character(sc)){
			itemNames.sc <- sc
			isCategorical <- TRUE
		}
		if (is.factor(sc)){
			itemNames.sc <- as.character(sc)
			isCategorical <- TRUE
		}

		#match the colors to names
		if (isCategorical){
			scores <- rep(0, length(gr)) #strangely the score needs to be 0 when itemRgb should be applied in IGV

			lvls <- sort(unique(itemNames.sc))
			if (!is.null(names(col.cat))){
				if (!all(lvls %in% names(col.cat))){
					stop("Not all categories have a mapped color")
				}
			} else {
				col.cat <- rep(col.cat, length.out=length(lvls))
				names(col.cat) <- lvls
			}
			col.cat.rgb.str <- apply(col2rgb(col.cat), 2, FUN=function(x){paste(as.integer(x), collapse=",")})
			names(col.cat.rgb.str) <- names(col.cat)
			itemColors <- col.cat.rgb.str[itemNames.sc]

			if (!is.null(nameCol)){
				itemNames <- paste0(itemNames, " (", itemNames.sc, ")")
			} else {
				itemNames <- itemNames.sc
			}
		}

		#score column is numeric
		if (is.numeric(sc)){
			#rescale to be between 0 and 1000
			# scores <- as.integer(round(rescale(sc, c(0, 1000))))
			scores <- rep(0, length(gr)) #strangely the score needs to be 0 when itemRgb should be applied in IGV

			#match the values to colors
			cr <- colorRamp(col.cont)
			itemColors <- cr(rescale(sc, c(0, 1)))
			if (!is.null(col.range)){
				#extend the scores with the min and max values to rescale the color scheme
				sc.padded <- c(col.range[1], sc, col.range[2])
				#truncate if values exceed the range
				sc.padded[sc.padded < col.range[1]] <- col.range[1]
				sc.padded[sc.padded > col.range[2]] <- col.range[2]
				itemColors <- cr(rescale(sc.padded, c(0, 1)))[2:(length(sc)+1),]
			}
			itemColors[is.na(sc),] <- col.na.rgb
			itemColors <- apply(itemColors, 1, FUN=function(x){paste(as.integer(x), collapse=",")})

			#write the core to the name column
			scStr <- signif(sc)
			if (!is.null(nameCol)){
				itemNames <- paste0(itemNames, " (", scStr, ")")
			} else {
				itemNames <- scStr
			}

		}
	}
	strands <- as.character(strand(gr))
	strands[strands=="*"] <- "."


	starts <- format(start(gr)-1, trim=TRUE, scientific=FALSE)
	ends   <- format(end(gr), trim=TRUE, scientific=FALSE)
	tt <- data.frame(
		chrom=seqnames(gr),
		start=starts,
		end=ends,
		itemName=itemNames,
		score=scores,
		strand=strands,
		thickStart=starts,
		thickEnds=ends,
		itemRgb=itemColors,
		stringsAsFactors=FALSE
	)

	#discard rows where the score is NA if demanded
	if (na.rm & !is.null(sc)){
		tt <- tt[!is.na(sc), ]
	}

	if (is.null(trackName)) trackName <- "GRanges track"

	trackLine <- paste0('track name="', trackName, '" description="', trackName, '" visibility=1 itemRgb="On"')
	write(trackLine, file=fn)
	write.table(tt, file=fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, na=".", append=TRUE)
	invisible(tt)
}
#granges2bed.igv(gr, "~/tmp/gr_conv.bed", trackName="blubb", scoreCol="cmp_LiHe_CTvST")
#granges2bed.igv(gr, "~/tmp/gr_conv_cat.bed", trackName="blubb", scoreCol="category")

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
#' @param adjNumChromNames Should numeric chromosome names be adjusted for by adding the prefix "chr"
#' @return \code{GRanges} object encapsulating of regions included in \code{df}.
#' 	       As GRanges, the coordinates will be 1-based right-inclusive.
#'         Columns other that the ones listed as parameters in this function are included as elementMetadata.
#'
#' @export
#' @examples
#' df <- data.frame(chrom=c(rep("chr5", 7), rep("chr21", 3)), start=1:10, end=seq(20, by=10, length.out=10), strand=rep(c("+","+", "-", "*"), length.out=10), letter=letters[1:10], score=rnorm(10))
#' df
#' df2granges(df, assembly="GRCh38_chr")
df2granges <- function(df, ids=rownames(df), chrom.col=1L, start.col=2L, end.col=3L, strand.col=NULL, coord.format="B1RI", assembly=NULL, doSort=FALSE, adjNumChromNames=FALSE) {
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
	if (adjNumChromNames){
		chroms <- sub("^([1-9XYMm]+$)", "chr\\1", chroms)
	}
	param.list <- list()
	param.list[["seqnames"]] <- chroms

	df <- df[!(is.na(df[, start.col]) | is.na(df[, end.col])), ]

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
		oo <- order(as.integer(seqnames(res)), start(res), end(res), as.integer(strand(res)))
		res <- res[oo]
	}
	return(res)
}
