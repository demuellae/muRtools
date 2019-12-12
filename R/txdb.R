#' getTxDb.gencode
#' 
#' Create a \code{TxDb} object by downloading the corresponding GTF file from Gencode
#' @param name	gencode identifier. Currently supported are: "gencode.v27", "gencode.v19", "gencode.vM16", "gencode.vM1"
#' @return \code{TxDb} object
#' @author Fabian Mueller
#' @export 
getTxDb.gencode <- function(name){
	require(GenomicFeatures)
	fileTab <- data.frame(
		name=c("gencode.v27", "gencode.v19", "gencode.vM16", "gencode.vM1"),
		version=c("27", "19", "M16", "M1"),
		fileURL=c(
			"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz",
			"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
			"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz",
			"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz"
		),
		organism = c("Homo sapiens", "Homo sapiens", "Mus musculus", "Mus musculus"),

		genomeAss=c("hg38", "hg19", "mm10", "mm9"),
		stringsAsFactors = FALSE
	)
	rownames(fileTab) <- fileTab$name
	if (!is.element(name, fileTab$name)){
		logger.error(c("Gencode name not supported:", name, "(must be one of", paste(fileTab$name, collapse=", "), ")"))
	}
	chrInfo <- getSeqlengths4assembly(fileTab[name, "genomeAss"])
	chrInfo <- Seqinfo(names(chrInfo), chrInfo, genome=fileTab[name, "genomeAss"])
	txdb <- makeTxDbFromGFF(fileTab[name, "fileURL"], format="gtf", dataSource=paste("Gencode version", fileTab[name, "version"]), circ_seqs=character(), organism=fileTab[name, "organism"], chrominfo=chrInfo)
	return(txdb)
}

#' getAnnotGrl.gencode
#' 
#' Create a \code{GRangesList} with element annotation by downloading the corresponding GTF file from Gencode
#' @param name	gencode identifier. Currently supported are: "gencode.v27", "gencode.v19", "gencode.vM16", "gencode.vM1"
#' @return \code{GRangesList} object with annotated elements for each element type (genes, transcripts, exons, ...)
#' @author Fabian Mueller
#' @export 
getAnnotGrl.gencode <- function(name){
	require(GenomicFeatures)
	require(rtracklayer)
	fileTab <- data.frame(
		name=c("gencode.v27", "gencode.v19", "gencode.vM16", "gencode.vM1"),
		version=c("27", "19", "M16", "M1"),
		fileURL=c(
			"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz",
			"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
			"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz",
			"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz"
		),
		organism = c("Homo sapiens", "Homo sapiens", "Mus musculus", "Mus musculus"),

		genomeAss=c("hg38", "hg19", "mm10", "mm9"),
		stringsAsFactors = FALSE
	)
	rownames(fileTab) <- fileTab$name
	if (!is.element(name, fileTab$name)){
		logger.error(c("Gencode name not supported:", name, "(must be one of", paste(fileTab$name, collapse=", "), ")"))
	}
	gr <- import.gff(fileTab[name, "fileURL"])
	gr <- setGenomeProps(gr, fileTab[name, "genomeAss"])
	gr <- sortGr(gr)
	grl <- split(gr, elementMetadata(gr)[,"type"])
	return(grl)
}


#' getGeneAnnotMap
#' 
#' Get a mapping (e.g. of identifiers) and automatically select the correct \code{AnnotationDbi} database for a given
#' assembly
#' @param assembly	character string specifying the assembly
#' @param from      the column name that will be used as the key for the resulting map
#' @param to        the column name that will be used as the result for the resulting map
#' @param multiMap  character string specifying what to do if multiple mappings are found for a key
#'                  by default (\code{multiMap="paste"}) the results will be pasted into a single character string (separated by ';').
#'                  Other options include \code{'first'} for just returning the first value or 'list' for returning a list of all values
#' @return a named vector (or list depending on how the \code{multiMap} argument is chosen) providing a mapping
#' @author Fabian Mueller
#' @export
getGeneAnnotMap <- function(assembly, from="ENSEMBL", to="SYMBOL", multiMap="paste"){
	aa <- NULL
	if (is.element(assembly, c("hg38", "hg19"))){
		require(org.Hs.eg.db)
		aa <- org.Hs.eg.db::org.Hs.eg.db
	} else if (is.element(assembly, c("mm9", "mm10"))){
		require(org.Mm.eg.db)
		aa <- org.Mm.eg.db::org.Mm.eg.db
	} else {
		stop(paste0("Unknown assembly:", assembly))
	}
	kk <- keys(aa, keytype=from)
	if (multiMap=="paste"){
		res <- mapIds(aa, keys=kk, column=to, keytype=from, multiVals="list")
		res <- sapply(res, FUN=function(x){paste(x, collapse=";")})
	} else {
		res <- mapIds(aa, keys=kk, column=to, keytype=from, multiVals=multiMap)
	}
	return(res)
}

