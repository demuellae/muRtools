#' getTxDb.gencode
#' 
#' Create a \code{TxDb} object by downloading the corresponding GTF file from Gencode
#' @param name	gencode identifier. Currently supported are: "gencode_v27", "gencode_v19", "gencode_vM16", "gencode_vM1"
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
