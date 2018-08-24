################################################################################
# Utilities for enrichment analysis
################################################################################
#' goEnrichment
#' 
#' Perform Gene Ontology (GO) enrichment using the \code{topGO} package
#' @param qids      character vector of query gene IDs
#' @param uids      character vector of universe gene IDs
#' @param ontology  character specifying the ontology to use (default: \code{"BP"})
#' @param idType    character specifying which universe the gene IDs come from (default: \code{"ensembl"}). Possible values are: \code{"entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene"}
#' @param assembly  character specifying the genome to use. Can either be the name of the package to be used for mapping the identifiers (e.g. \code{"org.Hs.eg.db"}; default)
#'                  or an identifier for a genome assembly ((e.g. \code{"hg38"})
#' @param algorithm algorithm employed by \code{topGO}. See \code{topgGO::runTest} for details.
#' @return a list (S3 class) object containing:
#'         - \code{$tgData}:    the used \code{topGOdata} object
#'         - \code{$resultObj}: the resulting \code{topGOresult} object
#'         - \code{$table}:     a summary table of statistics for each GO term
#' @author Fabian Mueller
#' @export 
goEnrichment <- function(qids, uids, ontology="BP", idType="ensembl", assembly="org.Hs.eg.db", algorithm="weight01"){
	require(topGO)
	if (!all(qids %in% uids)) stop("query ids must be a subset of universe ids")
	ass <- assembly
	if (is.element(assembly,c("hg19","hg38"))){
		ass <- "org.Hs.eg.db"
	} else if (is.element(assembly,c("mm9","mm10"))){
		ass <- "org.Mm.eg.db"
	} else {
		if(!is.element(assembly, c("org.Hs.eg.db", "org.Mm.eg.db"))) stop("Unsupported assembly")
	}
	geneList <- rep(0, length(uids))
	names(geneList) <- uids
	geneList[qids] <- 1

	tgData <- new("topGOdata",
		ontology=ontology,	allGenes=geneList, geneSel=function(x){x>0},
		annot=annFUN.org, mapping=ass, ID=idType
	)
	enrichRes <- runTest(tgData, algorithm=algorithm, statistic="fisher")

	pVals <- score(enrichRes)
	enrichTab <- GenTable(
		tgData, pVal=enrichRes,
		topNodes=length(pVals),
		numChar=1e6
	)
	rownames(enrichTab) <- enrichTab[,"GO.ID"]
	enrichTab[,"pVal"] <- pVals[rownames(enrichTab)]

	res <- list(
		tgData=tgData,
		resultObj=enrichRes,
		table=enrichTab
	)
	class(res) <- "GOenrichment [muRtools:topGO]"
	return(res)
}
