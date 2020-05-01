
#' chordDiagramFromContingencyTable
#' 
#' Plot a bipartite chord diagram (circular Sankey diagram with 2 categories) for a contingency matrix
#' the values of two matrices
#' @param contTab   \code{matrix} or \code{table} containing the contingency matrix
#' @param chordColorByCol color chords by column instead of by row
#' @param cs_rows   color scheme to use for the rows of the matrix
#' @param cs_columns   color scheme to use for the rows of the matrix
#' @return nothing of particular interest (include this function while plotting).
#' @author Fabian Mueller
#' @export
#' @examples
#' # contingency table of air quality quantile by month
#' contTab <- with(airquality, table(cut(Temp, quantile(Temp)), Month))
#' names(dimnames(contTab))[1] <- "quantile"
#' chordDiagramFromContingencyTable(contTab)
chordDiagramFromContingencyTable <- function(contTab, chordColorByCol=FALSE, cs_rows=colpal.mu.cat, cs_columns=colpal.mu.cat){
	require(circlize)

	if (is.character(cs_rows) && length(cs_rows)==1 && cs_rows=="[auto]"){
		cs_rows <- colpal.mu.cat
	}
	if (is.character(cs_columns) && length(cs_columns)==1 && cs_columns=="[auto]"){
		cs_columns <- colpal.mu.cat
	}

	# contTab <- t(contTab)
	dn1 <- names(dimnames(contTab))[1]
	dn2 <- names(dimnames(contTab))[2]

	# convert table to matrix (circlize cannot deal with table objects)
	ovMat <- matrix(contTab, nrow=nrow(contTab))
	rownames(ovMat) <- paste0(dn1, "_", rownames(contTab))
	colnames(ovMat) <- paste0(dn2, "_", colnames(contTab))


	# orient the segments to be bipartite
	# following the documentation in http://zuguang.de/circlize_book/book/advanced-usage-of-chorddiagram.html
	row_sum <- sum(rowSums(abs(ovMat)))
	col_sum <-sum(colSums(abs(ovMat)))
	small_gap <- 1
	big_gap   <- 20
	nr <- nrow(ovMat)
	nc <- ncol(ovMat)
	n_sector <- nr + nc
	row_sector_degree <- (360 - small_gap*(n_sector - 2) - big_gap*2) * (row_sum/(row_sum + col_sum)) + small_gap*(nr-1)
	# start_degree <- 90 - (180 - row_sector_degree)/2
	start_degree <- 270 - (180 - row_sector_degree)/2
	gaps <- c(rep(small_gap, nr - 1), big_gap, rep(small_gap, nc - 1), big_gap)

	# segment colors
	stateColors <- c(
		rep(cs_rows, length.out=nrow(contTab)),
		rep(cs_columns, length.out=ncol(contTab))
	)
	names(stateColors) <- c(rownames(ovMat), colnames(ovMat))

	circos.par(gap.after=gaps, start.degree=start_degree)
	# circos.par(start.degree = 90, clock.wise = FALSE)
	# circos.par(gap.after=c(rep(2, nrow(ovMat)-1), 20, rep(2, ncol(ovMat)-1), 20))
	column.col <- NULL
	if (chordColorByCol) {
		column.col <- rep(cs_columns, length.out=ncol(contTab))
	}
	chordDiagram(ovMat, order=c(rev(rownames(ovMat)), colnames(ovMat)), grid.col=stateColors, column.col=column.col, annotationTrack=c("name", "grid"))
	circos.clear()
	invisible(NULL)
}
