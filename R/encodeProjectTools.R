#' parse.encode.cv.file
#' 
#' parser for the ENCODE projects controlled vocabulary file
#' @param cvFile the file location for the ENCODE controlled vocabulary file. Defaults to the one provided by ENCODE
#' @return a list containing blank line seperated blocks in each element. Each element
#'         is a named list containing the content to each keyword. a keyword is the first word in a line.
#' @export
#' @examples 
#' cv.blocks <- parse.cv.file()
#}
parse.encode.cv.file <- function(cvFile="http://hgdownload.cse.ucsc.edu/goldenPath/encodeDCC/cv.ra"){
	#read the file
	cv.lines <- scan(file=cvFile, what="character",blank.lines.skip=FALSE, comment.char="", sep="\n")
	#separate into text blocks separated by one or multiple blank lines
	blank.i <- which(cv.lines=="")
	blank.i <- c(0,blank.i,length(cv.lines)+1) #add dummy blank lines before and after the file
	#remove comments (important to do this after the blank lines have been determined)
	cv.lines <- sub("#.*","",cv.lines)
	block.range.starts <- blank.i[1:(length(blank.i)-1)] + 1
	block.range.ends   <- blank.i[2:(length(blank.i))]   - 1
	block.ranges <- data.frame(start=block.range.starts,end=block.range.ends)
	#remove multiple consecutive blank lines from consideration as blocks
	block.ranges <- block.ranges[(block.ranges$end - block.ranges$start) > 0,]
	block.lines <- lapply(1:nrow(block.ranges),FUN=function(i){
		cv.lines[block.ranges$start[i]:block.ranges$end[i]]
	})
	#remove comment lines from blocks
	block.lines <- lapply(block.lines,FUN=function(x){x[x!=""]})
	#remove blocks consiting only of comments
	is.comment.block <- sapply(block.lines,FUN=function(x){length(x)<1})
	block.lines <- block.lines[!is.comment.block]
	#parse each block, considering the first word of each line as keyword
	blocks <- lapply(block.lines,FUN=function(x){
		keys <- sub("^(\\S+)\\s+.*$","\\1",x,perl=TRUE)
		vals <- sub("^\\S+\\s+(.*)$","\\1",x,perl=TRUE)
		res <- lapply(unique(keys),FUN=function(k){
			vals[keys==k]
		})
		names(res) <- unique(keys)
		return(res)
	})
	return(blocks)
}

#' get.encode.cell.table
#' 
#' given the ENCODE controlled vocabulary file, retrieves a table characterizing the ENCODE cells
#' @param cvFile the file location for the ENCODE controlled vocabulary file. Defaults to the one provided by ENCODE
#' @return a table containing ENCODE cell annotations
#' @export
#' @examples 
#' ect <- get.encode.cell.table()
#}
get.encode.cell.table <- function(cvFile="http://hgdownload.cse.ucsc.edu/goldenPath/encodeDCC/cv.ra") {
	sel.keys <- c("term","tag","description","karyotype","lineage","orderUrl","organism","protocol","sex","termId","termUrl","tier","tissue","vendorId","vendorName","color")
	blocks <- parse.encode.cv.file(cvFile)
	is.cell.block <- unlist(sapply(blocks,FUN=function(x){x[["type"]]=="Cell Line"}))
	is.cell.block[is.na(is.cell.block)] <- FALSE
	blocks <- blocks[is.cell.block]
	block.tabs <- lapply(blocks,FUN=function(x){
		res <- lapply(sel.keys,FUN=function(k){
			paste(x[[k]],collapse=";")
		})
		names(res) <- sel.keys
		return(as.data.frame(res))
	})
	res <- do.call("rbind",block.tabs)
	rownames(res) <- res$term
	return(res)
}
