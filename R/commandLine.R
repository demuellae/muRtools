#' parse.cl.args
#' 
#' parser for command line arguments
#' @return a named list with command line arguments
#' @export
#' @aliases ggplot2.heatmap
#' @examples 
#' cmd.args <- parse.cl.args()
#}
parse.cl.args <- function(){
	args <- commandArgs(trailingOnly = TRUE)
	args.split <- strsplit(args,"=")
	res <- list()
	if (length(args.split) < 1){
		return(res)
	}
	for (i in 1:length(args.split)){
		x <- args.split[[i]]
		if(length(x)>2){
			stop("Invalid command line arguments. Too many '='.")
		} else {
			if (length(x)==2){
				if (is.element(x[1],names(res))) {
					res[[x[1]]] <- c(res[[x[1]]],x[2])
				} else {
					res[[x[1]]] <- x[2]
				}
			} else {
				res <- c(res,list(NULL))
				names(res)[length(res)] <- x[1]
			}
		}
	}
	return(res)
}
