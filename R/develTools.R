#' reloadPackage
#' 
#' reloads a package without quitting the R session. Useful for developing packages interactively
#' @param package.name name (string) of the package to be reloaded
#' @return result of \code{library(package.name)}
#' @export 
#' @examples 
#' library(GenomicRanges)
#' reloadPackage("GenomicRanges")
reloadPackage <- function(package.name){
	detach(paste0("package:",package.name), unload=TRUE, character.only=TRUE)
	library(package.name,character.only=TRUE)
}
