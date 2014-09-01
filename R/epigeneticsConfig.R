#' Custom Color Paletes for epigenetic modifications
#' 
#' Named vectors of colors for different modifications
#' \describe{
#'   \item{\code{colpal.histone}}{
#'   	Histone modifications
#'   }
#' }
#'	
#' @export
#' @rdname colpal.epigenetics
#' @aliases colpal.histone,colpal.histone.ihec,colgrad.methylation.rb,colgrad.methylation.yb
#' @examples 
#' library(gplots)
#' pie(rep(1,length(colpal.histone)), labels=names(colpal.histone), col=colpal.histone)
colpal.histone <- c(
	"H3K4me1"  = "#7FCDBB",
	"H3K4me2"  = "#ADDD8E",
	"H3K4me3"  = "#238B45",
	"H3K9me1"  = "#FB6A4A",
	"H3K9me3"  = "#D94801",
	"H3K27me3" = "#A50F15",
	"H3K27ac"  = "#0570B0",
	"H3K36me3" = "#6A51A3",
	"Input"    = "#525252",
	"other"    = "#252525"
)

#' \describe{
#'   \item{\code{colpal.histone.ihec}}{
#'      Histone modifications. Colors defined by IHEC.
#'   }
#' }
#' @rdname colpal.epigenetics
#' @export
#' @examples
#' library(gplots)
#' pie(rep(1,length(colpal.histone.ihec)), labels=names(colpal.histone.ihec), col=colpal.histone.ihec)
colpal.histone.ihec <- c(
	"H3K4me1"  = "#B2DF8A",
	"H3K4me3"  = "#1F78B4",
	"H3K9me3"  = "#A6CEE3",
	"H3K27me3" = "#33A02C",
	"H3K27ac"  = "#E31A1C",
	"H3K36me3" = "#FB9A99",
	"Input"    = "#FF7F00",
	"other"    = "#FDBF6F"
)

#' \describe{
#'   \item{\code{colgrad.methylation.rb}}{
#'      3 point color palette for methylation: red-grey-blue
#'   }
#' }
#' @rdname colpal.epigenetics
#' @export
#' @examples
#' library(plotrix)
#' cp <- colorpanel(100,colgrad.methylation.rb["low"],colgrad.methylation.rb["mid"],colgrad.methylation.rb["high"])
#' plot.new()
#' gradient.rect(0,0,1,1,col=cp,nslices=length(cp),gradient="x",border=NA)
colgrad.methylation.rb <- c(low="#AD0021",mid="#909090",high="#39278C")
#' \describe{
#'   \item{\code{colgrad.methylation.rb}}{
#'      3 point color palette for methylation: red-grey-blue
#'   }
#' }
#' @rdname colpal.epigenetics
#' @export
#' @examples
#' library(plotrix)
#' cp <- colorpanel(100,colgrad.methylation.yb["low"],colgrad.methylation.yb["mid"],colgrad.methylation.yb["high"])
#' plot.new()
#' gradient.rect(0,0,1,1,col=cp,nslices=length(cp),gradient="x",border=NA)
colgrad.methylation.yb <- c(low="#EDF8B1",mid="#41B6C4",high="#081D58")
