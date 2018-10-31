
#' Custom Color Paletes
#' 
#' \describe{
#'   \item{\code{colpal.cb}}{
#'         color blind friendly color palettes (adapted from http://wiki.stdout.org/rcookbook/Graphs/Colors%20%28ggplot2%29/)
#'   }
#' }
#'	
#' @export
#' @rdname colpal
#' @aliases colpal,colpal.cb,colpal.colpal.bde
#' @examples 
#' plot.colpal(colpal.cb)
colpal.cb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000")
#' \describe{
#'   \item{\code{colpal.bde}}{
#'         Enhanced Color Brewer palette 'Dark2'
#'   }
#' }
#' @rdname colpal
#' @export
#' @examples
#' plot.colpal(colpal.bde)
colpal.bde <- c("#2166AC","#B2182B","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#00441B","#40004B","#053061")
#' \describe{
#'   \item{\code{colpal.nature}}{
#'         Color palette inspired by Nature journal color scheme's.
#'   }
#' }
#' @rdname colpal
#' @export
#' @examples
#' plot.colpal(colpal.nature)
colpal.nature <- c("#003D7C", "#D50911", "#0086A8", "#008136", "#7C68A4", "#8E1A47", "#E67800", "#709F28", "#008FB4", "#84486A", "#B5797F", "#7489A8",  "#6C9396", "#7D9FB1", "#84486A", "#7C698B", "#88A2C3")
#' \describe{
#'   \item{\code{colpal.nature}}{
#'         Color palette inspired by Nature journal color scheme's.
#'   }
#' }
#' @rdname colpal
#' @export
#' @examples
#' plot.colpal(colpal.mu.cat)
colpal.mu.cat <- c("#e69f00", "#56b4e9", "#74c476", "#cc79a7", "#d55e00", "#0072b2", "#009e73", "#6a51a3", "#f0e442", "#999999", "#000000")
#' \describe{
#'   \item{\code{colpal.corpid}}{
#'         Color palette inspired by coorporate identities I worked with
#'   }
#' }
#' @rdname colpal
#' @export
#' @examples
#' plot.colpal(colpal.corpid)
colpal.corpid <- c(
  "mpii.darkblue"="#1C1D3B",
  "mpii.lightblue"="#676D8D",
  "mpii.grey"="#C6C6B6",
  "deep.turq"="#018D9D",
  "deep.grey"="#C6C6C6",
  "stanford.cardinal"="#8c1515",
  "stanford.coolgrey"="#4d4f53",
  "stanford.black"="#2e2d29",
  "stanford.brightred"="#B1040E",
  "stanford.chocolate"="#2F2424",
  "stanford.fog"="#F4F4F4",
  "stanford.cloud"="#dad7cb",
  "stanford.lightsandstone"="#F9F6EF"
)
#' \describe{
#'   \item{\code{colpals.games}}{
#'         List of color palette inspired by board games
#'   }
#' }
#' @rdname colpal
#' @export
#' @examples
#' plot.colpal(colpals.games[["rollgalaxy"]])
colpals.games <- list(
  bruges=c("#3B8FCF", "#B52622", "#8D5718", "#572978", "#FFC856", "#28A742", "#256CBF", "#C4292C", "#E0EB36"),
  mombasa=c("#009FE3", "#8EC041", "#951B81", "#FFCC00", "#BE1716", "#DE7E00", "#7C83B3", "#671719", "#E0CDA6", "#775725", "#000000"),
  dominantspecies=c("#0093E9", "#D93123", "#19B245", "#F9ED00", "#E5BCCF", "#F7D58A", "#9FC6E0", "#ACD499", "#D8D088", "#59413B", "#000000"),
  terramystica=c("#3B76BB", "#A9122A", "#235E31", "#996F58", "#EDE980", "#D69A30", "#000000", "#B6B6B6"),
  gaiaproject=c("#005AA6", "#D52429", "#FFCB28", "#F37124", "#744F30", "#858A91", "#56A744", "#A52181", "#697F9B", "#0C9DD9"),
  rollgalaxy=c("#6DCFF6", "#ED1C24", "#8DC63F", "#BD8CBF", "#FFF200", "#DBAC78", "#00A651", "#2E3192", "#EC008C")
  #c("#", "#", "#", "#", "#", "#", "#", "#")
)

#' plot.colpal
#' 
#' Get a continuous color palette
#' @param cp   color palette, i.e. vector of colors
#' @param type pie chart or stripes
#' @return nothing of particular interest
#' @author Fabian Mueller
#' @export 
plot.colpal <- function(cp, type="pie"){
  if (type=="pie"){
    pie(rep(1,length(cp)), labels=names(cp), col=cp, clockwise=TRUE)
  } else if (is.element(type, c("bar", "stripes"))){
    require(plotrix)
    plot.new()
    gradient.rect(0,0,1,1,col=cp,nslices=length(cp),gradient="x",border=NA)
  } else {
    stop(c("Unknown type for color plotting:", type))
  }
}

#' colpal.cont
#' 
#' Get a continuous color palette
#' @param n   number of colors returned
#' @param name   name of the color palette
#' @param ...  arguments passed to other functions
#' @return a character vector containing n colors
#' @author Fabian Mueller
#' @export
#' @examples
#' plot.colpal(colpal.cont(5, "viridis"))
#' plot.colpal(colpal.cont(5, "cb.BrBG"))
colpal.cont <- function(n=3, name="viridis", ...){
  if (name=="viridis"){
    require(viridis)
    return(viridis(n, ...))
  } else if (grepl("^cb\\.", name)){
    require(RColorBrewer)
    name <- gsub("^cb\\.", "", name)
    return(brewer.pal(n, name, ...))
  } else {
    stop(c("Unknown name for color palette:", name))
  }
}
