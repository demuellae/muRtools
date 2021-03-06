% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lolaUtils.R
\name{loadLolaDbs}
\alias{loadLolaDbs}
\title{loadLolaDbs}
\usage{
loadLolaDbs(lolaDbPaths)
}
\arguments{
\item{lolaDbPaths}{vector of names of LOLA DB paths to be loaded}
}
\value{
LOLA DB list as returned by \code{LOLA::loadRegionDB}
}
\description{
Load LOLA databases from disk and merge them
}
\examples{
\donttest{
# download LOLA DB
lolaDest <- tempfile()
dir.create(lolaDest)
lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
}
}
\author{
Fabian Mueller
}
