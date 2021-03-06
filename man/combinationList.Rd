% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listOps.R
\name{combinationList}
\alias{combinationList}
\title{combinationList}
\usage{
combinationList(...)
}
\arguments{
\item{...}{vectors of elements. Ideally named}
}
\value{
a list containing all combinations of elements in the input. Each element contains a unique combination
}
\description{
get a list of all combinations of vectors. Basically a wrapper around \code{\link{expand.grid}}
}
\examples{
combinationList(a=letters[1:5], A=LETTERS[1:3], i=1:4) 
}
