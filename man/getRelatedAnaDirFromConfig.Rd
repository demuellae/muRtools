% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysisConfig.R
\name{getRelatedAnaDirFromConfig}
\alias{getRelatedAnaDirFromConfig}
\title{getRelatedAnaDirFromConfig}
\usage{
getRelatedAnaDirFromConfig(cfg, anaName, anaVersion = cfg[[".anaVersion"]])
}
\arguments{
\item{cfg}{configuration object as returned by \code{getConfig}}

\item{anaName}{name of the analysis}

\item{anaVersion}{version of the analysis. If \code{NULL}, it will look for the most recent one}
}
\value{
path of the related analysis directory
}
\description{
Get analysis configuration object from a JSON file. Automatically runs parsers
for config elements specified in the '.MU_ANA_CONFIG' section of the JSON file
}
\author{
Fabian Mueller
}
