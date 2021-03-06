% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomicRegions.R
\name{grGeneAnnot}
\alias{grGeneAnnot}
\title{grGeneAnnot}
\usage{
grGeneAnnot(
  gr,
  rsdb,
  geneSetName = "genes_protein_coding",
  geneSetCollection = "Gencode",
  maxDist = 1e+05
)
}
\arguments{
\item{gr}{\code{GRanges} object to liftOver}

\item{rsdb}{\code{RegionSetDB} object containing a region set database from which gene annotation can be retrieved}

\item{geneSetName}{Name of the region set containng gene annotation in the \code{RegionSetDB}}

\item{geneSetCollection}{Name of the region set collection containng gene annotation in the \code{RegionSetDB}}

\item{maxDist}{maximum distance for matching to nearest gene}
}
\value{
\code{data.frame} containing information on the nearest gene for each element in \code{gr}
}
\description{
get gene annotation for a \code{GRanges} object using a \code{RegionSetDB} region database object by linking to the nearest gene
}
