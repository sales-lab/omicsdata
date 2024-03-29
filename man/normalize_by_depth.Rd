% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression.R
\name{normalize_by_depth}
\alias{normalize_by_depth}
\title{Normalize expression profiles by sequencing depth.}
\usage{
normalize_by_depth(counts)
}
\arguments{
\item{counts}{The raw count matrix.}
}
\value{
Trasformed counts.
}
\description{
Normalize expression profiles by sequencing depth.
}
\examples{
library(airway)
data(airway)
normalize_by_depth(assays(airway)$counts)

}
