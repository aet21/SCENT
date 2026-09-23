\docType{data}
\name{dataLiver}
\alias{dataLiver}
\title{scRNA-Seq data from Yang et al}
\format{Normalized and log-transformed scRNA-Seq data representing a
  differentiation timecourse of hepatoblasts into hepatocytes and
  cholangiocytes. Genes have been annotated to unique human Entrez gene
  IDs. There are 447 cells in total and 7 different timepoints
  (embryonic stages).
}

\usage{
data(dataLiver)
}
\description{
The data has been normalized for library-size with cell-dependent
scale-factors, and has been log-transformed with a pseudocount of 1.1
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{scLiver.m}: the scRNA-Seq data matrix.
  \item \code{phenoLiver.v}: a vector of timepoints matching the
  columns of \code{scLiver.m}.
  \item \code{colorLiver.v}: a color vector for the timepoints.
}
}
\references{
Yang L et al. A single-cell transcriptomic analysis reveals precise
pathways and regulatory mechanisms underlying hepatoblast
differentiation. Hepatology 2017.
}
\keyword{datasets}
