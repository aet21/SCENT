\docType{data}
\name{dataChu}
\alias{dataChu}
\title{scRNA-Seq C1 Fluidigm data from Chu et al}
\format{Contains a subset of the 1018 single cells profiled in Chu et al,
  consisting of 374 hESCs and over 100 cells representing multipotent
  progenitors of endothelial cells, together with the cell-type
  information, and the SR-values.
  
}

\usage{
data(dataChu)
}
\description{
The data has been normalized for library-size with cell-dependent
scale-factors, but has not been log-transformed.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{scChu.m}: the scRNA-Seq data matrix
  \item \code{phenoChu.v}: a vector of cell-type labels matching the
  columns of \code{scChu.m}.
  \item \code{srChu.v}: the SR-values for the cells matching the columns
  of \code{scChu.m}.
}
}
\references{
Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of
differentiation potency from a cell’s transcriptome. Nat.Commun.2017
}
\keyword{datasets}
