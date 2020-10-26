\docType{data}
\name{net17Jan16}
\alias{net17Jan16}
\title{PPI network from Pathway Commons-v2}
\format{Adjacency matrix for the PPI network  
}

\usage{
data(net17Jan16)
}
\description{
Adjacency matrix describing a curated PPI network derived from Pathway
Commons on the 17th Jan 2016, containing over 10000 proteins. Proteins
have been annotated to the genes (Entrez gene ID) they encode. Diagonal
entries are 0. Other entries are zero or 1 if an interaction has been
reported. Please see references for detailed construction of the PPI network.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{net17Jan16.m}: the adjacency matrix
}
}
\references{
Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of
differentiation potency from a cell’s transcriptome. Nat.Commun.2017
}
\keyword{datasets}
