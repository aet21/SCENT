#' @title 
#' Correlation of Connectome And Transcriptome
#' 
#' @aliases CCAT
#'  
#' @description 
#' This function leverages Pearson correlation between gene expresion level
#' and gene connectome derived from PPI network to fastly estimate signaling
#' entropy rate.
#' 
#' 
#' @param exp.m
#' A scRNA-Seq data matrix with rows labeling genes and columns 
#' labeling single cells and with rownames annotated to a gene-identifier
#' identical to that used in the \code{ppiA.m} argument. scRNA-Seq data matrix
#' should have undergone prior QC to remove poor quality cells and each cell
#' normalized by library size. If data has not been log-transformed, the
#' function will log2-transform with a pseudocount of +1.
#' 
#' @param ppiA.m
#' The adjacency matrix of a user-given PPI network with rownames and 
#' colnames labeling genes (same gene identifier as in \code{exp.m}). Diagonal
#' entries should be zero and the number of genes in the network should be large,
#' i.e. at least over 8000, to ensure a reasonably large overlap of genes with
#' those in the expression data matrix.
#' 
#' 
#' 
#' @return CCAT
#' The estimated CCAT values as a vector
#'  
#' @references 
#' Teschendorff Andrew E., Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' 
#' 
#' @examples 
#' 
#' 
#' @import Matrix
#' @importFrom qlcMatrix corSparse
#'
#' @export
#'     

CompCCAT <- function(exp.m, ppiA.m){

    if(max(exp.m) > 100){ ### check if data has been log2-transformed or not and if not, then log2-transform with a pseudcount of 1, as for CCAT we can allow 0s.
        exp.m <- log2(exp.m+1);
    }
    # get input class of data matrix
    classMATRIX <- class(exp.m);
    
    # set common gene IDs
    commonEID.v <- intersect(rownames(ppiA.m),rownames(exp.m));
    
    # check the consistency of gene IDs and size of overlap
    if ( length(commonEID.v) < 5000 ){
      stop(paste("The overlap of common genes between PPI and expression matrix is only ",length(commonEID.v)," so check that gene identifiers are correct. We don't recommend running CCAT on less than 5000 overlapping genes.",sep=""));
    }
    
    # compute degrees
    k.v <- rowSums(ppiA.m[match(commonEID.v, rownames(ppiA.m)),]);

    if(classMATRIX=="matrix"){ ## ordinary matrix
      ccat.v <- as.vector(cor(exp.m[match(commonEID.v,rownames(exp.m)),],k.v));
    }
    else if (classMATRIX=="dgCMatrix"){
      ccat.v <- as.vector(corSparse(exp.m[match(commonEID.v,rownames(exp.m)),],Matrix(matrix(k.v,ncol=1))));
    }
    
    # subsample data set
  #  if (parallelMode) {
#      chunk <- round(length(col_names)/subsamplesize)
#      subsamples <- split(1:length(col_names), sample(factor(1:length(col_names) %% chunk)))
#    }else if (length(col_names) > 10000) {
#      chunk <- round(length(col_names)/1000)
#      subsamples <- split(1:length(col_names), sample(factor(1:length(col_names) %% chunk)))
#    }else{
  #    subsamples <- NULL
#    }
    
 return(ccat.v);   
    
} ## EOF




