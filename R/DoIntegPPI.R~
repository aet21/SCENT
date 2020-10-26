#' @title 
#' Integration of gene expression matrix and PPI network
#' 
#' @aliases DoIntegPPI
#'  
#' @description 
#' This function finds the common genes between the scRNA-Seq data matrix 
#' and the genes present in the PPI network, and constructs the maximally 
#' connected subnetwork and reduced expression matrix for the computation 
#' of signaling entropy.
#' 
#' @param exp.m
#' Can be three major kinds of input:
#' One is a scRNA-Seq data matrix with rows labeling genes and columns 
#' labeling single cells. And it can be either a log-transformed data
#' matrix with minimal value around 0.1 (recommended), or an 
#' nonlog-transformed data matrix with minimal value 0.
#' The other two kinds of input can be either a "SingleCellExperiment"
#' class object or a "CellDataSet" class object
#' 
#' @param ppiA.m
#' The adjacency matrix of a user-given PPI network with rownames and 
#' colnames labeling genes (same gene identifier as in \code{exp.m})
#' 
#' @param log_trans
#' A logical. Whether to do log-transformation on the input data
#' matrix or not. Default is FALSE
#' 
#' @return A list of two or four objects:
#' 
#' @return expMC
#' Reduced expression matrix with genes in the maximally connected subnetwork
#' 
#' @return adjMC
#' Adjacency matrix of the maximally connected subnetwork
#' 
#' @return data.sce/data.cds
#' Orginal input sce/cds data objects
#' 
#' @return data
#' Normalized data matrix
#' 
#' @return degree.v
#' The nodes(gene) degree in the integrated network
#' 
#' @return dgC.m
#' An additional dgCMatrix object of the input data, for the convenience 
#' of later on calculation
#'  
#' @references 
#' Chen, Weiyan, et al.
#' \emph{Single-cell landscape in mammary epithelium reveals 
#' bipotent-like cells associated with breast cancer risk 
#' and outcome.}
#' Communications Biology 2 (2019): 306.
#' doi:\href{https://doi.org/10.1038/s42003-019-0554-8}{
#' 10.1038/s42003-019-0554-8}. 
#' 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Teschendorff AE, Banerji CR, Severini S, Kuehn R, Sollich P. 
#' \emph{Increased signaling entropy in cancer requires the scale-free 
#' property of protein interaction networks.}
#' Scientific reports 5 (2015): 9646.
#' doi:\href{https://doi.org/10.1038/srep09646}{
#' 10.1038/srep09646}.
#' 
#' Banerji, Christopher RS, et al. 
#' \emph{Intra-tumour signalling entropy determines clinical outcome 
#' in breast and lung cancer.}
#' PLoS computational biology 11.3 (2015): e1004115.
#' doi:\href{https://doi.org/10.1371/journal.pcbi.1004115}{
#' 10.1371/journal.pcbi.1004115}.
#' 
#' Teschendorff, Andrew E., Peter Sollich, and Reimer Kuehn.
#' \emph{Signalling entropy: A novel network-theoretical framework 
#' for systems analysis and interpretation of functional omic data.}
#' Methods 67.3 (2014): 282-293.
#' doi:\href{https://doi.org/10.1016/j.ymeth.2014.03.013}{
#' 10.1016/j.ymeth.2014.03.013}.
#' 
#' Banerji, Christopher RS, et al. 
#' \emph{Cellular network entropy as the energy potential in 
#' Waddington's differentiation landscape.}
#' Scientific reports 3 (2013): 3039.
#' doi:\href{https://doi.org/10.1038/srep03039}{
#' 10.1038/srep03039}.
#' 
#' @examples 
#' ### define a small network
#' ppiA.m <- matrix(0,nrow=10,ncol=10);
#' ppiA.m[1,] <- c(0,1,1,1,1);
#' for(r in 2:nrow(ppiA.m)){
#'   ppiA.m[r,1] <- 1;
#' }
#' rownames(ppiA.m) <- paste("G",1:10,sep="");
#' colnames(ppiA.m) <- paste("G",1:10,sep="");
#' 
#' ### define a positively valued expression matrix (20 genes x 10 samples)
#' exp.m <- matrix(rpois(20*10,8),nrow=20,ncol=10);
#' colnames(exp.m) <- paste("S",1:10,sep="");
#' rownames(exp.m) <- paste("G",1:20,sep="");
#' 
#' ### run integration function
#' Integration.l <- DoIntegPPI(exp.m,ppiA.m);
#' print(dim(Integration.l$expMC));
#' print(dim(Integration.l$adjMC));
#' 
#' @import Biobase
#' @import Matrix
#' @import scater
#' @importFrom BiocGenerics estimateSizeFactors
#' @importFrom scater normalize
#' @importFrom scater librarySizeFactors
#' @importFrom igraph graph.adjacency
#' @importFrom igraph clusters
#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedArray isEmpty
#' @importFrom Matrix colSums
#' @importFrom methods as
#' @import monocle
#' @export
#'     
DoIntegPPI <- function(exp.m, 
                       ppiA.m,
                       log_trans = FALSE)
{
    # set input data matrix class
    data.class <- class(exp.m)[1]
    
    ## set row_names & col_names
    row_names <- rownames(exp.m)
    col_names <- colnames(exp.m)
    
    if (is.null(col_names)) {
        warning(paste0("No cell names been specified, forcedly set as cell_1 to cell_", 
                       ncol(exp.m), "."))
        colnames(exp.m) <- paste0("cell_", 1:ncol(exp.m))
        col_names <- paste0("cell_", 1:ncol(exp.m))
    }
    
    # set common gene IDs
    commonEID.v <- intersect(rownames(ppiA.m), row_names)
    
    # check the consistency of gene IDs
    if (DelayedArray::isEmpty(commonEID.v) == TRUE) {
        stop("scRNA-seq data should have the same gene identifier with the network!")
    }
    
    # get maximum connected network
    match(commonEID.v,rownames(ppiA.m)) -> map2.idx
    adj.m <- ppiA.m[map2.idx, map2.idx]
    
    gr.o <- igraph::graph.adjacency(adj.m,mode="undirected")
    comp.l <- igraph::clusters(gr.o)
    cd.v <- summary(factor(comp.l$member))
    mcID <- as.numeric(names(cd.v)[which.max(cd.v)])
    maxc.idx <- which(comp.l$member==mcID)
    adjMC.m <- adj.m[maxc.idx, maxc.idx]
    
    # select log-normalization methods based on data class
    if (data.class == "SingleCellExperiment") {
        sizeFactors(exp.m) <- scater::librarySizeFactors(exp.m)
        exp.m <- scater::normalize(exp.m, log_exprs_offset = 1.1)
        data.m <- Matrix::as.matrix(SummarizedExperiment::assay(exp.m, i = "logcounts"))
    }else if (data.class == "CellDataSet") {
        exp.m <- BiocGenerics::estimateSizeFactors(exp.m)
        data.m <- Matrix::as.matrix(t(t(Biobase::exprs(exp.m)) / 
                                          Biobase::pData(exp.m)[, 'Size_Factor']))
        data.m <- log2(data.m + 1.1)
    }else{
        
        if (log_trans) {
            
            if (data.class == "dgCMatrix") {
                temp.m <- exp.m
                rm(exp.m)
                gc(verbose = FALSE)
            }else{
                ### set data matrix as dgCMatrix
                temp.m <- as(exp.m, "dgCMatrix")
                rm(exp.m)
                gc(verbose = FALSE)
            }
            
            ### column normalization
            # pre-compute maximum library size
            lib_max <- max(Matrix::colSums(temp.m))
            temp.m <- wordspace::normalize.cols(temp.m, method = "manhattan")
            # save additional dgCMatrix
            dgC.m <- temp.m
            dgC.m@x <- log2((dgC.m@x * lib_max) + 1)
            temp.m@x <- log2((temp.m@x * lib_max) + 1.1)
            
            temp.dgT <- as(temp.m, "dgTMatrix")
            data.m <- as.matrix(temp.dgT)
            data.m[data.m == 0] <- log2(1.1)
            rm(temp.m, temp.dgT)
            gc(verbose = FALSE)
        }else{
            data.m <- exp.m
            dgC.m <- data.m - min(data.m)
            dgC.m <- as(dgC.m, "dgCMatrix")
            rm(exp.m)
            gc(verbose = FALSE)
        }
        
    }
    
    if (min(data.m) == 0) {
        stop(paste0("Input matrix must have non-zero minimal value, please set", 
                    " 'log_trans' = TRUE!"))
    }
    
    match(rownames(adjMC.m), row_names) -> map1.idx
    expMC.m <- data.m[map1.idx ,]
    dgC.m <- dgC.m[map1.idx ,]
    
    if (!identical(rownames(dgC.m), rownames(expMC.m))) {
        stop("Non identical!!!")
    }
    
    # compute node degree
    gr.o <- igraph::graph.adjacency(adjMC.m, mode="undirected")
    degree.v <- igraph::degree(gr.o)
    degree.v <- degree.v[rownames(expMC.m)]
    
    if (data.class == "SingleCellExperiment") {
        return(list(data.sce = exp.m, expMC = expMC.m, adjMC = adjMC.m, 
                    data = data.m, degree.v = degree.v, dgC.m = dgC.m))
    }else if (data.class == "CellDataSet") {
        return(list(data.cds = exp.m, expMC = expMC.m, adjMC = adjMC.m, 
                    data = data.m, degree.v = degree.v, dgC.m = dgC.m))
    }else{
        return(list(expMC = expMC.m, adjMC = adjMC.m, 
                    data = data.m, degree.v = degree.v, dgC.m = dgC.m))
    }
}
