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
#' @param Integration.l
#' A list object from \code{DoIntegPPI} function.
#' 
#' @param data.m
#' A scRNA-Seq data matrix with rows labeling genes and columns 
#' labeling single cells. And it can be either a log-transformed data
#' matrix with minimal value around 0.1 (recommended), or an 
#' nonlog-transformed data matrix with minimal value 0.
#' 
#' @param ppiA.m
#' The adjacency matrix of a user-given PPI network with rownames and 
#' colnames labeling genes (same gene identifier as in \code{exp.m})
#' 
#' @param log_trans
#' A logical. Whether to do log-transformation on the input data
#' matrix or not. Default is FALSE
#' 
#' @param parallelMode
#' A logical. Indicating whether or not to run CCAT in parallel. 
#' It will be disable if datasets are < 5,000 cells. Parallel mode 
#' uses a subsampling approach to reduce runtime. Default is FALSE
#' 
#' @param mcores
#' A integer. Indicating the number of cores to use when parallelMode = TRUE
#' 
#' @param subsamplesize
#' A integer. Indicating the number of cells to subsample when parallelMode = TRUE
#' 
#' @return A list incorporates the input list and CCAT velues or CCAT values 
#' itself, depending on the input object(s):
#' 
#' @return CCAT
#' The estimated signaling entropy rate using Pearson correlation coefficient
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
#' Teschendorff Andrew E., Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Teschendorff Andrew E., Banerji CR, Severini S, Kuehn R, Sollich P. 
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
#' ### load example data & network matrix
#' data(Example.m)
#' data(net13Jun12.m)
#' 
#' ### integrate expr matrix and PPI network
#' Integration.l <- DoIntegPPI(exp.m = Example.m, ppiA.m = net13Jun12.m)
#' 
#' ### estimate SR with PCC
#' ### get it with the integration list
#' Integration.l <- CCAT(Integration.l)
#' 
#' ### or get CCAT directly from data matrix
#' CCAT.v <- CCAT(data.m = Example.m, ppiA.m = net13Jun12.m)
#' 
#' 
#' @import Biobase
#' @import Matrix
#' @import wordspace
#' @importFrom wordspace normalize.cols
#' @importFrom igraph graph.adjacency
#' @importFrom igraph clusters
#' @importFrom DelayedArray isEmpty
#' @importFrom ccaPP corPearson
#' @importFrom Matrix colSums
#' @importFrom methods as
#' @export
#'     
CCAT <- function(Integration.l = NULL,
                 data.m = NULL, 
                 ppiA.m = NULL,
                 log_trans = FALSE,
                 parallelMode = FALSE,
                 mcores = 1,
                 subsamplesize = 1000)
{
  
  if (is.null(Integration.l)) {
    
    # set matrix as dgCMatrix
    data.m <- as(data.m, "dgCMatrix")
    gc(verbose = FALSE)
    
    # get rownames/colnames of matrix
    row_names <- rownames(data.m)
    col_names <- colnames(data.m)
    
    if (is.null(col_names)) {
      warning(paste0("No cell names been specified, forcedly set as cell_1 to cell_", 
                     ncol(data.m), "."))
      colnames(data.m) <- paste0("cell_", 1:ncol(data.m))
      col_names <- paste0("cell_", 1:ncol(data.m))
    }
    
    # set common gene IDs
    commonEID.v <- intersect(rownames(ppiA.m), row_names)
    
    # check the consistency of gene IDs
    if (DelayedArray::isEmpty(commonEID.v) == TRUE) {
      stop("scRNA-seq data should have the same gene identifier with the network!")
    }
    
    # set parallelMode = FALSE if cell number is small
    if (parallelMode) {
      if (length(col_names) < 5000) {
        parallelMode <- FALSE
        warning(paste0("Number of cells is lower than 5K, ",
                       "'parallelMode' is disabled!"))
      }
    }
    
    # subset PPI network
    match(commonEID.v, rownames(ppiA.m)) -> map2.idx
    adj.m <- ppiA.m[map2.idx,map2.idx]
    
    # get maximum connected sub-network
    gr.o <- igraph::graph.adjacency(adj.m, mode="undirected")
    comp.l <- igraph::clusters(gr.o)
    cd.v <- summary(factor(comp.l$member))
    mcID <- as.numeric(names(cd.v)[which.max(cd.v)])
    maxc.idx <- which(comp.l$member==mcID)
    adjMC.m <- adj.m[maxc.idx, maxc.idx]
    
    if (log_trans) {
      # pre-compute maximum library size
      lib_max <- max(Matrix::colSums(data.m))
      data.m <- wordspace::normalize.cols(data.m, method = "manhattan")
      data.m@x <- log2((data.m@x * lib_max) + 1)
    }
    
    # subset data matrix
    match(rownames(adjMC.m), row_names) -> map1.idx
    expMC.m <- data.m[map1.idx ,]
    
    # compute node degree
    gr.o <- igraph::graph.adjacency(adjMC.m, mode="undirected")
    degree.v <- igraph::degree(gr.o)
    degree.v <- degree.v[rownames(expMC.m)]
    
    # subsample data set
    if (parallelMode) {
      chunk <- round(length(col_names)/subsamplesize)
      subsamples <- split(1:length(col_names), sample(factor(1:length(col_names) %% chunk)))
    }else if (length(col_names) > 10000) {
      chunk <- round(length(col_names)/1000)
      subsamples <- split(1:length(col_names), sample(factor(1:length(col_names) %% chunk)))
    }else{
      subsamples <- NULL
    }
    
    CCAT <- Comp_PCCSR(expMC.m = expMC.m,
                       col_names = col_names,
                       degree.v = degree.v,
                       subsamples = subsamples,
                       parallelMode = parallelMode,
                       mcores = mcores)
    
    return(CCAT = CCAT)
    
  }else{
    
    degree.v <- Integration.l$degree.v
    expMC.m <- Integration.l$dgC.m
    
    # get rownames/colnames of matrix
    col_names <- colnames(expMC.m)
    
    # set parallelMode = FALSE if cell number is small
    if (parallelMode) {
      if (length(col_names) < 5000) {
        parallelMode <- FALSE
        warning(paste0("Number of cells is lower than 5K, ",
                       "'parallelMode' is disabled!"))
      }
    }
    
    # subsample data set
    if (parallelMode) {
      chunk <- round(length(col_names)/subsamplesize)
      subsamples <- split(1:length(col_names), sample(factor(1:length(col_names) %% chunk)))
    }else if (length(col_names) > 10000) {
      chunk <- round(length(col_names)/1000)
      subsamples <- split(1:length(col_names), sample(factor(1:length(col_names) %% chunk)))
    }else{
      subsamples <- NULL
    }
    
    CCAT <- Comp_PCCSR(expMC.m = expMC.m,
                       col_names = col_names,
                       degree.v = degree.v,
                       subsamples = subsamples,
                       parallelMode = parallelMode,
                       mcores = mcores)
    
    Integration.l$CCAT <- CCAT
    return(Integration.l)
    
  }
  
}

Comp_PCCSR <- function(expMC.m = NULL,
                       col_names = NULL,
                       degree.v = NULL,
                       subsamples = NULL,
                       parallelMode = FALSE,
                       mcores = 1){
  
  if (is.null(subsamples)) {
    
    PCCSR <- apply(expMC.m, 2, function(x) ccaPP::corPearson(x, degree.v))
    names(PCCSR) <- col_names
    return(PCCSR = PCCSR)
    
  }else{
    
    if (length(subsamples) == 1) {
      parallelMode <- FALSE
      warning(paste0("Parallel Mode disabled, Run the computation on single core!\n"
                     ,"If you want to run on multiple cores, try to reduce subsamplesize."))
    }
    
    if (parallelMode) {
      
      parallel.l <- parallel::mclapply(subsamples, mc.cores = mcores, function(subsample){
        
        exp.m <- expMC.m[, subsample]
        PCCSR <- apply(exp.m, 2, function(x) ccaPP::corPearson(x, degree.v))
        names(PCCSR) <- colnames(exp.m)
        return(list(PCCSR = PCCSR))
        
      })
      PCCSR <- do.call(c, c(lapply(parallel.l, function(x) x$PCCSR), use.names = FALSE))
      PCCSR <- PCCSR[col_names]
      return(PCCSR = PCCSR)
      
    }else{
      
      PCCSR <- vector(length = ncol(expMC.m))
      for (i in 1:length(subsamples)) {
        gc(verbose = FALSE)
        subsample <- subsamples[[i]]
        exp.m <- expMC.m[, subsample]
        PCCSR.tmp <- apply(exp.m, 2, function(x) ccaPP::corPearson(x, degree.v))
        PCCSR[subsample] <- PCCSR.tmp
      }
      names(PCCSR) <- col_names
      return(PCCSR = PCCSR)
      
    }
    
  }
  
}




