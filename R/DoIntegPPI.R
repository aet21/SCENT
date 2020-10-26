#' @title 
#' Integration of gene expression matrix and PPI network
#' 
#' @aliases DoIntegPPI
#'  
#' @description 
#' This function finds the common genes between the scRNA-Seq data matrix 
#' and the genes present in the PPI network, and constructs the maximally 
#' connected subnetwork and associated expression matrix for the computation 
#' of signaling entropy.
#' 
#' @param exp.m
#' The scRNA-Seq data matrix normalized for library size and log2-transformed
#' with a pseudocount of 1.1
#' 
#' @param ppiA.m
#' The adjacency matrix of a user-given PPI network with rownames and 
#' colnames labeling genes (same gene identifier as in \code{exp.m})
#' 
#' 
#' @return A list of two or four objects:
#' 
#' @return expMC
#' Reduced expression matrix with genes in the maximally connected subnetwork
#' 
#' @return adjMC
#' Adjacency matrix of the maximally connected subnetwork
#' 
#' 
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
#'
#' 
#' @import Matrix
#' @importFrom igraph graph.adjacency
#' @importFrom igraph clusters
#' @export
#'     
DoIntegPPI <- function(exp.m, 
                       ppiA.m)
{

   if(max(exp.m) > 100){ ### check if data has been log2-transformed or not and if not, then log2-transform with a pseudcount of 1.1, as for SR computation later we are not allowed 0's.
        exp.m <- log2(exp.m+1.1);
    }

    # set common gene IDs
    commonEID.v <- intersect(rownames(ppiA.m),rownames(exp.m));
    
    # check the consistency of gene IDs and size of overlap
    if ( length(commonEID.v) < 5000 ){
      stop(paste("The overlap of common genes between PPI and expression matrix is only ",length(commonEID.v)," so check that gene identifiers are correct. We don't recommend running CCAT on less than 5000 overlapping genes.",sep=""));
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
    
    match(rownames(adjMC.m),rownames(exp.m)) -> map1.idx
    expMC.m <- exp.m[map1.idx ,]

    return(list(expMC=expMC.m,adjMC=adjMC.m));        

} ### EOF
