#' @title 
#' Computes potency estimates of single-cells 
#' using the signaling entropy rate
#' 
#' @aliases CompSRana
#'  
#' @description 
#' This is the main user function for computing signaling entropy of 
#' single cells. It takes as input the gene expression profile of 
#' single cells and the adjacency matrix of a connected network. These 
#' inputs will be typically the output of the \code{DoIntegPPI} function.
#' 
#' @param integ.l
#' The output from \code{DoIntegPPI} function.
#' 
#' @param local
#' A logical (default is FALSE). If TRUE, function returns the normalized 
#' local signaling entropies of each gene in the network alongside the unnormalized
#' local entropies. If FALSE, only unnormalized local entropies are returned.
#' 
#' @param mc.cores
#' The number of cores to use, i.e. at most how many child processes will 
#' be run simultaneously. The option is initialized from environment variable 
#' MC_CORES if set. Must be at least one, and parallelization requires at 
#' least two cores.
#' 
#' @return A list with four elements:
#' 
#' @return SR
#' The global signaling entropy rate. It is normalized by the 
#' maximum rate, hence a value between 0 and 1
#' 
#' @return inv
#' The stationary distribution of every sample
#' 
#' @return locS
#' The unnormalised local entropies of each gene in every cell
#' 
#' @return nlocS
#' The normalised local entropies of each gene, so that each value is 
#' between 0 and 1
#' 
#' @references 
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
#' @import parallel
#' 
#' @export
#'
CompSRana <- function(integ.l, local = FALSE, mc.cores=1){

   integ.l$adjMC <- as(integ.l$adjMC, "dgCMatrix");

 ### compute maxSR for SR normalization
    maxSR <- CompMaxSR(integ.l);

    idx.l <- as.list(seq_len(ncol(integ.l$expMC)));
    out.l <- mclapply(idx.l, CompSRanaPRL, 
                      exp.m=integ.l$expMC, 
                      adj.m=integ.l$adjMC,
                      local=local,
                      maxSR=maxSR,
                      mc.cores=mc.cores)
    SR.v <- sapply(out.l, function(v) return(v[[1]]))
    invP.v <- sapply(out.l, function(v) return(v[[2]]))
    S.v <- sapply(out.l, function(v) return(v[[3]]))
    NS.v <- sapply(out.l, function(v) return(v[[4]]))
    
    return(list(SR=SR.v,inv=invP.v,locS=S.v,nlocS=NS.v));
}

#### Auxilliary functions
#' @import igraph
#' 
CompMaxSR <- function(integ.l){
    
    adj.m <- integ.l$adjMC
    if (!inherits(adj.m, "Matrix")){
        adj.m <- Matrix::Matrix(adj.m, sparse = TRUE)
    }
    # find right eigenvector of adjacency matrix
    fa <- function(x,extra=NULL) {
        as.vector(adj.m %*% x)
    }
    ap.o <- igraph::arpack(fa,options=list(n=nrow(adj.m),nev=1,which="LM"), sym=TRUE)
    v <- ap.o$vectors
    lambda <- ap.o$values
    
    # maximum entropy
    maxSR <- log(lambda)
    
    return(maxSR);
}

CompSRanaPRL <- function(idx,
                         exp.m,
                         adj.m,
                         local=TRUE,
                         maxSR=NULL)
{
    n <- nrow(adj.m);

    e <- exp.m[,idx];
    if (inherits(e, "Matrix")) {e <- as.vector(e);}
    s <- as.vector(adj.m %*% e);


    invP.v <- e*s;
    nf <- sum(invP.v);
    invP.v <- invP.v/nf;
    B <- adj.m;

    B@x <- B@x * rep(e,diff(B@p));
    Blog <- B;
    nz <- Blog@x > 0   ;
    Blog@x[nz]  <- Blog@x[nz] * log(Blog@x[nz]);
    Blog@x[!nz] <- 0 ;

    rB    <- Matrix::rowSums(B) ;
    rBlog <- Matrix::rowSums(Blog) ;

    S.v <- numeric(n);
    good <- which(s>0);
    S.v[good] <- log(s[good]) - rBlog[good]/s[good] ;

    NS.v <- NULL;
    if(local){
        deg <- Matrix::rowSums(adj.m > 0);
        NS.v <- numeric(n);
        g2 <- which(deg > 1);
        NS.v[g2] <- S.v[g2] / log(deg[g2]);
    
    }
    SR <- sum(invP.v * S.v);
    if(!is.null(maxSR)){ SR <- SR/maxSR ;}

    return(list(sr=SR,inv=invP.v,locS=S.v,nlocS=NS.v));
}

CompNS <- function(p.v){
    
    tmp.idx <- which(p.v>0);
    if(length(tmp.idx)>1){
        NLS <- -sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
    }
    else {
        # one degree nodes have zero entropy, avoid singularity.
        NLS <- 0;
    }
    return(NLS);
}

CompS <- function(p.v){
    
    tmp.idx <- which(p.v>0);
    LS <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
    return(LS);
}
