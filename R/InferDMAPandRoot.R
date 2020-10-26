#' @title 
#' Infers Diffusion Map, Root-state and root-cell
#' 
#' @aliases InferDMAPandRoot
#'  
#' @description 
#' This function constructs the diffusion map, Markov transition matrix and uses the provided potency 
#' estimates to infer the root-state and root-cell.
#' 
#' @param pot.v
#' A vector of potency estimates for all cells, e.g. a vector of SR or CCAT values.
#' 
#' @param exp.m
#' The normalized and log-transformed scRNA-Seq data matrix
#'
#' @param avth
#' The threshold on the average expression to use for selecting genes in the diffusion map. Default values is 1.
#'
#' @param sdth
#' The threshold on the standard deviation to use for selecting genes in the diffusion map. Default value is 0.25.
#' 
#' @param kDMAP
#' An integer specifying the input argument to the diffusion map, representing the number of
#' nearest neighbors to consider. By default this value is 30.
#'
#' @param pctop
#' The proportion of top-ranked cells (with cells ranked by SR or CCAT) to consider when inferring teh root-state. By default this values is 0.05.
#' 
#' @return A list with the following elements:
#'
#' @return dmap
#' The diffusion map object, as given by the `destiny` package.
#' 
#' @return dc
#' The inferred diffusion map component matrix, with rows labeling cells and columns labeling ranked diffusion components.
#' 
#' @return transM
#' The Markov transition matrix associated with the inferred diffusion map. Rows and columns label cells.
#' @return root
#' The column index of the root-cell.
#' 
#' 
#' @references 
#' Teschendorff AE, Tariq Enver. 
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
#' @import destiny
#' @import igraph
#' @import marray
#' 
#' @export
#'
#' 
InferDMAPandRoot <- function(pot.v,exp.m,avth=1,sdth=0.25, kDMAP=30, pctop=0.05){

  spot.o <- sort(pot.v,decreasing=TRUE,index.return=TRUE);
  ntop <- floor(0.05*length(pot.v));
  print(paste("The root-state will be inferred from the top ",ntop," cells with highest potency",sep=""));
  if(ntop < 10){
      stop("We recommend increasing pctop to ensure a larger number of cells to select root-state");
  }
  cand.idx <- spot.o$ix[1:ntop];

  sd.v <- apply(exp.m,1,sd);
  mean.v <- apply(exp.m,1,mean);
  selG.idx <- intersect(which(mean.v>avth),which(sd.v>sdth));
  print(paste("The diffusion map will be constructed from ",length(selG.idx)," variable genes",sep=""));
  tmp.m <- exp.m[selG.idx,];

  dmap.o <- DiffusionMap(data=t(tmp.m),k=kDMAP,verbose=TRUE);
  dc.m <- eigenvectors(dmap.o);
  transM.m <- as.matrix(dmap.o@transitions); ## symmetric matrix
  colnames(transM.m) <- paste("Cell",1:length(pot.v),sep="");
  rownames(transM.m) <- paste("Cell",1:length(pot.v),sep="");

  gr.o <- graph_from_adjacency_matrix(transM.m[cand.idx,cand.idx], mode = "undirected", weighted =TRUE,  diag = TRUE);
  nsteps <- floor(0.25*length(cand.idx));
  print(paste("Running walk-trap community detection with max.step number= ",nsteps,sep=""));
  walk.o <- cluster_walktrap(gr.o, steps = nsteps);
  distr.v <- summary(factor(walk.o$member));
  max.idx <- which.max(distr.v);
  maxclID <- max.idx;
  selcand.idx <- cand.idx[which(walk.o$member==max.idx)];
  print(paste("Walk-trap identified a root-state consisting of ",length(selcand.idx)," cells",sep=""));
    
  pos.m <- dc.m[selcand.idx,];
  med.m <- apply(pos.m,2,median);
  mad.v <- apply(abs(pos.m-med.m),1,mean)
  root.idx <- selcand.idx[which.min(mad.v)];

  q.v <- quantile(pot.v,probs=seq(0.05,0.95,0.1));
  macolor.v <- maPalette(low="skyblue",high="blue",k=10);
  color.v <- rep(NA,length(pot.v));
  color.v[which(pot.v > q.v[length(q.v)])] <- "black";
  for(i in length(q.v):1){
   color.v[which(pot.v <= q.v[i])] <- macolor.v[i]
  }
    
  return(list(dmap=dmap.o,dc=dc.m,transM=transM.m,root=root.idx,color=color.v));
} ## EOF
