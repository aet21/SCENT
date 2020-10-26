#' @title 
#' Infer distinct potency states from the single cell potency estimates
#' 
#' @aliases InferPotencyStates
#'  
#' @description 
#' This function infers the discrete potency states of single cells and 
#' its distribution across the single cell population.
#' 
#' @param potest.v
#' A vector of potency estimates for all cells, e.g. a vector of SR or CCAT values.
#' 
#' @param type
#' The type of potency estimate used (SR or CCAT).
#' 
#' @param pheno.v
#' A phenotype vector for the single cells.
#' 
#' @param diffvar
#' A logical. Default is TRUE.
#' Specifies whether the Gaussian mixture model to be fit assumes components 
#' to have different (default) or equal variance.
#' In the latter case, use *modelNames = c("E")*.
#' 
#' @param maxPS
#' Maximum number of potency states, when inferring discrete potency 
#' states of single cells. Default value is 5.
#' 
#' @return A list with the following elements:
#' 
#' @return potS
#' Inferred discrete potency states for each single cell. It is indexed so 
#' that the index increases as the signaling entropy of the state decreases
#' 
#' @return distr
#' If phenotype information provided, it will be a table giving the 
#' distribution of single-cells across potency states and phenotypes
#' 
#' @return prob
#' Table giving the probabilities of each potency state per phenotype value
#' 
#' @return het
#' The normalised Shannon Index of potency per phenotype value
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
#' 
#' @examples 
#' 
#' 
#' @import mclust
#' 
#' 
#' @export
#'
#' 
InferPotencyStates <- function(potest.v, type=c("SR","CCAT"), pheno.v = NULL, diffvar = TRUE,maxPS = 5){


    if(type=="SR"){
    sr.v <- potest.v;
    ### fit Gaussian Mixture Model for potency inference
    print("Fit Gaussian Mixture Model to logit-transformed SR values")
    logitSR.v <- log2(sr.v / (1 - sr.v))
    if(diffvar == TRUE){ 
        ## default assumes different variance for clusters
        mcl.o <- Mclust(logitSR.v, G = seq_len(maxPS))
    }
    else {
        mcl.o <- Mclust(logitSR.v, G = seq_len(maxPS), modelNames = c("E"))
    }
    mu.v <- mcl.o$param$mean
    sd.v <- sqrt(mcl.o$param$variance$sigmasq)
    avPS.v <- (2^mu.v)/(1+2^mu.v)
    }
    else if (type=="CCAT"){
      ccat.v <- potest.v;
      print("Fit Gaussian Mixture Model to Z-transformed CCAT values")    
      zccat.v <- log2((1+ccat.v)/(1-ccat.v));
      if(diffvar == TRUE){ 
        ## default assumes different variance for clusters
         mcl.o <- Mclust(zccat.v, G = seq_len(maxPS))
      }
      else {
         mcl.o <- Mclust(zccat.v, G = seq_len(maxPS), modelNames = c("E"))
      }
      mu.v <- mcl.o$param$mean
      sd.v <- sqrt(mcl.o$param$variance$sigmasq)
      avPS.v <- (2^mu.v - 1)/(2^mu.v + 1);
    }
    potS.v <- mcl.o$class
    nPS <- length(levels(as.factor(potS.v)))
    print(paste("Identified ",nPS," potency states",sep=""))
    for (i in seq_len(nPS)) {
        names(potS.v[which(potS.v == i)]) <- rep(paste("PS",i,sep=""), 
                                                 times = length(which(potS.v == i)))
    }
    
  
    savPS.s <- sort(avPS.v, decreasing=TRUE, index.return=TRUE)
    spsSid.v <- savPS.s$ix
    ordpotS.v <- match(potS.v,spsSid.v)
    
    if(!is.null(pheno.v)){
        nPH <- length(levels(as.factor(pheno.v)))
        distPSph.m <- table(pheno.v,ordpotS.v)
        print("Compute Shannon (Heterogeneity) Index for each Phenotype class")
        probPSph.m <- distPSph.m/apply(distPSph.m,1,sum)
        hetPS.v <- vector()
        for(ph in seq_len(nPH)){
            prob.v <- probPSph.m[ph,]
            sel.idx <- which(prob.v >0)
            hetPS.v[ph] <- 
                - sum(prob.v[sel.idx]*log(prob.v[sel.idx]))/log(nPS)
        }
        names(hetPS.v) <- rownames(probPSph.m)
        print("Done")
    }
    else {
        distPSph.m=NULL
        probPSph.m=NULL
        hetPS.v=NULL
    }
    

    

    
    return(list(potS=ordpotS.v,distr=distPSph.m,prob=probPSph.m,het=hetPS.v));
}
