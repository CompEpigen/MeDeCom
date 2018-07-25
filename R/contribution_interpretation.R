########################################################################################################################################
## contribution_interpretation.R
## -------------------------------------------------------------------------------------------------------------------------------------
## This scripts contains functions for assocating LMC contributions with sample groupings and assessing statistical significance
## for the differences.
######################################################################################################################################## 

######################################################################################################################################## 
## GLOBALS

######################################################################################################################################## 
## FUNCTIONS

#' run.trait.association
#' 
#' Computes test statistics for all possible group assignments of samples defined in \code{medecom.set} and \code{rnb.set} and stores
#' heatmaps of p-values on the given location
#' 
#' @param medecom.set An object of type \code{\link{MeDeComSet}} as the result of \code{\link{runMeDeCom}} containing LMCs and their
#'                     proportions in the samples. The Set can contain multiple runs for different values of K and lambda.
#' @param rnb.set An object of type \code{\link[RnBeads]{RnBSet}} containing methylation data and metadata for the same samples for which 
#'                 \code{medecom.set} was computed.
#' @param test.fun Test statistic used to compute p-values of differences between LMC contributions in pairwise sample comparisons.
#'                  Defaults to \code{t.test}.
#' @param plot.path Path to store the p-value heatmaps.
#' @param figure.format Character describing the format in which plots should be stored on disk. Either \code{"pdf"} or \code{"png"}.
#' 
#' @details This function creates a new folder names \code{pdfs} at the location given by \code{plot.path} and stores a heatmap for
#'           all possible Ks and lambdas defined in \code{medecom.set}. The p-values are produced by comparing the LMC contributions
#'           in all sample comparisons defined by \code{\link[RnBeads]{rnb.sample.groups}} on \code{rnb.set}. The employed test statistic for
#'           pariwise comparison can be specified by \code{test.fun}, for groups defining more than one group \code{\link{kruskal.test}}
#'           is employed. P-values lower than 0.01 are added to the heatmap.
#'           
#' @author Michael Scherer
#' 
#' @export
run.trait.association <- function(medecom.set,rnb.set,test.fun=t.test,plot.path=getwd(),figure.format="pdf"){
  Ks <- medecom.set@parameters$Ks
  lambdas <- medecom.set@parameters$lambdas
  if(!file.exists(plot.path)){
    stop(paste("Location",plot.path,"does not exist."))
  }
  if(!figure.format %in% c("pdf","png")){
    stop(paste("Invalid value for figure.format, needs to be 'pdf' of 'png'"))
  }
  plot.path <- file.path(plot.path,figure.format)
  if(!file.exists(plot.path)){
    dir.create(plot.path)
  }
  for(K in Ks){
    for(lambda in lambdas){
      p.vals <- link.to.traits(medecom.set=medecom.set,K=K,lambda=lambda,rnb.set=rnb.set,test.fun=test.fun)
      plot <- plot.p.val.heatmap(p.vals)
      fname <- paste0("trait_association_heatmap_K",K,"_lambda",lambda,figure.format)
      ggsave(filename=file.path(plot.path,fname),plot=plot,device = figure.format)
    }
  }
}

#' link.to.traits
#' 
#' This routine performs a statistical test to determine if the difference in LMC contributions is different for the sample groups
#' defined in \code{rnb.set}.
#' 
#' @param medecom.set An object of type \code{\link{MeDeComSet}} as the result of \code{\link{runMeDeCom}} containing LMCs and their
#'                     proportions in the samples. The Set can contain multiple runs for different values of K and lambda.
#' @param K The K parameter, determining the number of LMCs to extract from \code{medecom.set}.
#' @param K The lambda parameter, determining the regularization used for the LMCs in \code{medecom.set}.
#' @param rnb.set An object of type \code{\link{RnBSet}} containing methylation data and metadata for the same samples for which 
#'                 \code{medecom.set} was computed.
#' @param test.fun Test statistic used to compute p-values of differences between LMC contributions in pairwise sample comparisons.
#'                  Defaults to \code{t.test}.
#'                  
#' @return A list with an element for each sample grouping defined by \code{\link{rnb.sample.groups}}.
#' 
#' @details Each element in the returned list is of left \code{K}, displaying the p-value of the statistical assocation of the
#'           contributions of the corresponding LMC to the sample grouping. The p-values are produced by comparing the LMC 
#'           contributions in all sample comparisons defined by \code{\link{rnb.sample.groups}} on \code{rnb.set}. The employed test
#'           statistic for pariwise comparison can be specified by \code{test.fun}, for groups defining more than one group
#'           \code{\link{kruskal.test}} is employed.
#' @author Michael Scherer
#' 
#' @noRd

link.to.traits <- function(medecom.set,K,lambda,rnb.set,test.fun=t.test){
  require("RnBeads")
  sample.grps <- rnb.sample.groups(rnb.set)
  props <- getProportions(medecom.set,K=K,lambda=lambda)
  res <- list()
  names.grps <- names(sample.grps)
  for(i in 1:length(sample.grps)){
    grp <- sample.grps[[i]]
    names.traits <- names(grp)
    if(length(grp)==2){
      p.vals <- apply(props,1,function(x){
        test.fun(x[grp[[1]]],x[grp[[2]]])$p.val
      })
    }else if(length(grp)>2){
      vec <- rep(NA,length(samples(rnb.set)))
      for(name in names.traits){
        vec[grp[[name]]] <- name
      }
      vec <- as.factor(vec)
      p.vals <- apply(props,1,function(x){
        kruskal.test(x,vec)$p.val
      })
    }else{
      p.vals <- NA
    }
    res[[names.grps[i]]] <- p.vals
  }
  return(res)
}

#' plot.p.val.heatmap
#' 
#' Produces a heatmap for the p-values produced in \code{\link{link.to.traits}}.
#' 
#' @param trait.res A list as the result of \code{\link{link.to.traits}}.
#' 
#' @return A plot object displaying the heatmap of the input p-values. P-values lower than 0.01 are added to the plot as numbers,
#'          otherwise the decadic logarithm of the p-value resembles the shading of the tiles.
#'          
#' @author Michael Scherer
#' @noRd

plot.p.val.heatmap <- function(trait.res){
  require("ggplot2")
  require("reshape2")
  to.plot <- as.data.frame(trait.res)
  to.plot$LMC <- row.names(to.plot)
  to.plot <- melt(to.plot)
  colnames(to.plot)[2:3] <- c("Trait","PValue")
  to.plot$LogPValue <- log(to.plot$PValue)
  plot <- ggplot(to.plot,aes(x=LMC,y=Trait,fill=LogPValue))+geom_tile()+theme_bw()+scale_fill_gradient(low="red",high = "white")+
    geom_text(aes(label=ifelse(round(PValue,2)< 0.01,format(PValue,digits = 2),"")))
  return(plot)
}
