###############################################################################################################
# LMC_interpreation.R
# functions to interpret LMCs and to assign potential functions to the LMCs
###############################################################################################################

#' lmc.lola.enrichment
#' 
#' This routine computes LOLA enrichment results for LMC-specifically hypomethylated sites. 
#' 
#' @param medecom.result An object of type \code{\link{MeDeComSet-class}} or the location of an .RData file, where
#'                 such an object is stored.
#' @param annotation.filter A numeric vector specifyfing the sites that have been removed from \code{rnb.set} in a
#'                 preprocessing step (e.g. coverage filtering) or a path to an .RData file.
#' @param rnb.set The original \code{\link[RnBeads]{RnBSet-class}} object containing methylation, sample meta and annotation
#'                 information or a path to a directory stored by \code{\link[RnBeads]{save.rnb.set}}.
#' @param K The number of LMCs specified for the MeDeCom run.
#' @param lambda The lambda parameter selected.
#' @param cg_subset The index of the selection strategy employed (e.g. most variable CpGs).
#' @param diff.threshold The difference cutoff between median methylation in the remaining LMCs and the LMC of interest
#'                  used to call a CpG differentially methylated. The higher this value, the more conservative the
#'                  selection.
#' @param region.type Region type used to annotate CpGs to potentially regulatory regions (see \url{https://rnbeads.org/regions.html})
#'                  for a list of available region types.
#' @param temp.dir Path to a directory used to store temporary files.
#' @return A list with K elements. One element is the result of the corresponding LMC-specific hypomethylated CpG sites.
#' @export
#' @details This function employs LOLA on the CpG sites that are LMC-specifically hypomethylated, after annotating 
#'                  the sites to the closest region defined by \code{region.type}. The sites are selected by computing
#'                  the median methylation value of the other LMCs and then selecting those sites that are more than 
#'                  \code{diff.threshold} away from the median in the LMC of interest. This is done for all LMCs from
#'                  1 to K.
#' @author Michael Scherer                                

lmc.lola.enrichment <- function(medecom.result,
                                annotation.filter=NULL,
                                rnb.set,
                                K,
                                lambda,
                                cg_subset=NULL,
                                diff.threshold=0.5,
                                region.type="ensembleRegBuildBPall",
                                temp.dir=tempdir()){
  require("RnBeads")
  require("LOLA")
  if(is.character(medecom.result)){
    new.envi <- new.env()
    load(medecom.result,envir = new.envi)
    medecom.result <- get(ls(envir = new.envi),envir = new.envi)
  }
  if(is.character(annotation.filter)){
    new.envi <- new.env()
    load(annotation.filter,envir = new.envi)
    annotation.filter <- get(ls(envir = new.envi),envir = new.envi)
    
  }
  if(is.character(rnb.set)){
    options(fftempdir=temp.dir)
    rnb.set <- load.rnb.set(rnb.set)
  }
  if(!region.type %in% rnb.region.types()){
    rnb.load.annotation.from.db(region.type,assembly=assembly(rnb.set))
  }
  agg.type <- unlist(rnb.get.annotation(region.type,assembly=assembly(rnb.set)))
  anno <- annotation(rnb.set)
  if(!is.null(annotation.filter)){
    anno <- anno[annotation.filter,]
  }
  if(!is.null(medecom.result@parameters$GROUP_LISTS)){
    sset <- medecom.result@parameters$GROUP_LIST[[cg_subset]]
    anno <- anno[sset,]
  }
  anno <- GRanges(Rle(anno$Chromosome),IRanges(start=anno$Start,end=anno$End))
  op <- findOverlaps(anno,agg.type)
  agg.type <- agg.type[subjectHits(op)]
  db.dir <- temp.dir
  dir <- downloadLolaDbs(tempdir(),"LOLACore")
  lola.db <- loadRegionDB(dir$hg38)
  lmcs <- getLMCs(medecom.result,cg_subset=cg_subset,K=K,lambda=lambda)
  lola.results <- list()
  for(i in 1:K){
    first.lmc <- lmcs[,i]
    ref.lmc <- rowMedians(lmcs[,-i])
    lmc.diff <- ref.lmc - first.lmc
    is.hypo <- lmc.diff > diff.threshold
    first.sites <- anno[is.hypo]
    op <- findOverlaps(first.sites,agg.type)
    first.sites <- agg.type[subjectHits(op)]
    lola.res <- NULL
    if(length(first.sites)>0){
      lola.res <- runLOLA(userSets=unique(first.sites),userUniverse=agg.type,regionDB=lola.db)
    }
    comp <- paste0("LMC",i)
    lola.results[[comp]] <- lola.res
  }
  return(lola.results)
}