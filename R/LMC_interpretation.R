###############################################################################################################
# LMC_interpreation.R
# functions to interpret LMCs and to assign potential functions to the LMCs
###############################################################################################################

#' lmc.lola.enrichment
#' 
#' This routine computes LOLA enrichment results for LMC-specifically hypo- or hypermethylated sites. 
#' 
#' @param medecom.result An object of type \code{\link{MeDeComSet-class}} or the location of an .RData file, where
#'                 such an object is stored.
#' @param annotation.filter A numeric vector specifying the sites that have been removed from \code{rnb.set} in a
#'                 preprocessing step (e.g. coverage filtering) or a path to an .RData file.
#' @param rnb.set The original \code{\link[RnBeads]{RnBSet-class}} object containing methylation, sample meta and annotation
#'                 information, a path to a directory stored by \code{\link[RnBeads]{save.rnb.set}} or a data.frame containing
#'                 CpG annotations (ann_C)
#' @param K The number of LMCs specified for the MeDeCom run.
#' @param lambda The lambda parameter selected.
#' @param cg_subset The index of the selection strategy employed (e.g. most variable CpGs).
#' @param diff.threshold The difference cutoff between median methylation in the remaining LMCs and the LMC of interest
#'                  used to call a CpG differentially methylated. The higher this value, the more conservative the
#'                  selection.
#' @param region.type Region type used to annotate CpGs to potentially regulatory regions (see \url{https://rnbeads.org/regions.html})
#'                  for a list of available region types.
#' @param temp.dir Path to a directory used to store temporary files.
#' @param type Which direction is to be tested for enrichment. Can be one of "hypo", "hyper", or "differential"
#' @param assembly The assembly used. Needs to be one of "hg19", "hg38" or "mm10". Does not need to be specified, if rnb.set is a
#'                 \code{\link{RnBSet-class}}
#' @param lola.db A loaded LOLA database as loaded with \code{\link{loadRegionDB}}. If this value is NULL, the database is loaded
#'                 automatically and stored in the temporary directory.
#' @return A list with K elements. One element is the enrichment result of the corresponding LMC-specific hypomethylated CpG sites.
#' @export
#' @details This function employs LOLA on the CpG sites that are LMC-specifically hypomethylated, after annotating 
#'                  the sites to the closest region defined by \code{region.type}. The sites are selected by computing
#'                  the median methylation value of the other LMCs and then selecting those sites that are more than 
#'                  \code{diff.threshold} away from the median in the LMC of interest. This is done for all LMCs from
#'                  1 to K.
#' @seealso lmc.lola.plot.tables
#' @author Michael Scherer                                

lmc.lola.enrichment <- function(medecom.result,
                                annotation.filter=NULL,
                                rnb.set,
                                K=NULL,
                                lambda=NULL,
                                cg_subset=NULL,
                                diff.threshold=0.5,
                                region.type="ensembleRegBuildBPall",
                                temp.dir=tempdir(),
                                type="hypo",
                                assembly="hg19",
                                lola.db=NULL){
  require("RnBeads")
  require("LOLA")
  rnb.mode <- F
  if(is.null(cg_subset)){
    cg_subset <- medecom.result@parameters$cg_subsets[1]
  }else if(!(cg_subset %in% medecom.result@parameters$cg_subsets)){
    stop("Invalid value for cg_subset; not available in medecom.result")
  }
  if(is.null(K)){
    K <- medecom.result@parameters$Ks[1]
  }else if(!(K %in% medecom.result@parameters$Ks)){
    stop("Invalid value for K; not available in medecom.result")
  }
  if(is.null(lambda)){
    lambda <- medecom.result@parameters$lambdas[1]
  }else if(!(lambda %in% medecom.result@parameters$lambdas)){
    stop("Invalid value for lambdas; not available in medecom.result")
  }
  if(inherits(rnb.set,"RnBSet")){
    rnb.mode <- T
    assembly <- assembly(rnb.set)
  }
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
    rnb.load.annotation.from.db(region.type,assembly=assembly)
  }
  agg.type <- unlist(rnb.get.annotation(region.type,assembly=assembly))
  if(rnb.mode){
    anno <- annotation(rnb.set)
  }else{
    anno <- rnb.set
    rm(rnb.set)
  }
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
  if(is.null(lola.db)){
    lola.db <- load.lola.for.medecom(temp.dir)
  }
  lmcs <- getLMCs(medecom.result,cg_subset=cg_subset,K=K,lambda=lambda)
  lola.results <- list()
  for(i in 1:K){
    first.lmc <- lmcs[,i]
    ref.lmc <- rowMedians(lmcs[,-i,drop=F])
    lmc.diff <- ref.lmc - first.lmc
    if(type=="hypo"){
      is.hypo <- lmc.diff > diff.threshold
    }else if (type=="hyper"){
      is.hypo <- lmc.diff < -diff.threshold
    }else if (type=="differential"){
      is.hypo <- (lmc.diff < -diff.threshold) | (lmc.diff > diff.threshold)
    }else{
      stop("Invalid value for type, needs to be one of 'hypo', 'hyper' or 'differential'")
    }
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

#' load.lola.for.medecom
#' 
#' This functions downloads and loads the LOLA database in the specified directory. Should only be called once per session to save time.
#' 
#' @param dir.path A path to a directory, where the LOLA database is to be downloaded. Defaults to the temporary directory.
#' @return The loaded LOLA database
#' @export
#' @author Michael Scherer
load.lola.for.medecom <- function(dir.path=tempdir()){
  require("LOLA")
  dir <- downloadLolaDbs(dir.path,"LOLACore")
  lola.db <- loadRegionDB(dir$hg38)
  return(lola.db)
}

#' lmc.lola.plots.tables
#' 
#' This functions calls \link{lmc.lola.enrichment} and returns plots representing those results, as well as the tables with LOLA
#' enrichment results.
#' 
#' @param medecom.result An object of type \code{\link{MeDeComSet-class}} or the location of an .RData file, where
#'                 such an object is stored.
#' @param annotation.filter A numeric vector specifying the sites that have been removed from \code{rnb.set} in a
#'                 preprocessing step (e.g. coverage filtering) or a path to an .RData file.
#' @param rnb.set The original \code{\link[RnBeads]{RnBSet-class}} object containing methylation, sample meta and annotation
#'                 information, a path to a directory stored by \code{\link[RnBeads]{save.rnb.set}} or a data.frame containing
#'                 CpG annotations (ann_C)
#' @param K The number of LMCs specified for the MeDeCom run.
#' @param lambda The lambda parameter selected.
#' @param cg_subset The index of the selection strategy employed (e.g. most variable CpGs).
#' @param diff.threshold The difference cutoff between median methylation in the remaining LMCs and the LMC of interest
#'                  used to call a CpG differentially methylated. The higher this value, the more conservative the
#'                  selection.
#' @param region.type Region type used to annotate CpGs to potentially regulatory regions (see \url{https://rnbeads.org/regions.html})
#'                  for a list of available region types.
#' @param temp.dir Path to a directory used to store temporary files.
#' @param type Which direction is to be tested for enrichment. Can be one of "hypo", "hyper", or "differential"
#' @param assembly The assembly used. Needs to be one of "hg19", "hg38" or "mm10". Does not need to be specified, if rnb.set is a
#'                 \code{\link{RnBSet-class}}
#' @param lola.db A loaded LOLA database as loaded with \code{\link{loadRegionDB}}. If this value is NULL, the database is loaded
#'                 automatically and stored in the temporary directory.
#' @return A list with two elements, one of them containing the plots for each LMC and the other for the corresponding LOLA
#'         enrichment tables
#' @seealso lmc.lola.enrichment
#' @author Michael Scherer                 
#' 
#' @export

lmc.lola.plots.tables <- function(medecom.result,
                           annotation.filter=NULL,
                           rnb.set,
                           K=NULL,
                           lambda=NULL,
                           cg_subset=NULL,
                           diff.threshold=0.5,
                           region.type="ensembleRegBuildBPall",
                           temp.dir=tempdir(),
                           type="hypo",
                           assembly="hg19",
                           lola.db=NULL){
  enrichment.results <- lmc.lola.enrichment(medecom.result,
                                          annotation.filter=NULL,
                                          rnb.set,
                                          K=NULL,
                                          lambda=NULL,
                                          cg_subset=NULL,
                                          diff.threshold=0.5,
                                          region.type="ensembleRegBuildBPall",
                                          temp.dir=tempdir(),
                                          type="hypo",
                                          assembly="hg19",
                                          lola.db=NULL)
  lola.plots <- lapply(enrichment.results,do.lola.plot,lola.db)
  return(list(Plots=lola.plots,Tables=enrichment.results))
}

#' do.lola.plot
#' 
#' This functions creates a single LOLA enrichment plot
#' 
#' @param enrichment.result LOLA enrichment result for a single LMC
#' @param lola.db LOLA database that was used
#' @param pvalCut P-value cutoff
#' @return An object of type ggplot, containing the enrichment plot
#' @author Michael Scherer
#' @noRd

do.lola.plot <- function(enrichment.result,lola.db,pvalCut=0.01){
  plot <- lolaBarPlot(lolaDb = lola.db,lolaRes=enrichment.result,pvalCut=pvalCut)+theme_bw()+theme(axis.text.x = element_text(angle=45,hjust = 1))
  return(plot)
}

#' lmc.go.enrichment
#' 
#' This routine computes GO enrichment results for LMC-specifically hypo- or hypermethylated sites. 
#' 
#' @param medecom.result An object of type \code{\link{MeDeComSet-class}} or the location of an .RData file, where
#'                 such an object is stored.
#' @param annotation.filter A numeric vector specifying the sites that have been removed from \code{rnb.set} in a
#'                 preprocessing step (e.g. coverage filtering) or a path to an .RData file.
#' @param rnb.set The original \code{\link[RnBeads]{RnBSet-class}} object containing methylation, sample meta and annotation
#'                 information or a path to a directory stored by \code{\link[RnBeads]{save.rnb.set}} or a data.frame containing
#'                 CpG annotations (ann_C)
#' @param K The number of LMCs specified for the MeDeCom run.
#' @param lambda The lambda parameter selected.
#' @param cg_subset The index of the selection strategy employed (e.g. most variable CpGs).
#' @param diff.threshold The difference cutoff between median methylation in the remaining LMCs and the LMC of interest
#'                  used to call a CpG differentially methylated. The higher this value, the more conservative the
#'                  selection.
#' @param region.type Region type used to annotate CpGs to potentially regulatory regions (see \url{https://rnbeads.org/regions.html})
#'                  for a list of available region types. Here, only "genes" "promoters" and their gencode versions
#'                  are available.
#' @param temp.dir Path to a directory used to store temporary files.
#' @param type Which direction is to be tested for enrichment. Can be one of "hypo", "hyper", or "differential"
#' @param assembly The assembly used. Needs to be one of "hg19", "hg38" or "mm10". Does not need to be specified, if rnb.set is a
#'                 \code{\link{RnBSet-class}}
#' @return A list with K elements. One element is the enrichment result of the corresponding LMC-specific hypomethylated CpG sites.
#' @export
#' @details This function employs GO enrichment analysis with the GOstats package on the CpG sites that are LMC-specifically
#'                  hypomethylated, after annotating the sites to the closest promotor/gene defined by \code{region.type}.
#'                  The sites are selected by computing the median methylation value of the other LMCs and then selecting
#'                  those sites that are more than \code{diff.threshold} away from the median in the LMC of interest.
#'                  This is done for all LMCs from 1 to K.
#' @author Michael Scherer                                

lmc.go.enrichment <- function(medecom.result,
                                  annotation.filter=NULL,
                                  rnb.set,
                                  K=NULL,
                                  lambda=NULL,
                                  cg_subset=NULL,
                                  diff.threshold=0.5,
                                  region.type="genes",
                                  temp.dir=tempdir(),
                                  type="hypo",
                                  assembly="hg19"){
  require("RnBeads")
  require("GOstats")
  rnb.mode <- F
  if(is.null(cg_subset)){
    cg_subset <- medecom.result@parameters$cg_subsets[1]
  }else if(!(cg_subset %in% medecom.result@parameters$cg_subsets)){
    stop("Invalid value for cg_subset; not available in medecom.result")
  }
  if(is.null(K)){
    K <- medecom.result@parameters$Ks[1]
  }else if(!(K %in% medecom.result@parameters$Ks)){
    stop("Invalid value for K; not available in medecom.result")
  }
  if(is.null(lambda)){
    lambda <- medecom.result@parameters$lambdas[1]
  }else if(!(lambda %in% medecom.result@parameters$lambdas)){
    stop("Invalid value for lambdas; not available in medecom.result")
  }
  if(inherits(rnb.set,"RnBSet")){
    rnb.mode <- T
    assembly <- assembly(rnb.set)
  }
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
    rnb.load.annotation.from.db(region.type,assembly=assembly)
  }
  agg.type <- unlist(rnb.get.annotation(region.type,assembly=assembly))
  entrez.id <- values(agg.type)$entrezID
  longer <- lengths(strsplit(entrez.id,";"))>0
  entrez.id[longer] <- unlist(lapply(strsplit(entrez.id,";"),function(x)x[1]))[longer]
  if(rnb.mode){
    anno <- annotation(rnb.set)
  }else{
    anno <- rnb.set
  }
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
  entrez.id <- entrez.id[subjectHits(op)]
  lmcs <- getLMCs(medecom.result,cg_subset=cg_subset,K=K,lambda=lambda)
  go.results <- list()
  for(i in 1:K){
    first.lmc <- lmcs[,i]
    ref.lmc <- rowMedians(lmcs[,-i,drop=F])
    lmc.diff <- ref.lmc - first.lmc
    if(type=="hypo"){
      is.hypo <- lmc.diff > diff.threshold
    }else if (type=="hyper"){
      is.hypo <- lmc.diff < -diff.threshold
    }else if (type=="differential"){
      is.hypo <- (lmc.diff < -diff.threshold) | (lmc.diff > diff.threshold)
    }else{
      stop("Invalid value for type, needs to be one of 'hypo', 'hyper' or 'differential'")
    }
    first.sites <- anno[is.hypo]
    op <- findOverlaps(first.sites,agg.type)
    first.ids <- entrez.id[subjectHits(op)]
    go.res <- NULL
    if(length(first.ids)>0){
      params <- new("GOHyperGParams",annotation="org.Hs.eg.db",geneIds = first.ids, universeGeneIds = entrez.id, ontology = "BP",conditional = TRUE, testDirection = "over")
      test.result.hyper <- tryCatch(
        hyperGTest(params),
        error = function(e) {
          print(e)
        }
      )
      if(!inherits(test.result.hyper,"error")){
        go.res <- summary(test.result.hyper)
        if(!is.null(dim(go.res)) && !is.null(go.res) && nrow(go.res)>0){
            go.res$p.val.adj.fdr <- p.adjust(go.res$Pvalue,method="fdr",n=length(test.result.hyper@pvalue.order))
        }
      }else{
        go.res <- NA
      }
    }
    comp <- paste0("LMC",i)
    go.results[[comp]] <- go.res
  }
  return(go.results)
}
