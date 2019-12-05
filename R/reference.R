########################################################################################################################################
#' reference.R
#' -------------------------------------------------------------------------------------------------------------------------------------
#' This scripts contains functions for assocating MeDeCom with known reference cell type profiles
######################################################################################################################################## 

######################################################################################################################################## 
#' GLOBALS
AVAIL.REFS <- c("reinius","local")

######################################################################################################################################## 
#' FUNCTIONS

#' run.refbased
#' 
#' Function to link MeDeCom's output to reference methylation profiles. The function returns both the MeDeCom result and a methyalation
#' matrix with reference profiles from the source specified.
#' 
#' @param rnb.set An Object of type \code{RnBSet} from the \pkg{RnBeads} package containing methylation information on which
#'                 the deconvolution is to be performed
#' @param Ks Numeric vector containing the number of components to be computed by \pkg{MeDeCom}
#' @param lambdas Numeric vector specifying the regularzation parameters to be explored
#' @param cg_subsets List of numeric vectors specifying the rows that are to be used in the analysis.
#' @param opt.method Optimization method to be employed. For further information see \code{\link{runMeDeCom}}
#' @param temp.dir Optional temporary directory to store intermediate results
#' @param ref.base Reference profile data base to be used. Supported are 
#'                 \itemize{
#'                   \item \code{"reinius"} A blood cell type reference methylome data set from Reinius et.al. (Reference to be added)
#'                 }
#'                 \itemize{
#'                   \item \code{"local"} The reference data set is provided by the user. If this option is selected, \code{ref.set} 
#'                                      must not be empty.
#'                 }
#' @param most.var Number specifying the number of most variable to be selected from \code{rnb.set}
#' @param NCORES Number of cores to be used for analysis.
#' @param cluster.settings Setting for the environment of a high performance compute cluster. Passed to \code{\link{runMeDeCom}}
#' @param ref.set A \code{RnBSet} object containing a reference data set to be used besides \code{rnb.set}. Is only compatible with
#'                 \code{ref.base="local"}.
#' @param id.col The name of the column in the sample annotation sheet of \code{ref.set} containing the reference cell type. Is only
#'                 compatible with \code{ref.base="local"}.
#' @param save.restricted.sites Flag indicating if \code{rnb.set} restricted to the \code{most.var} sites is to be saved in \code{temp.dir}
#'                             for potential downstream analysis.                 
#' 
#' @return A list object containing two elements \itemize{
#'            \item \code{"MeDeComSet"} Results of applying MeDeCom with the setting above
#'            \item \code{"RefMeth"} A matrix containing reference profiles from the specified data set. The number of rows in 
#'                        this matrix has been reduced according to the most variable sites in \code{rnb.set}.
#'         }
#' 
#' @details This function applied MeDeCom to the specified data set and only support \code{RnBSet} objects as inputs. The function
#'  internally manipulated the object by selecting the most variable sites according to \code{most.var}. This leads to a decrease in
#'  the number of rows in the reference profiles to this number. 
#'  
#'  Please note that an active internet connection is required, since this routine downloads data through the world wide web.
#'  
#' @author Michael Scherer
#' 
#' @export
run.refbased <- function(rnb.set, 
                         Ks, 
                         lambdas,
                         cg_subsets=NULL,
                         opt.method="MeDeCom.cppTAfact",
                         temp.dir=NULL,
                         ref.base="reinius",
                         most.var=NULL,
                         NCORES=1,
                         cluster.settings=NULL,
                         ref.set=NULL,
                         id.col=NULL,
                         save.restricted.sites=FALSE){
  #require("RnBeads")
  res <- load.ref.set(ref.base,temp.dir)
  if(!is.null(res)){
    if(!is.null(ref.set)){
      warning("ref.set specified although ref.base is local. It will be ignored.")
    }
    if(!is.null(id.col)){
      warning("id.col specified although ref.base is local. It will be ignored.")
    }
    id.col <- res$ID
    ref.set <- res$rnb.set
  }else{
    if(is.null(ref.set)||is.null(id.col)){
      stop("Please specify ref.set and ref.id.col if.")
    }
  }

  rnb.set <- rnb.execute.imputation(rnb.set)
  anno.set <- annotation(rnb.set)
  anno.ref <- annotation(ref.set)
  anno.set <- GRanges(Rle(anno.set$Chromosome),IRanges(start=anno.set$Start,end=anno.set$End))
  anno.ref <- GRanges(Rle(anno.ref$Chromosome),IRanges(start=anno.ref$Start,end=anno.ref$End))
  op <- findOverlaps(anno.set,anno.ref)
  #anno.set <- anno.set[queryHits(op)]
  rem.sites <- !(1:nsites(rnb.set) %in% queryHits(op))
  rnb.set <- remove.sites(rnb.set,rem.sites)
  meth.data <- meth(rnb.set)
  if(!is.null(cg_subsets)){
    new.anno <- anno.set[!rem.sites]
    op2 <- findOverlaps(new.anno,anno.set)
    cg_subsets <- lapply(cg_subsets,function(x){
      matchi <- match(x,subjectHits(op2))
      queryHits(op2)[matchi[!is.na(matchi)]]
    })
  }
  if(!is.null(most.var)){
    vars <- apply(meth.data,1,function(x)var(x,na.rm=T))
    ordered <- order(vars,decreasing=T)
    rem.sites <- rep(TRUE,nsites(rnb.set))
    if(most.var>nsites(rnb.set)){
      warning("most.var bigger than number of sites, reduced to maximum")
      most.var <- nsites(rnb.set)
    }
    rem.sites[queryHits(op)][ordered[1:most.var]] <- FALSE
    #rem.sites[sample(1:nsites(rnb.set),most.var)] <- FALSE
    rnb.set <- remove.sites(rnb.set,rem.sites)
    if(save.restricted.sites){
      save.rnb.set(rnb.set,file.path(temp.dir,"restrictedRnBSet"))
    }
    anno.set <- anno.set[!rem.sites]
    op <- findOverlaps(anno.set,anno.ref)
  }
  meth.ref <- meth(ref.set)[subjectHits(op),]
  meth.ref <- rnb.execute.imputation(meth.ref)
  colnames(meth.ref) <- pheno(ref.set)[,id.col]
  medecom.result <- runMeDeCom(rnb.set,
                               Ks=Ks,
                               lambdas=lambdas,
                               cg_subsets=cg_subsets,
                               opt.method=opt.method,
                               temp.dir=temp.dir,
                               NCORES=NCORES,
                               cluster.settings=cluster.settings)
  if(!is.null(cg_subsets)){
    medecom.result@parameters$GROUP_LISTS <- cg_subsets
  }
  return(list("MeDeComSet"=medecom.result,"RefMeth"=meth.ref))
}

#' cluster.refbased
#' 
#' This routine performs general hierachical clustering as in the \code{\link{plotLMCs}} and adds reference methylome profiles to this
#' clustering. 
#' 
#' @param ref.run A result of \code{\link{run.refbased}} with the results of \code{\link{runMeDeCom}} and the reference data matrix
#' @param K Selected number of components
#' @param lambda Selected regularization parameter
#' @param cg_subset Numeric vector representing the CpG sites selected for analysis from the orignial MeDeComSet
#' @param plot.type Type of plot to be created, is passed to \code{\link{plotLMCs}}
#' @param ... Further arguments passed to plotLMCs
#' 
#' @return A plot object displaying the clustering
#' 
#' @author Michael Scherer
#' 
#' @export
cluster.refbased <- function(ref.run,
                             K,
                             lambda,
                             cg_subset=1,
                             plot.type="dendrogram",
                             ...){
  if(!is.list(ref.run) || length(ref.run)<2){
    stop("Argument needs to be the results obtained from 'run.refbased'")
  }
  medecom.result <- ref.run$MeDeComSet
  meth.ref <- ref.run$RefMeth
  if(!is.null(medecom.result@parameters$GROUP_LISTS)){
    sset <- medecom.result@parameters$GROUP_LISTS[[cg_subset]]
    meth.ref <- meth.ref[sset,]
  }
  plot <- plotLMCs(medecom.result,K=K,lambda=lambda,cg_subset=cg_subset,type=plot.type,Tref=meth.ref,...)
  return(plot)
}

#' plotLMC.reference
#' 
#' This routine uses a reference methylome, together with an annotation, and calls plotLMCs for the matching postions
#' 
#' @param medecom.set An object of type MeDeComSet
#' @param ref.meth Reference methylome in the form of a matrix of a RnBSet
#' @param ann.md Genomic annotation for the CpG sites present in medecom.set
#' @param ann.ref Genomic annotation for the CpG sites in the referenc methylome. Can be omitted, if ref.meth is an RnBSet object
#' @param cg_subset The cg_subset of interest
#' @param K K value
#' @param lambda lambda value
#' @param type Plot type, see \code{\link{plotLMCs}}
#' @param chrom.col The chromosome column name in ann.md and ann.ref
#' @param start.col The start column name in ann.md and ann.ref
#' @param end.col The end column name in ann.md and ann.ref
#' @param ct.color.column A column name in the phenotypic information of ref.meth to be shown in the clustering
#' @return NULL
#' @author Michael Scherer
#' @export

plotLMC.reference <- function(medecom.set,
                              ref.meth,
                              ann.md,
                              cg_subset=medecom.set@parameters$cg_subsets[1],
                              K=medecom.set@parameters$Ks[1],
                              lambda=medecom.set@parameters$lambdas[1],
                              ann.ref=NULL,
                              type="dendrogram",
                              chrom.col="Chromosome",
                              start.col="Start",
                              end.col="End",
                              ct.color.column="cell_type"){
  #require("RnBeads")
  if(!inherits(medecom.set,"MeDeComSet")){
    stop("Invalid value for medecom.set")
  }
  if(!inherits(ref.meth,"RnBSet") && !is.matrix(ref.meth)){
    stop("Invalid value for ref.meth; needs to be RnBSet or matrix")
  }
  if(is.matrix(ref.meth) && is.null(ann.ref)){
    stop("Genomic annotation needs to be provided for the reference methylome")
  }
  if(inherits(ref.meth,"RnBSet")){
    ann.ref <- annotation(ref.meth)
    meth.ref <- meth(ref.meth)
    colnames(meth.ref) <- pheno(ref.meth)[,ct.color.column]
    ref.meth <- meth.ref
    rm(meth.ref)
  }
  ann.md <- GRanges(Rle(ann.md[,chrom.col]),IRanges(start=ann.md[,start.col],end=ann.md[,end.col]))
  ann.ref <- GRanges(Rle(ann.ref[,chrom.col]),IRanges(start=ann.ref[,start.col],end=ann.ref[,end.col]))
  op <- findOverlaps(ann.md,ann.ref)
  if(any(!(1:length(ann.md) %in% queryHits(op)))){
    stop("Insuffiecient information present in reference methylome")
  }
  ref.meth <- ref.meth[subjectHits(op),]
  ref.meth <- ref.meth[medecom.set@parameters$GROUP_LISTS[[cg_subset]],]
  ref.meth <- rnb.execute.imputation(ref.meth,method="knn")
  plotLMCs(medecom.set,Tref=ref.meth,K=K,lambda=lambda,cg_subset=cg_subset,type = type)
}
#' load.ref.set
#' 
#' This functions loads a reference data base and return the corresponding \code{RnBSet}.
#' 
#' @param ref.base Reference base to be used. See \code{\link{run.refbased}} for further information.
#' @param temp.dir Temporary directory to store the object.
#' 
#' @return List of two elements \itemize{
#'           \item ID sample identifier column of the refernce data set used for plotting
#'           \item rnb.set \code{RnBSet} object containing the reference methylation profiles
#' }
#' 
#' @author Michael Scherer
#' 
#' @noRd

load.ref.set <- function(ref.base,temp.dir=NULL){
  #require("RnBeads")
  if(!ref.base %in% AVAIL.REFS){
    stop(paste("Unsupported reference data base, must be one of",AVAIL.REFS))
  }
  if(is.null(temp.dir)){
    temp.dir <- tempdir()
  }
  if(ref.base=="local") return(NULL)
  id.col <- "sample_id"
  if(ref.base=="reinius"){
    location <- file.path(temp.dir,"Reinius_Blood_Reference.zip")
    if(!file.exists(location)){
      cat("Downloading Reinius reference set \n")
      downloaded <- tryCatch(download.file("http://rnbeads.mpi-inf.mpg.de/publication/Reinius_Blood_Reference.zip",destfile = location),error=function(e){
        if(inherits(e,"Error")){
          stop("Failed to download reference data set. Check internet connection.")
        }
      })
    }
    ref.set <- load.rnb.set(location)
    id.col <- "Cell/Tissue"
  }
  return(list(ID=id.col,rnb.set=ref.set))
}
