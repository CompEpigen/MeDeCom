#' MeDeComSet Class
#'
#' Stores the result of a methylation deconvolution experiment
#'
#' @section Slots:
#' \describe{
#'   \item{\code{dataset_info}}{\code{list} with information about the input data set.}
#'   \item{\code{parameters}}{\code{list} containing parameters of the deconvolution experiment.}
#'   \item{\code{outputs}}{\code{list} of deconvolution products for each combination of MeDeCom parameters.}
#' }
#' @section Methods and functions:
#' \describe{
#'   \item{\code{\link[=getStatistics,MeDeComSet-method]{getStatistics}}}{Returns the value of one goodness-of-fit statistics for a given parameter combination.}
#'   \item{\code{\link[=getLMCs,MeDeComSet-method]{getLMCs}}}{Returns a matrix of LMCs for a given parameter combination.}
#'   \item{\code{\link[=getProportions,MeDeComSet-method]{getProportions}}}{Returns a matrix of mixing proportions for a given parameter combination.}
#'   \item{\code{\link{plotParameters}}}{Create a parameter selection plot.}
#'   \item{\code{\link{plotLMCs}}}{Visualize the LMCs.}
#'   \item{\code{\link{plotProportions}}}{Visualize the mixing proportions.}
#' }
#' 
#'
#' @name MeDeComSet-class
#' @rdname MeDeComSet-class
#' @author Pavlo Lutsik
#' @exportClass MeDeComSet
setClass("MeDeComSet",
		representation(
				dataset_info="list",
				parameters="list",
				outputs="list"
		),
		prototype(
				dataset_info=list(),
				parameters=list(),
				outputs=list()				
		),
		package = "MeDeCom")
########################################################################################################################
#setMethod("initialize", "MeDeComSet",
#		function(.Object,
#				outputs=list(),
#				parameters=list()
#		) {
#			
#			.Object@outputs<-outputs
#			.Object@parameters<-parameters
#			.Object
#})
########################################################################################################################
#' MeDeComSet
#' 
#' Wrapper function MeDeComSet
#'
#'
#' @param parameters       	\code{list} of MeDeCom parameters with elements \code{K}, integer vector of k values, \code{lambdas}, numeric vector of lambda values.
#' @param outputs			\code{list} of MeDeCom resutls with one element per each used CpG subset.
#' @param dataset_info   \code{list} with information about the input data set.
#' 
#' @return an object of class MeDeComSet
#' 
#' @name MeDeComSet
#' @rdname MeDeComSet-class
#' @aliases initialize,MeDeComSet-method
#' @export
MeDeComSet<-function(
		parameters,
		outputs,
		dataset_info=list()){
	
	object<-new("MeDeComSet",
			dataset_info = dataset_info,
			parameters = parameters,
			outputs = outputs)
	object
}
########################################################################################################################
if(!isGeneric("getStatistics")) setGeneric("getStatistics",
			function(object, ...) standardGeneric("getStatistics"))

#' getStatistics-methods
#'
#' Methylation sites object information for which is present in the \code{RnBSet} object.
#'
#' @param object   object returned by \link{runMeDeCom}
#' @param Ks			   numbers of LMCs
#' @param lambdas	   regularlization parameters
#' @param cg_subset	   used CpG subset, defaults to the full data set
#' @param statistic   \code{character} of length 1 specifying goodness of fit statistics
#' 
#' @details 
#' Currently the following values for \code{statistics} can be supplied: \code{objective}, \code{RMSE}, \code{CVE}.
#'
#' @return A numeric \code{matrix} or \code{vector} with the requested statistics
#'
#' @rdname getStatistics-methods
#' @docType methods
#' @aliases getStatistics
#' @aliases getStatistics,MeDeComSet-method
#' @export
#' @examples
#' \donttest{
#' data(example.data)
#' getStatistics(example_MeDeComSet, K=2, lambda=0.001)
#' }
setMethod("getStatistics", signature(object="MeDeComSet"),
		function(object, Ks=object@parameters$Ks, lambdas=object@parameters$lambdas, cg_subset=1, statistic="cve"){
		check_inputs(object, cg_subset, Ks, lambda=lambdas)
		elt<-c(
				"objective"="Fval",	"Fval"="Fval", 
				"rmse"="rmse", "RMSE"="rmse",
				"CVE"="cve", "cve"="cve",
				"MAEA"="maeA", "maeA"="maeA",
				"RMSET"="rmseT","rmseT"="rmseT",
				"deviance"="Deviance", "Deviance"="Deviance",
				"rss"="RSS","RSS"="RSS"
		)[statistic]
		return(as.numeric(object@outputs[[cg_subset]][[elt]][match(Ks, object@parameters$Ks), match(lambdas, object@parameters$lambdas)]))
})
########################################################################################################################
if(!isGeneric("getLMCs")) setGeneric("getLMCs",
			function(object, ...) standardGeneric("getLMCs"))
#'
#' getLMCs-methods
#' 
#' Return a matrix of LMCs
#' 
#' @param object   object returned by \link{runMeDeCom}
#' @param K			   number of LMCs
#' @param lambda	   regularlization parameter
#' @param cg_subset	   used CpG subset, defaults to the full data set
#' @param statistic    statistic to be used in returning
#' 
#' @rdname getLMCs-methods
#' @docType methods
#' @aliases getLMCs
#' @aliases getLMCs,MeDeComSet-method
#' @export
#' @examples
#' \donttest{
#' data(example.data)
#' getLMCs(example_MeDeComSet, K=2, lambda=0.001)
#' }
setMethod("getLMCs", signature(object="MeDeComSet"),
		function(object, K=object@parameters$Ks[1], lambda=object@parameters$lambdas[1], cg_subset=1, statistic="cve"){
	check_inputs(object, cg_subset, K, lambda)
	return(object@outputs[[cg_subset]]$T[[match(K, object@parameters$Ks), match(lambda, object@parameters$lambdas)]])
})
########################################################################################################################
if(!isGeneric("getProportions")) setGeneric("getProportions",
			function(object, ...) standardGeneric("getProportions"))
#'
#' getProportions-methods
#' 
#' Return a matrix of LMCs
#' 
#' @param object   object returned by \link{runMeDeCom}
#' @param K			   number of LMCs
#' @param lambda	   regularlization parameter
#' @param cg_subset	   used CpG subset, defaults to the full data set
#' @param statistic    statistic to be used in returning
#' 
#' @rdname getProportions-methods
#' @docType methods
#' @aliases getProportions
#' @aliases getProportions,MeDeComSet-method
#' @export
#' @examples
#' \donttest{
#' data(example.data)
#' getProportions(example_MeDeComSet, K=2, lambda=0.001)
#' }
setMethod("getProportions", signature(object="MeDeComSet"),
		function(object, K=object@parameters$Ks[1], lambda=object@parameters$lambdas[1], cg_subset=1, statistic="cve"){
	check_inputs(object, cg_subset, K, lambda)
	Ahat<-object@outputs[[cg_subset]]$A[[match(K, object@parameters$Ks), match(lambda, object@parameters$lambdas)]]
	if(!is.null(dim(Ahat))){
	  rownames(Ahat)<-sprintf("LMC%d", 1:nrow(Ahat))
  	if(!is.null(object@dataset_info$sample_names)){
	  	colnames(Ahat)<-object@dataset_info$sample_names
  	}
	}
	return(Ahat)
})
########################################################################################################################
#
# Check inputs for the get* methods
#
check_inputs<-function(MeDeComSet, cg_subset, K, lambda){
	if(!all(cg_subset %in% MeDeComSet@parameters$cg_subsets)){
		stop("wrong cg subset supplied")
	}
	
	if(!all(K %in% MeDeComSet@parameters$Ks)){
		stop("wrong K value supplied")
	}
	
	if(!all(lambda %in% MeDeComSet@parameters$lambdas)){
		stop("wrong lambda value supplied")
	}
}
#' as.MeDeComSet
#' 
#' Function to convert object of type RefFreeCellMix to MeDeComSet
#' 
#' @param object An object of class \code{RefFreeCellMix} containing cell type deconvolution information, or a list of such objects.
#' @param cg_subsets The indeces of the CpG subsets used in the analysis.
#' @param Ks The values of K used in the analysis. If NULL, K is determined by the size of the matrices.
#' @param deviances Optional argument specifying the deviances as computed with \code{\link{RefFreeCellMixArrayDevianceBoots}}.
#' @param m.orig The original number of rows (CpGs) in the methylation matrix.
#' @param n.orig The original number of columns (samples) in the methylation matrix.
#' @return An object of type \code{MeDeComSet}
#' @details Since \code{RefFreeCellMix} only contains information on a single value for K, and does not contain any regularization
#'           (lambda), the corresponding parameters in the MeDeComSet are set to single numeric values. Furthermore, no information
#'           on goodness of fit (CVE, Fval) can be stored. If \code{cg_subsets} is not of length 1, an object containing multiple
#'           subsets is creared. 
#' @export
as.MeDeComSet <- function(object,cg_subsets=1,Ks=NULL,deviances=NULL,rss=NULL,m.orig=NULL,n.orig=NULL){
  c.obj <- class(object)
  if(c.obj=="list"){
    c.obj <- class(object[[1]])
    if(c.obj == "list"){
      c.obj <- class(object[[1]][[1]])
    }
  }  
  if(!(c.obj=="RefFreeCellMix" | c.obj=="list")){
      stop(paste("Cannot convert object of type",c.obj,"to MeDeComSet"))
  }
  if(c.obj == "RefFreeCellMix"){
    output <- list()
    if(is.null(Ks)){
      Ks <- "1"
      all.Ks <- ncol(object$Omega)
      object <- list("1"=object)
      if(!is.null(deviances)){
        deviances <- list("1"=deviances)
      }
    }else{
      all.Ks <- Ks
      Ks <- as.character(Ks)
    }
    if(length(cg_subsets)==1){
      object <- list(object)
      if(!is.null(deviances)){
        deviances <- list(deviances)
      }
    }
    for(ssets in cg_subsets){
      sel.sset <- object[[ssets]]
      lambda <- 0
      T.all <- list()
      A.all <- list()
      for(i in 1:length(Ks)){
        K <- Ks[i]
        sel.object <- sel.sset[[K]]
        A <- sel.object$Omega
        T <- sel.object$Mu
        K <- ncol(A)
        T.all[[i]] <- T
        A.all[[i]] <- t(A)
      }
      T.all <- matrix(T.all,nrow=length(Ks))
      row.names(T.all) <- paste("K",Ks,sep="_")
      colnames(T.all) <- paste("lambda",lambda,sep = "_")
      A.all <- matrix(A.all,nrow = length(Ks))
      row.names(A.all) <- paste("K",Ks,sep="_")
      colnames(A.all) <- paste("lambda",lambda,sep = "_")
      if(is.null(deviances)){
        output[[ssets]] <- list(T=T.all,A=A.all)
      }else{
        deviances.all <- matrix(deviances[[ssets]],nrow = length(Ks))
        row.names(deviances.all) <- paste("K",Ks,sep="_")
        colnames(deviances.all) <- paste("lambda",lambda,sep = "_")
        output[[ssets]] <- list(T=T.all,A=A.all,Deviance=deviances.all)
      }
    }
    parameters <- list(cg_subsets=cg_subsets,
                       Ks=all.Ks,
                       lambdas=0)
    if(is.null(m.orig)){
      m <- nrow(T)
    }else{
      m <- m.orig
    }
    if(is.null(n.orig)){
      n <- ncol(A)
    }else{
      n <- n.orig
    }
    d.info <- list(m=m,n=n,TYPE="RefFreeCellMix")
    new.obj <- MeDeComSet(parameters = parameters,
                          outputs = output,
                          dataset_info = d.info)
  }else if(c.obj=="list"){
    output <- list()
    if(is.null(Ks)){
      Ks <- "1"
      all.Ks <- ncol(object$T$T.all$T.list[[1]])
      object <- list("1"=object)
      if(!is.null(rss)){
        rss <- list("1"=rss)
      }
    }else{
      all.Ks <- Ks
      Ks <- as.character(Ks)
    }
    if(length(cg_subsets)==1){
      object <- list(object)
      if(!is.null(rss)){
        rss <- list(rss)
      }
    }
    for(ssets in cg_subsets){
      sel.sset <- object[[ssets]]
      lambda <- 0
      T.all <- sel.sset$T
      A.all <- sel.sset$A
      T.all <- matrix(T.all,nrow=length(Ks))
      row.names(T.all) <- paste("K",Ks,sep="_")
      colnames(T.all) <- paste("lambda",lambda,sep = "_")
      A.all <- matrix(A.all,nrow = length(Ks))
      row.names(A.all) <- paste("K",Ks,sep="_")
      colnames(A.all) <- paste("lambda",lambda,sep = "_")
      if(is.null(rss)){
        output[[ssets]] <- list(T=T.all,A=A.all)
      }else{
        rss.all <- matrix(rss[[ssets]],nrow = length(Ks))
        row.names(rss.all) <- paste("K",Ks,sep="_")
        colnames(rss.all) <- paste("lambda",lambda,sep = "_")
        output[[ssets]] <- list(T=T.all,A=A.all,RSS=rss.all)
      }
    }
    parameters <- list(cg_subsets=cg_subsets,
                       Ks=all.Ks,
                       lambdas=0)
    if(is.null(m.orig)){
      m <- nrow(T)
    }else{
      m <- m.orig
    }
    if(is.null(n.orig)){
      n <- ncol(A)
    }else{
      n <- n.orig
    }
    d.info <- list(m=m,n=n,TYPE="EDec")
    new.obj <- MeDeComSet(parameters = parameters,
                          outputs = output,
                          dataset_info = d.info)
    
  }
  return(new.obj)
}
########################################################################################################################
setMethod("show", "MeDeComSet", function(object){
			cat("An object of class MeDeComSet\n")
			cat("Input data set:\n")
			cat(sprintf("\t%d CpGs\n", object@dataset_info$m))
			cat(sprintf("\t%d methylomes\n", object@dataset_info$n))
			cat("Experimental parameters:\n")
			#cat(sprintf("\tCpG subsets: %s\n", paste(object@parameters$lambdas, collapse=", ")))
			cat(sprintf("\tk values: %s\n", paste(object@parameters$Ks, collapse=", ")))
			cat(sprintf("\tlambda values: %s\n", paste(object@parameters$lambdas, collapse=", ")))
			
		})
########################################################################################################################