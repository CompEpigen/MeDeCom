#' MeDeCom: Methylome DeComposition using regularized constrained matrix factorization
#'
#' MeDeCom is an R-package that discovers and quantifies latent components
#' in the DNA methylomes of heterogeneous samples
#' 
#' @references TBA
#' @docType package  
#' @name MeDeCom-package
#' @useDynLib MeDeCom
NULL

#' Example CG annotation
#' 
#' An data frame containing CpGs with their corresponding annotation in the genome to be used for the analysis. You can provide a similar
#' \code{data.frame} for your own analysis.
#' 
#' @docType data
#' @keywords datasets
#' @name example.cg.annotation
#' @format \code{cg.ann} is a \code{data.frame} containing identifiers and postitions of CpG sites to be analyzed
#' @author Michael Scherer
NULL

#' Example dataset
#' 
#' Contains a data set of methylation values as a \code{matrix} with 10,000 rows (CpGs) and 100 columns (samples), a reference of
#' contributions of 5 cell types in the samples, and the reference methylomes of the samples as a matrix with 10,000 rows and 5
#' columns.
#'  
#' \itemize{
#'           \item D, a \code{matrix} with 10,000 rows (CpGs) and 100 columns (samples) representing a potential input methylation
#'           matrix
#'           \item Aref, a reference of contributions of 5 cell types in the samples
#'           \item Tref, reference methylomes of the samples as a matrix with 10,000 rows and 5 columns
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name example.dataset
#' @usage data(example.dataset)
#' @author Michael Scherer
NULL

#' Example result object (MeDeComSet)
#' 
#' An example result of applying \pkg{MeDeCom} for a deconvolution experiment.
#'  
#' @docType data
#' @keywords datasets
#' @name example.MeDeComSet
#' @format \code{MeDeComSet}, the result of a deconvolution experiment and input to any follow-up analysis or to retrieve latent 
#' methylation components (LMCs)
#' @author Michael Scherer
NULL