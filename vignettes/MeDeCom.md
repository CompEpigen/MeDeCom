---
title: "MeDeCom: Methylome Decomposition via Constrained Matrix Factorization"
author: "Pavlo Lutsik, Martin Slawski, Gilles Gasparoni, Nikita Vedeneev, Matthias Hein and Joern Walter"
date: "2020-03-12"
output:
  rmarkdown::html_document:
    mathjax: default
    toc: true
    number_sections: false
    fig_width: 5
    fig_height: 5
vignette: >
  %\VignetteIndexEntry{MeDeCom}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Introduction

*MeDeCom* is an R-package for reference-free decomposition of heterogeneous DNA methylation profiles. 
It uses matrix factorization enhanced by constraints and a specially tailored regularization. 
*MeDeCom* represents an input $m\times n$ data matrix ($m$ CpGs measured in $n$ samples) as a product of two other matrices. 
The first matrix has $m$ rows, just as the input data, but the number of columns is equal to $k$. 
The columns of this matrix can be interpreted as methylomes of the $k$ unknown 
cell populations underlying the samples and will be referred to as **latent methylation components** or **LMCs**. 
The second matrix has $k$ rows and $n$ columns, and can be interpreted as a matrix of relative contributions (mixing proportions) 
of each LMC to each sample.

*MeDeCom* starts with a set of related DNA methylation profiles, e.g. a series of Infinium microarray measurements
from a population cohort, or several bisulfite sequencing-based methylomes. The key requirement is that
 the input data represents absolute DNA methylation measurements in a population of cells 
and contains values between 0 and 1. *MeDeCom* implements an alternating scheme which iteratively 
updates randomly initialized factor matrices until convergence or until the maximum number of iterations 
has been reached. This is repeated for multiple random initializations and the best solution is returned.

MeDeCom features two tunable parameters. The first one is the number of LMCs $k$, an approximate choice for which 
should be known from prior information. To enforce the distribution properties of a methylation profile 
upon LMCs *MeDeCom* uses a special for of regularization controlled by the parameter $\lambda$. 
A typical *MeDeCom* experiment includes testing a grid of values for $k$ and $\lambda$. For each combination 
of parameter values *MeDeCom* estimates a cross-validation error. The latter helps select the optimal number 
of LMCs and the strength of regularization.

# Installation

*MeDeCom* can be installed directly from github using the package `devtools`:


```r
devtools:::install_github("lutsik/MeDeCom")
```

Currently only *nix-like platforms with a C++11-compatible compiler are supported.
MeDeCom uses stack model for memory to accelerate factorization for smaller ranks. 
This requires certain preparation during compilation, therefore, please, note the extended
compilation time (15 to 20 minutes).

# Data preparation

*MeDeCom* accepts DNA methylation data in several forms. Preferably the user may load and preprocess the data
using a general-purpose DNA methylation analysis package [RnBeads](http://rnbeads.mpi-inf.mpg.de). A resulting RnBSet object
can be directly supplied to *MeDeCom*. Alternatively, *MeDeCom* runs on any matrix of type `numeric` with valid methylation values.

*MeDeCom* comes with a small example data set obtained by mixing reference profiles of blood cell methylomes *in silico*.
The example data set can be loaded in a usual way:


```r
## load the package
suppressPackageStartupMessages(library(MeDeCom))
##  load the example data sets
data(example.dataset, package="MeDeCom")
## you should get objects D, Tref and Aref
## in your global R environment
ls()
```

```
## [1] "Aref" "D"    "Tref"
```

Loaded numeric matrix `D` contains 100 *in silico* mixtures and serves as an example input. Columns of matrix `Tref` contains the methylomes 
of 5 blood cell types used to generate the mixtures, while matrix `Aref` provides the mixing proportions.

```r
## matrix D has dimension 10000x100
str(D)
```

```
##  num [1:10000, 1:100] 0.0354 0.0491 0.3411 0.873 0.1781 ...
```

```r
## matrix Tref has dimension 10000x5
str(Tref)
```

```
##  num [1:10000, 1:5] 0.1 0.0847 0.3167 0.834 0.0455 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : NULL
##   ..$ : chr [1:5] "Neutrophils" "CD4+ T-cells" "CD14+ Monocytes" "CD8+ T-cells" ...
```

```r
## matrix Aref has dimension 5x100
str(Aref)
```

```
##  num [1:5, 1:100] 6.54e-01 3.80e-02 3.06e-01 1.56e-13 2.45e-03 ...
```
# Performing a methylome decomposition experiment

*MeDeCom* can be run directly on matrix `D`. 

It is crucial to select the values of parameters $k$ and $\lambda$ to test. 
A choice of $k$ is often dictated by prior knowledge about the methylomes.
Precise value of lambda has to be selected for each data set independently. A good start 
is a logarithmic grid of lambda values. It is important to include $\lambda=0$ into the 
grid, as this particular case the regularization is effectively absent making *MeDeCom* similar 
to other NMF-based deconvolution algorithms.

```
medecom.result<-runMeDeCom(D, Ks=2:10, lambdas=c(0,10^(-5:-1)))
```

*MeDeCom* is based upon an alternating optimization heuristic and requires a lot 
of computation. The processing of the data matrix can take several hours.
One can speed up the run by decreasing the number of cross-validation folds and random initializations, and 
increasing the number of computational cores.


```r
medecom.result<-runMeDeCom(D, 2:10, c(0,10^(-5:-1)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=9)
```


```
## 
## [Main:] checking inputs
## [Main:] preparing data
## [Main:] preparing jobs
## [Main:] 3114 factorization runs in total
## [Main:] runs 2755 to 2788 complete
## [Main:] runs 2789 to 2822 complete
## [Main:] runs 2823 to 2856 complete
## [Main:] runs 2857 to 2890 complete
## ......
## [Main:] finished all jobs. Creating the object
```

This can, however, lead to decomposition slightly different from the one presented given below.

The results of a decomposition experiment are saved to an object of class `MeDeComSet`.
The contents of an object can be conveniently displayed using the `print` functionality.


```r
medecom.result
```

```
## An object of class MeDeComSet
## Input data set:
## 	10000 CpGs
## 	100 methylomes
## Experimental parameters:
## 	k values: 2, 3, 4, 5, 6, 7, 8, 9, 10
## 	lambda values: 0, 1e-05, 1e-04, 0.001, 0.01, 0.1
```

# Exploring the decomposition results

## Parameter selection 

The first key step is parameter selection. It is important to carefully explore the obtained results and make a decision about 
the most feasible parameter values, or about extending the parameter value grids to be tested in refinement experiments.

*MeDeCom* provides a **cross-validation error** (CVE) for each tested parameter combination.


```r
plotParameters(medecom.result)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

A lineplot helping to select parameter $\lambda$ can be produced by specifying a fixed value for $k$:


```r
plotParameters(medecom.result, K=5, lambdaScale="log")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

Cross-validation error has a minimum at $\lambda=10^{-2}$ so this value is preferred.

## Latent methylation components (LMCs)

A matrix of LMCs can be extracted using `getLMCs`:


```r
lmcs<-getLMCs(medecom.result, K=5, lambda=0.01)
str(lmcs)
```

```
##  num [1:10000, 1:5] 0 0.0184 0.2898 0.7856 0 ...
```

LMCs can be seen as measured methylation profiles of purified cell populations. 
*MeDeCom* provides for several visualization methods for LMCs using the function `plotLMCs` 
which operates directly on `MeDeComSet` objects.

### Clustering

For instance, standard hierarchical clustering can be visualized using:

```r
plotLMCs(medecom.result, K=5, lambda=0.01, type="dendrogram")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

A two-dimensional embedding with MDS is also obtainable:


```r
plotLMCs(medecom.result, K=5, lambda=0.01, type="MDS")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

Input data can be included into the MDS plot to enhance the interpretation.


```r
plotLMCs(medecom.result, K=5, lambda=0.01, type="MDS", D=D)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

### Matching LMCs to reference profiles

In many cases reference methylomes exists, which are relevant for the data set in question.
For our example analysis matrix `Tref` contains the reference type 
profiles which were *in silico* mixed. *MeDeCom* offers several ways to visualize the 
resulting LMCs together with the reference methylation profiles. The reference methylomes 
can be included into a joint clustering analysis: 


```r
plotLMCs(medecom.result, K=5, lambda=0.01, type="dendrogram", Tref=Tref, center=TRUE)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

Furthermore, a similarity matrix of LMCs vs reference profiles can be visualized as a heatmap.


```r
plotLMCs(medecom.result, K=5, lambda=0.01, type="heatmap", Tref=Tref)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

Correlation coefficient values and asterisks aid the interpretation.
The values are displayed in the cells which contain maximal values column-wise.
The asterisks mark cells which have the highest correlation value in the respective rows.
Thus, a value with asterisk corresponds to a mutual match, i.e. LMC unambiguously 
matching a reference profile.

In this example analysis each LMC uniquely matches one of the reference 
profiles. The matching of 

Function `matchLMCs` offers several methods for 
matching LMCs to reference profiles.


```r
perm<-matchLMCs(lmcs, Tref)
```

### LMC enrichment analysis

MeDeCom provides functions to perform enrichment analysis on the sites that are particularly hypo-/hypermethylated in an LMC. These sites can then be used for GO and LOLA enrichment analysis. Importantly, genomic annotations of the LMC sites is required to be specified. We thus recommend to use the DecompPipelie [https://github.com/CompEpigen/DecompPipeline](https://github.com/CompEpigen/DecompPipeline) for processing, but the annotation can also be specified manually using a ```data.frame``` that looks as follows:


```r
      Chromosome   Start     End Strand CpG GC CGI Relation SNPs
30365       chr1 1036375 1036376      +   2 59        Shelf <NA>
42681       chr1 1184537 1184538      +   2 58     Open Sea <NA>
45091       chr1 1218625 1218626      +   5 66       Island <NA>
51615       chr1 1292773 1292774      +   3 64       Island <NA>
52001       chr1 1295504 1295505      +   9 65        Shore <NA>
52003       chr1 1295507 1295508      +   9 65        Shore <NA>
```

The required columns are `Chromosome`, `Start`, `End`, and `Strand`. Using this ```data.frame``` (called `df` in the following), enrichment analysis can be performed using:


```r
lmc.lola.enrichment(medecom.result,anno.data=df,K=5,lambda=0.001,diff.threshold = 0.5, region.type = "tiling")
```

Please note that `df` needs to have the same number of rows than the methylation matrix used as input to MeDeCom. CpGs are first aggregated over the `region.type` specified, then the regions are selected that have a difference larger than `diff.threshold`. The list of available region types is published here [https://rnbeads.org/regions.html](https://rnbeads.org/regions.htm).

## Mixing proportions

A matrix of mixing proportions is obtained using `getProportions`:


```r
prop<-getProportions(medecom.result, K=5, lambda=0.001)
str(prop)
```

```
##  num [1:5, 1:100] 0.3411 0 0.065 0.0213 0.5726 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:5] "LMC1" "LMC2" "LMC3" "LMC4" ...
##   ..$ : NULL
```

### Visualization of the complete proportion matrix

A complete matrix of propotions can be visualized as a stacked barplot:

```r
plotProportions(medecom.result, K=5, lambda=0.01, type="barplot")
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

or a heatmap:


```r
plotProportions(medecom.result, K=5, lambda=0.01, type="heatmap")
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png)

The heatmap can be enhanced by clustering the columns:


```r
plotProportions(medecom.result, K=5, lambda=0.01, type="heatmap", heatmap.clusterCols=TRUE)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

or adding color code for the samples:


```r
sample.group<-c("Case", "Control")[1+sample.int(ncol(D))%%2]
plotProportions(medecom.result, K=5, lambda=0.01, type="heatmap", sample.characteristic=sample.group)
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

### Visualization of selected LMC proportions




```r
plotProportions(medecom.result,  K=5, lambda=0.01, type="lineplot", lmc=2, Aref=Aref, ref.profile=2)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

# Advanced usage

## Running *MeDeCom* on a compute cluster

*MeDeCom* experiments require a lot of computational time. On the other hand most of the factorization runs are 
independent and, therefore, can be run in parallel. Thus, a significant speedup can be achieved when running *MeDeCom* 
in an HPC environment. *MeDeCom* can be easily adapted to most of the popular schedulers. There are, however, several prerequisites:

 * the scheduler provides the standard utilities `qsub` for the submission of the cluster jobs and `qstat` for obtaining the job statistics;
 * the cluster does not have a low limit on the number of submitted jobs;
 * the R installation (location of the R binary and the package library) is consistent across the execution nodes.

The example below 
is for the cluster operated by *Son of Grid Engine* (SoGE). To be able to run on a SoGE cluster *MeDeCom* needs to know:

 * location of the R executable (directory);
 * an operating memory limit per each factorization job;
 * a pattern for the names of cluster nodes to run the jobs on.

These settings should be stored in a `list` object:

```r
sge.setup<-list(
R_bin_dir="/usr/bin",
host_pattern="*",
mem_limit="5G"
)
```
This object should be supplied to *MeDeCom* as the argument `cluster.settings`. It is also important to specify a valid temporary 
directory, which is available to all execution nodes.

```r
medecom.result<-runMeDeCom(D, Ks=2:10, lambdas=c(0,10^(-5:-1)), N_COMP_LAMBDA=1, NFOLDS=5, NINIT=10, 
temp.dir="/cluster_fs/medecom_temp",
cluster.settings=sge.setup)
```
*MeDeCom* will start the jobs and will periodically monitor the number of remaining ones.

```
## 
## [Main:] checking inputs
## [Main:] preparing data
## [Main:] preparing jobs
## [Main:] 3114 factorization runs in total
## [Main:] 3114 jobs remaining
## ....
## [Main:] finished all jobs. Creating the object
```

# R session
Here is the output of `sessionInfo()` on the system on which this document was compiled:

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 10 (buster)
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
## LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] MeDeCom_1.0.0                           RnBeads_2.4.0                          
##  [3] plyr_1.8.5                              methylumi_2.32.0                       
##  [5] minfi_1.32.0                            bumphunter_1.28.0                      
##  [7] locfit_1.5-9.1                          iterators_1.0.12                       
##  [9] foreach_1.4.7                           Biostrings_2.54.0                      
## [11] XVector_0.26.0                          SummarizedExperiment_1.16.1            
## [13] DelayedArray_0.12.2                     BiocParallel_1.20.1                    
## [15] FDb.InfiniumMethylation.hg19_2.2.0      org.Hs.eg.db_3.10.0                    
## [17] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 GenomicFeatures_1.38.1                 
## [19] AnnotationDbi_1.48.0                    reshape2_1.4.3                         
## [21] scales_1.1.0                            Biobase_2.46.0                         
## [23] illuminaio_0.28.0                       matrixStats_0.55.0                     
## [25] limma_3.42.1                            gridExtra_2.3                          
## [27] ggplot2_3.2.1                           fields_10.2                            
## [29] maps_3.3.0                              spam_2.5-1                             
## [31] dotCall64_1.0-0                         ff_2.2-14                              
## [33] bit_1.1-15.1                            cluster_2.1.0                          
## [35] MASS_7.3-51.5                           GenomicRanges_1.38.0                   
## [37] GenomeInfoDb_1.22.0                     IRanges_2.20.2                         
## [39] S4Vectors_0.24.3                        BiocGenerics_0.32.0                    
## [41] RUnit_0.4.32                            gplots_3.0.1.1                         
## [43] gtools_3.8.1                            pracma_2.2.9                           
## [45] Rcpp_1.0.3                              knitr_1.27                             
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.4-1         siggenes_1.60.0          mclust_5.4.5             base64_2.0              
##  [5] bit64_0.9-7              xml2_1.2.2               splines_3.6.3            codetools_0.2-16        
##  [9] scrime_1.3.5             Rsamtools_2.2.1          annotate_1.64.0          dbplyr_1.4.2            
## [13] HDF5Array_1.14.2         readr_1.3.1              compiler_3.6.3           httr_1.4.1              
## [17] assertthat_0.2.1         Matrix_1.2-18            lazyeval_0.2.2           prettyunits_1.1.1       
## [21] tools_3.6.3              gtable_0.3.0             glue_1.3.1               GenomeInfoDbData_1.2.2  
## [25] dplyr_0.8.4              rappdirs_0.3.1           doRNG_1.8.2              vctrs_0.2.2             
## [29] multtest_2.42.0          nlme_3.1-144             preprocessCore_1.48.0    gdata_2.18.0            
## [33] rtracklayer_1.46.0       DelayedMatrixStats_1.8.0 xfun_0.12                stringr_1.4.0           
## [37] lifecycle_0.1.0          rngtools_1.5             XML_3.99-0.3             beanplot_1.2            
## [41] zlibbioc_1.32.0          hms_0.5.3                GEOquery_2.54.1          rhdf5_2.30.1            
## [45] RColorBrewer_1.1-2       curl_4.3                 memoise_1.1.0            biomaRt_2.42.0          
## [49] reshape_0.8.8            stringi_1.4.5            RSQLite_2.2.0            genefilter_1.68.0       
## [53] highr_0.8                caTools_1.17.1.1         rlang_0.4.4              pkgconfig_2.0.3         
## [57] bitops_1.0-6             nor1mix_1.3-0            evaluate_0.14            lattice_0.20-40         
## [61] purrr_0.3.3              Rhdf5lib_1.8.0           GenomicAlignments_1.22.1 tidyselect_1.0.0        
## [65] magrittr_1.5             R6_2.4.1                 DBI_1.1.0                pillar_1.4.3            
## [69] withr_2.1.2              survival_3.1-8           RCurl_1.98-1.1           tibble_2.1.3            
## [73] crayon_1.3.4             KernSmooth_2.23-16       BiocFileCache_1.10.2     progress_1.2.2          
## [77] data.table_1.12.8        blob_1.2.1               digest_0.6.23            xtable_1.8-4            
## [81] tidyr_1.0.2              openssl_1.4.1            munsell_0.5.0            quadprog_1.5-8          
## [85] askpass_1.1
```

 
