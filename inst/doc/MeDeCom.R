## ---- eval=FALSE---------------------------------------------------------
#  devtools:::install_github("lutsik/MeDeCom")

## ------------------------------------------------------------------------
## load the package
library(MeDeCom)
##  load the example data sets
data(example.dataset, package="MeDeCom")
## you should get objects D, Tref and Aref
## in your global R environment
ls()

## ------------------------------------------------------------------------
## matrix D has dimension 10000x100
str(D)
## matrix Tref has dimension 10000x5
str(Tref)
## matrix Aref has dimension 5x100
str(Aref)

## ---- eval=FALSE---------------------------------------------------------
#  medecom.result<-runMeDeCom(D, 2:10, c(0,10^(-5:-1)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=9)

## ---- echo=FALSE---------------------------------------------------------
cat("
[Main:] checking inputs
[Main:] preparing data
[Main:] preparing jobs
[Main:] 3114 factorization runs in total
[Main:] runs 2755 to 2788 complete
[Main:] runs 2789 to 2822 complete
[Main:] runs 2823 to 2856 complete
[Main:] runs 2857 to 2890 complete
......
[Main:] finished all jobs. Creating the object
")
data(example.MeDeComSet)

## ------------------------------------------------------------------------
medecom.result

## ---- fig.width=7--------------------------------------------------------
plotParameters(medecom.result)

## ---- fig.width=5.5, fig.height=6----------------------------------------
plotParameters(medecom.result, K=5, lambdaScale="log")

## ------------------------------------------------------------------------
lmcs<-getLMCs(medecom.result, K=5, lambda=0.01)
str(lmcs)

## ------------------------------------------------------------------------
plotLMCs(medecom.result, K=5, lambda=0.01, type="dendrogram")

## ------------------------------------------------------------------------
plotLMCs(medecom.result, K=5, lambda=0.01, type="MDS")

## ------------------------------------------------------------------------
plotLMCs(medecom.result, K=5, lambda=0.01, type="MDS", D=D)

## ------------------------------------------------------------------------
plotLMCs(medecom.result, K=5, lambda=0.01, type="dendrogram", Tref=Tref, center=TRUE)

## ------------------------------------------------------------------------
plotLMCs(medecom.result, K=5, lambda=0.01, type="heatmap", Tref=Tref)

## ------------------------------------------------------------------------
perm<-matchLMCs(lmcs, Tref)

## ------------------------------------------------------------------------
prop<-getProportions(medecom.result, K=5, lambda=0.001)
str(prop)

## ------------------------------------------------------------------------
plotProportions(medecom.result, K=5, lambda=0.01, type="barplot")

## ---- fig.width=8, fig.height=6------------------------------------------
plotProportions(medecom.result, K=5, lambda=0.01, type="heatmap")

## ---- fig.width=8, fig.height=6------------------------------------------
plotProportions(medecom.result, K=5, lambda=0.01, type="heatmap", heatmap.clusterCols=TRUE)

## ---- fig.width=8, fig.height=6------------------------------------------
sample.group<-c("Case", "Control")[1+sample.int(ncol(D))%%2]
plotProportions(medecom.result, K=5, lambda=0.01, type="heatmap", sample.characteristic=sample.group)

## ---- echo=FALSE---------------------------------------------------------
rownames(Aref)<-colnames(Tref)

## ------------------------------------------------------------------------
plotProportions(medecom.result,  K=5, lambda=0.01, type="lineplot", lmc=2, Aref=Aref, ref.profile=2)

## ------------------------------------------------------------------------
sge.setup<-list(
R_bin_dir="/usr/bin",
host_pattern="*",
mem_limit="5G"
)

## ---- eval=FALSE---------------------------------------------------------
#  medecom.result<-runMeDeCom(D, Ks=2:10, lambdas=c(0,10^(-5:-1)), N_COMP_LAMBDA=1, NFOLDS=5, NINIT=10,
#  temp.dir="/cluster_fs/medecom_temp",
#  cluster.settings=sge.setup)

## ---- echo=FALSE---------------------------------------------------------
cat("
[Main:] checking inputs
[Main:] preparing data
[Main:] preparing jobs
[Main:] 3114 factorization runs in total
[Main:] 3114 jobs remaining
....
[Main:] finished all jobs. Creating the object
")

## ---- echo=FALSE---------------------------------------------------------
sessionInfo()

