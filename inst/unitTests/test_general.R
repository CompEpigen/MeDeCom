test.reinius.reference <- function(){
  data("small.example.RnBSet")
  Ks <- 2
  lambdas <- c(0.01)
  res <- run.refbased(rnb.set = rnb.set.example,Ks = Ks,lambdas = lambdas)
  plot <- cluster.refbased(res,plot.type="dendrogram",K=2,lambda=0.01)
  passes <- is.null(plot)
  checkTrue(passes)
}

test.general <- function(){
  data("example.dataset")
  input.data <- D[sample(1:nrow(D),1000),sample(1:ncol(D),5)]
  cg_subsets <- list("var"=sample(1:nrow(D),250),"random"=sample(1:nrow(D),500))
  Ks <- 2
  lambdas <- c(0.01)
  res <- runMeDeCom(data = D,Ks = Ks,lambdas = lambdas,cg_subsets = cg_subsets)
  passes <- inherits(res,"MeDeComSet")
  checkTrue(passes)
}

test.contribution.interpretation <- function(){
  data("example.MeDeComSet")
  anno.frame <- data.frame(Sex=sample(c("M","F"),100,replace=T),Age=sample(1:100,100,replace=T),Ethnicity=sample(c("A","B","C"),100,replace = T))
  res <- run.trait.association.single(medecom.result,pheno.data=anno.frame)
  passes <- all(names(res) %in% c("linear model","qualitative","quantitative"))
  checkTrue(passes)
}

# test.enrichment <- function(){
#   require("RnBeads")
#   data("example.MeDeComSet")
#   anno.frame <- rnb.annotation2data.frame(rnb.get.annotation("probes450"))[sample(1:460000,10000),]
#   res <- lmc.lola.plots.tables(medecom.result,anno.data=anno.frame)
#   passes <- all(names(res) %in% c("Plots","Tables"))
#   res <- lmc.go.enrichment(medecom.result,anno.data=anno.frame)
#   passes <- passes && (class(res) == "list")
#   checkTrue(passes)
# }

test.routine <- function(){
  require("RUnit")
  require("MeDeCom")
  cat("STARTED testing general function \n")
  test.general()
  cat("COMPLETED testing general function \n")
  cat("STARTED testing contribution interpretation \n")
  test.contribution.interpretation()
  cat("COMPLETED testing contribution interpretation \n")
  # cat("STARTED testing enrichment functions \n")
  # test.enrichment()
  # cat("COMPLETED testing enrichment functions \n")
}

test.routine()