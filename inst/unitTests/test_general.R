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
  Ks <- 2
  lambdas <- c(0.01)
  res <- runMeDeCom(data = D,Ks = Ks,lambdas = lambdas)
  passes <- inherits(res,"MeDeComSet")
  checkTrue(passes)
}

test.routine <- function(){
  require("RUnit")
  require("MeDeCom")
  cat("STARTED testing general function \n")
  test.general()
  cat("COMPLETED testing general function \n")
}

test.routine()