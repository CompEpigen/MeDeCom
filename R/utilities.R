###############################################################################
#  
# Utility routines 
#  
# created 2012-10-11
# 
# Author: Pavlo Lutsik
###############################################################################

#
# rmse
#
# Root mean square error
#
# D data matrix
# T recovered profiles 
# A recovered proportions
#
#
rmse<-function(D,T,A){
	
	sqrt(sum((D-T%*%A)^2)/nrow(D)/ncol(D))
	
}

###############################################################################

# logLik2
#
# Calculate log-likelihood of a fitted linear regression model
# 
# @param vector of the model residuals
# @return log-likelihood value
# 
# @author Pavlo Lutsik
logLik2<-function(res){
	N<-length(res)
	val <- 0.5 * (- N * (log(2 * pi) + 1 - log(N) + log(sum(res^2))))
	return(val)
}

###############################################################################
#wilcoxTestTwins<-function(betas, twin1, twin2){
#	
#	deltas<-betas[,twin1]-betas[,twin2]
#	wilcox<-apply(deltas, 1, wilcox.test)
#	results<-data.frame(meanDeltas=rowMeans(deltas), pval=sapply(wilcox, simplify=T, function(t) return(t$p.value)))
#	
#}
#
###############################################################################
#tTestTwins<-function(betas, twin1, twin2){
#	
#	deltas<-betas[,twin1]-betas[,twin2]
#	tt<-apply(deltas, 1, t.test)
#	results<-data.frame(meanDeltas=rowMeans(deltas), pval=sapply(tt, simplify=T, function(t) return(t$p.value)))
#	
#}
###############################################################################
# fetch.annotations
#
# collects annotations for all probe ids in rownames(dataset) from IlluminaHumanMethylation450k.db.
# @author Pavlo Lutsik
# @param dataset a matrix/data.frame object all row names of which are valid Infinium probe ids (cgXXXXXXXX).
# @param tables a character vector of IlluminaHumanMethylation450k.db tables. If NULL all existing annotations are fetched
# @return a data.frame object with annotations
# @seealso \code{\link{IlluminaHumanMethylation450k.db}}
# @export

fetch.annotations<-function(dataset, ids, tables=NULL, verbose=F)
{
	require(IlluminaHumanMethylation450k.db)
	if(is.null(tables)){
		all.tables<-ls("package:IlluminaHumanMethylation450k.db")
		meta.tables<-all.tables[c(14:65)]
	}else{
		meta.tables<-tables
	}
	
	#metaTables<-allTables[grep("IlluminaHumanMethylation450k[A-Z]", allTables, perl=T)]
	if(!missing(dataset)) ids<-rownames(dataset)
	all.results<-data.frame(ID=ids)
	cn<-colnames(all.results)
	
	
	for (table.name in meta.tables) {
#		table2 <- get(table.name)

		do.call("<-", list(as.name("table2"), as.name(table.name)))
		if(!class(table2) %in% c("function","integer") ) {
			if(length(table2)>1){
				if(verbose) print(paste("Extracting data from:", table.name))
				present<-intersect(ids,keys(table2))
				result<-data.frame(ID=present, V=as.character(as.list(table2[present])))
				cn<-c(cn,strsplit(table.name,"IlluminaHumanMethylation450k")[[1]][2])
				all.results<-merge(all.results, result, by="ID", all.x=T, sort=F)
				colnames(all.results)<-cn
			}
		}
		
	}
	
	
	return(all.results[match(ids,all.results$ID),])
	
}
###############################################################################
fetch.by.interval<-function(chromosome, start, end, strand="*", type="probes450"){
	
	if(!type %in% c("probes450", "hg19", "mm9")){
		stop("unsupported annotation type")
	}
	
	ann<-rnb.get.annotation(type)
	ann<-ann[match(chromosome,names(ann))]
	
	interval<-GRanges(chromosome,IRanges(start=start, end=end), strand=strand)
	
	olap<-findOverlaps(ann[[1]], interval)
	
	return(rownames(rnb.annotation2data.frame(ann))[queryHits(olap)])
		
}
###############################################################################
combin<-function(V,k){
#	
#	if(k>1){
#		vv<-rep(V,k)
#		Posss<-apply(combn(length(vv),k),2,function(idx) vv[idx])
#		Posss<-Posss[,!duplicated(t(Posss)), drop=FALSE]
#	}else{
#		return(matrix(V,ncol=length(V)))
#	}
#	
#	print("##################### DIMENSIONS ####################")
#	print(dim(Posss))
#	return(Posss)
	
	t(expand.grid(list(V)[rep(1,k)]))
		
}

###############################################################################

randsplxmat<-function(m,n){
	
	U = -log(matrix(runif(m*n), nrow=m))
	S = colSums(U) # probably only 1 row, so specify 1 explicitly 
	dump<-sapply(1:n, function(j) U[,j] <<- U[,j]/S[j])   
	return(U)
	
}

###############################################################################
projectV<-function(TP,V){
	
	TP_proj<-apply(TP, 2, function(col){
				dist<-abs(repmat(matrix(V, ncol=length(V)), n=length(col), m=1)-col)								
				proj<-V[apply(dist,1,which.min)]
			})
	TP_proj
	
}
###############################################################################

opt.factorize.exact<-function(affine=TRUE,
		nonnegative=TRUE,
		aggr="no",
		replace=1L,
		nsamples=0L,
		naggsets=5L,
		chunksize=10L,
		verbose=FALSE,
		varargin=NULL){
	
	options<-list(
			affine=affine,
			nonnegative=nonnegative,
			aggr=aggr,
			replace=replace,
			nsamples=nsamples,
			naggsets=naggsets,
			chunksize=chunksize,
			verbose=verbose)
	
	if(is.null(varargin)){
		
		return(options)
		
	}else{
		print("not supported yet")
		
		return(options)
	}
	
}
###############################################################################

licols<-function(X,tol=1e-10){
	
	if(all(X==0)){ #X has no non-zeros and hence no independent columns
		
		return(list(r<-integer(), Xsubs=matrix(0), idx<-integer()))
	}
	
	colnames(X)<-as.character(1:ncol(X))
	
	qr.res<-qr(X, tol = tol, LAPACK=T)
	Q<-qr.Q(qr.res)
	R<-qr.R(qr.res)
	E<-match(colnames(R), colnames(X))
	
	diagr <- abs(diag(R));
	
	r <-length(which(diagr >=tol*diagr[1])) #Rank estimation
	
	idx<-sort(E[1:r])
	
	Xsub<-X[,idx]
	
	return(list(r=r, idx=idx, Xsub=Xsub))
}

adjust.nmf.coefs<-function(coefs){
	t(t(coefs)/colSums(coefs))
}

###############################################################################
#
# Down-rank SNP-affected probes using a one-dimentional k-means
# with k=3
#

rank.snp<-function(D){
	
	require(Ckmeans.1d.dp)
	
	R3 <- numeric(nrow(D))
	for(i in 1:nrow(D)){
		res <- Ckmeans.1d.dp(D[i,], 3)
		R3[i] <- sum(res$withinss)
	}
	
	R3
}

###############################################################################
get.cpg.intervals<-function(cpg.coords, ids=NULL, offset=50){
	
	#if(annot.full<-rnb.get.annotation("probes450")
	intervals<-GRanges(seq=cpg.coords$Chromosome, 
			IRanges(start=cpg.coords$Start, width=1), 
			strand=rep("*", nrow(cpg.coords)))
	
	
	intervals<-flank(intervals, width=offset, both=TRUE)
	intervals<-reduce(intervals)
	as.data.frame(intervals)
	
}
###############################################################################
#
# generateExample
#
# Examples for testing factorization methods
# 
# m	number of genomic features
# n number of profiles
# k hidden dimension
# t.method method used to generate matrix T: "integer", "uniform" or "beta" 
# a.method method used to generate matrix A: "uniform" or "dirichlet" 
# V a vector of possible values for \cs{method} integer or the upper and lower 
# 			bounds for the \cs{method} uniform 
# beta1 first beta-distribution parameter for \cs{t.method} beta
# beta2 second beta-distribution parameter for \cs{t.method} beta
# proportion.prior numeric vector of length \code{r}
# noise.sd a standard deviation for additive Gaussian noize
# digits desired precision
#
# return \cs{list} with elements \cs{D}, \cs{T} and \cs{A}
#
#
generateExample<-function(
		m, 
		n, 
		k, 
		t.method="beta",
		a.method="dirichlet",
		e.method="gaussian",
		V=c(0,1), 
		beta1=0.5, 
		beta2=0.5,
		proportion.prior=NULL,
		proportion.var.facror=1,
		Alower=rep(0, k),
		Aupper=rep(1, k),
		noise.sd=0,
		digits=12
){
	
	if(t.method=="integer"){
		
		Tt<-matrix(Inf,ncol=k, nrow=m)
		ri<-1
		it<-0
		while(any(Tt==Inf) && it<100){
			it<-it+1
			vertex<-V[sample.int(length(V), m, replace=T)]
			if(all(colSums(abs(Tt-vertex))!=0)){
				Tt[1:m,ri]<-vertex
				ri<-ri+1
			}
		}
		
	}else if(t.method=="uniform"){
		
		Tt<-matrix(runif(m*k, min=min(V), max=max(V)), ncol=k)
		
	}else if(t.method=="beta"){
		
		Tt<-matrix(rbeta(m*k, shape1=beta1, shape2=beta2), ncol=k)
		
	}else{
		stop("this method for generating T is not implemented")
	}
	
	if(a.method=="uniform"){
		
		if(max(Alower)==0 && min(Aupper)==1){
			
			A<-matrix(-log(runif(k*n)), ncol=n, nrow=k)
			A<-t(t(A)/colSums(A))
			A[1,]<-A[1,]+(rep(1,n)-colSums(A))
			
		}else{
			A <- t(sapply(1:k, function(pri){
								-log(runif(n, min=Alower[pri], max=Aupper[pri]))
							}))
			A<-t(t(A)/colSums(A))
			A[1,]<-A[1,]+(rep(1,n)-colSums(A))
	#		for(kk in 1:n){
	#			A[,kk] <- RProjSplxBox(A[,kk,drop=FALSE], Alower, Aupper);
	#		
		}
		
	}else if(a.method=="dirichlet"){
		
		if(is.null(proportion.prior)){
			proportion.prior<-rep(1/k,k)
		}
		A<-t(rdirichlet(n, proportion.prior*proportion.var.facror))
		
	}else{
		stop("this method for generating A is not implemented")
	}
	
	if(e.method=="gaussian"){
		if(noise.sd>0){
			E=matrix(rnorm(m*n, sd=noise.sd), ncol=n)
		}else{
			E=matrix(0, nrow=m, ncol=n)
		}
	}else{
		stop("this method for generating additive noise is not implemented")
	}
	
	# get the data matrix
	Tt<-round(Tt, digits=digits)
	A<-round(A, digits=digits)
	E<-round(E, digits=digits)
	
	D<-Tt%*%A+E
	
	D[D<min(V)]<-min(V)
	D[D>max(V)]<-max(V)
	
	return(list("T"=Tt,"A"=A,"D"=D, "E"=E))
}

###############################################################################

RMSE_T<-function(That, Tstar, perm){
	
	if(length(unique(perm))==length(perm)){
		rmseT<-sqrt(sum((That-Tstar[,perm])^2)/ncol(Tstar)/nrow(Tstar))
	}else{
		rmseT<-sqrt(mean((
					sapply(unique(perm), function(comp){
								abs(Tstar[,comp]-rowMeans(That[,perm==comp,drop=FALSE]))
					}))^2))
	}
	
	rmseT
}
###############################################################################
MAE_A<-function(Ahat, Astar, perm){
	
	if(length(unique(perm))==length(perm)){
		maeA<-sum(abs(Ahat-Astar[perm,]))/ncol(Astar)/nrow(Astar)
	}else{
#		maeA<-sum(abs(
#						sapply(unique(perm), function(comp){
#								abs(Astar[comp,,drop=FALSE]/colSums(Astar[unique(perm),,drop=FALSE])-colSums(Ahat[perm==comp,,drop=FALSE]))
#						})))/ncol(Astar)/nrow(Astar)

		aggrAhat<-t(sapply(unique(perm), function(comp){
					colSums(Ahat[perm==comp,,drop=FALSE])
				}))
		
		aggrAstar<-sweep(Astar[unique(perm),,drop=FALSE], 2, colSums(Astar[unique(perm),,drop=FALSE]), "/")
		
		maeA<-sum(abs(aggrAhat-aggrAstar))/nrow(aggrAstar)/ncol(aggrAstar)
	}
	maeA
}
###############################################################################
estimate.accuracy<-function(fr, trueT, trueA, check=FALSE){
	
	perm<-MeDeCom:::match.components(fr$T, trueT, check=check)
	
	if(!is.null(perm)){
		if(length(unique(perm))==length(perm)){
			rmseT<-sqrt(sum((trueT-fr$T[,perm])^2)/ncol(trueT)/nrow(trueT))
			maeA<-sum(abs(trueA-fr$A[perm,]))/ncol(trueA)/nrow(trueA)
		}else{
			rmseT<-sqrt(sum((
										sapply(unique(perm), function(comp){
													abs(trueT[,comp]-rowMeans(fr$T[,perm==comp,drop=FALSE]))
												}))^2)/ncol(trueT)/nrow(trueT))
			
			maeA<-sum(
					sapply(unique(perm), function(comp){
								abs(trueA[comp,]-colMeans(fr$A[perm==comp,,drop=FALSE]))
							}))/ncol(trueA)/nrow(trueA)
		}
	}else{
		rmseT<-NA
		maeA<-NA
	}
	
	return(list(rmseT=rmseT, maeA=maeA))
}
###############################################################################
get.distance.matrix<-function(mdd, measure, centered=FALSE){
	
	if(centered){
		mdd<-lapply(mdd, function(mm) sweep(mm, 1, rowMeans(mm)))
	}
	mdd<-do.call("cbind", mdd)
	
	if(measure=="euclidean"){
		d <- dist(t(mdd))
	}else if(measure=="angular"){
		dm<-matrix(NA, ncol=ncol(mdd), nrow=ncol(mdd))
		colnames(dm)<-colnames(mdd)
		rownames(dm)<-colnames(mdd)
		for(ri in 1:ncol(mdd)){
			for(ci in 1:ncol(mdd)){
				dm[ri,ci]<-sum(mdd[,ci]*mdd[,ri])/sqrt(sum(mdd[,ci]^2))/sqrt(sum(mdd[,ri]^2))
			}
		}
		d <- as.dist(1-dm)
	}else if(measure=="correlation"){
		d <- as.dist(1-cor(mdd, method="pearson"))
	}
}
#####
RQuadHC_dummy<-function(G, W, Tk, tol, lower, upper){
	outp<-Tk+matrix(runif(ncol(Tk)*nrow(Tk))*tol*10, nrow(Tk),ncol(Tk))
	outp[outp>1]<-1;
	dummy_loss<-sqrt(mean((Tk-outp)^2))
	return(list(outp, dummy_loss))
}

###############################################################################
# Windows version of mclapply from 
# https://www.r-bloggers.com/implementing-mclapply-on-windows-a-primer-on-embarrassingly-parallel-computation-on-multicore-systems-with-r/
mclapply_win<-function(...){
    ## Create a cluster
    ## ... How many workers do you need?
    ## ... N.B. list(...)[[1]] returns the first 
    ##          argument passed to the function. In
    ##          this case it is the list to iterate over
    #size.of.list <- length(list(...)[[1]])
    print(str(list(...)))
    ncores<-list(...)[["NCORES"]]
    cl <- makeCluster(min(ncores, detectCores()) )
    
    ## Find out the names of the loaded packages 
    loaded.package.names <- c(
            ## Base packages
            sessionInfo()$basePkgs,
            ## Additional packages
            names( sessionInfo()$otherPkgs ))
    
    ## N.B. tryCatch() allows us to properly shut down the 
    ##      cluster if an error in our code halts execution
    ##      of the function. For details see: help(tryCatch)
    tryCatch( {
                
                ## Copy over all of the objects within scope to
                ## all clusters. 
                ## 
                ## The approach is as follows: Beginning with the 
                ## current environment, copy over all objects within
                ## the environment to all clusters, and then repeat
                ## the process with the parent environment. 
                ##
                this.env <- environment()
                while( identical( this.env, globalenv() ) == FALSE ) {
                    clusterExport(cl,
                            ls(all.names=TRUE, env=this.env),
                            envir=this.env)
                    this.env <- parent.env(environment())
                }
                ## repeat for the global environment
                clusterExport(cl,
                        ls(all.names=TRUE, env=globalenv()),
                        envir=globalenv())
                
                ## Load the libraries on all the clusters
                ## N.B. length(cl) returns the number of clusters
                parLapply( cl, 1:length(cl), function(xx){
                            lapply(loaded.package.names, function(yy) {
                                        ## N.B. the character.only option of 
                                        ##      require() allows you to give the 
                                        ##      name of a package as a string. 
                                        require(yy , character.only=TRUE)})
                        })
                
                ## Run the lapply in parallel 
                print("start_executing")
                return( parLapply( cl, ...) )
            }, finally = {        
                ## Stop the cluster
                stopCluster(cl)
            })
}

### END