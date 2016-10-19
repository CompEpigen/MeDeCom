#######################################################################################################################
##
##  Matching the recovered latent components to the true (experimentally obtained) ones 
##
#######################################################################################################################

match.matrix<-function(TT, Tref, method="pearson"){
	
	if(method=="pearson"){
		
		matchmat<-cor(TT,Tref, use="pairwise.complete.obs")
		
	}else if(method=="spearman"){
		
		matchmat<-cor(TT,Tref, method="spearman", use="pairwise.complete.obs")
		
	}else if(method=="anova"){
		
		fitsmat<-matrix(0L, nrow=ncol(TT), ncol=ncol(Tref))
		signsmat<-matrix(1L, nrow=ncol(TT), ncol=ncol(Tref))
		
		for(tt in 1:ncol(TT)){
			fits<-list()
			signs<-integer()
			for(i in 1:ncol(Tref)){
				signsmat[tt,i]<-sign(mean(Tref[TT[,tt]==range(TT)[2],i])-mean(Tref[TT[,tt]==range(TT)[1],i]))
				fitsmat[tt,i]<-summary(aov(ct~comp1, data=data.frame(ct=Tref[,i], comp1=TT[,tt])))[[1]]$`F value`[1]	
			}
		}
		
		matchmat<-signsmat*fitsmat
	}
}

#corrmatch<-function(TT, Tref, method="pearson",return.vals=F){
#	
#	if(method=="pearson"){
#		
#		values<-apply(cor(TT,Tref, use="pairwise.complete.obs"),1,max)
#		indices<-apply(cor(TT,Tref, use="pairwise.complete.obs"),1,which.max)
#		
#	}else if(method=="spearman"){
#		
#		values<-apply(cor(TT,Tref, use="pairwise.complete.obs"),1,max)
#		indices<-apply(cor(TT,Tref, method="spearman", use="pairwise.complete.obs"),1,which.max)
#		
#	}else if(method=="anova"){
#		
#		indices<-integer()
#		values<-numeric()
#		for(tt in 1:ncol(TT)){
#			fits<-list()
#			signs<-integer()
#			for(i in 1:ncol(Tref)){
#				signs<-c(signs, sign(mean(Tref[TT[,tt]==range(TT)[2],i])-mean(Tref[TT[,tt]==range(TT)[1],i])))
#				fits[[length(fits)+1]]<-aov(ct~comp1, data=data.frame(ct=Tref[,i], comp1=TT[,tt]))	
#			}
#			fvals<-sapply(fits, function(fit) summary(fit)[[1]]$`F value`[1])
#			fvals<-fvals*signs
#			indices<-c(indices,(which.max(fvals)))
#			values<-c(values,(max(fvals)))
#		}
#		
#	}
#	if(return.vals){
#		return(values)
#	}else{
#		return(indices)
#	}
#	
#}


corrmatch<-function(TT, Tref, method="pearson",return.vals=F){
	
	mmat<-match.matrix(TT,Tref,method)
	## a very bad solution for the degenerate cases
	mmat[is.na(mmat)]<-0 
	if(return.vals){
		return(apply(mmat, 1, max))
	}else{
		return(apply(mmat, 1, which.max))
	}
}

#######################################################################################################################
#'
#' greedymatch
#'
#' Matching of latent components 
#'
#'
#' @details suppose Tstar contains K* topics (columns)
#' 1. select the K* most popular topics of That
#' 2. use a greedy method to match them to those of Tstar (w.r.t L2 norm)
#'    and correspondingly permute rows of Ahat
#' 
#' @author Martin Slawski
#' @author R port by Pavlo Lutsik
#' 
#' @export
#' 
greedymatch<-function(Tstar,That,Ahat){
	
	### 1. select K* most popular topics of That, get submatrix T, A
	k <- size(That,2) ;
	kstar <- size(Tstar,2);
	
	if (k < kstar){
		# pad zeros rows and columns to T, A when k < K*
		Delta <- kstar - k;
		That <- cbind(That, matrix(0, nrow=nrow(That), ncol=Delta));
		Ahat <- rbind(Ahat, matrix(0, nrow=Delta, ncol=ncol(Ahat)));
	}
	
	pop <- rowSums(Ahat) / sum(Ahat)
	poporder <- order(pop, decreasing=TRUE)
	TT <- That[,poporder[1:kstar],drop=FALSE]
	A <- Ahat[poporder[1:kstar],,drop=FALSE]
	
	### matching submatrix
	idx <- 1:kstar
	Tm <- zeros(nrow(TT), ncol(TT))
	Am <- zeros(nrow(A), ncol(A))
	for(i in 1:kstar){
		
		topic <- Tstar[,i,drop=FALSE]
		dist <- zeros(length(idx),1)
		for(j in 1:length(idx)){
			dist[j,1] = norm(topic-TT[,idx[j],drop=FALSE],"2")
		}
		ii = which.min(dist)
		Tm[,i] = TT[,idx[ii]]
		Am[i,] = A[idx[ii],]
		idx<-idx[-ii]
		
	}
	
	return(list(Tm=Tm, Am=Am)) 
}

#######################################################################################################################
#'
#' match.components
#'
#' Matching the recovered components to the reference
#' 
#' @param That 		recovered components from factorization 
#' @param Tref		reference components 
#' @param method	one of "corrmatch" and "greedymatch"
#' @param Ahat 		for method "greedymatch": recovered mixing proportions
#' @param check		for method "corrmatch": a flag specifying whether to check uniquenes of the match
#' 
#' 
#' @details Wrapper function for component matching methods 
#' 
#' @return 			a vector of indices. The length and the order of the vector corresponds to the columns of \code{TT}
#' 					and the indices specify the columns of Tref
#' 
#' @export 
#' 
matchLMCs<-function(
		That, 
		Tref,
		method="corrmatch",
		check=TRUE,
		Ahat=NULL
		){
	
	if(method=="corrmatch"){
	
		mr<-corrmatch(That, Tref, method="pearson", return.vals=F)
		
		if(check && length(unique(mr))<ncol(Tref)){
			return(NULL)
		}else{
			return(mr)
		}
		
	}else if(method=="greedymatch"){
		
		if(is.null(Ahat)){
			stop("need estimated mixture components matrix for method greedymatch")
		}
		
		greedymatch(Tref, That, Ahat)
		
	}else{
		stop("match.components: this method has not been implemented yet")
	}
}

#######################################################################################################################

#
#  Legacy matching routine
#

assign.components<-function(TT, Tref, method="pearson",return.vals=F){
	
	if(method=="pearson"){
		values<-apply(cor(TT,Tref, use="pairwise.complete.obs"),1,max)
		indices<-apply(cor(TT,Tref, use="pairwise.complete.obs"),1,which.max)
	}else if(method=="spearman"){
		values<-apply(cor(TT,Tref, use="pairwise.complete.obs"),1,max)
		indices<-apply(cor(TT,Tref, method="spearman", use="pairwise.complete.obs"),1,which.max)
	}else if(method=="anova"){
		indices<-integer()
		values<-numeric()
		for(tt in 1:ncol(TT)){
			fits<-list()
			signs<-integer()
			for(i in 1:ncol(Tref)){
				signs<-c(signs, sign(mean(Tref[TT[,tt]==range(TT)[2],i])-mean(Tref[TT[,tt]==range(TT)[1],i])))
				fits[[length(fits)+1]]<-aov(ct~comp1, data=data.frame(ct=Tref[,i], comp1=TT[,tt]))	
			}
			fvals<-sapply(fits, function(fit) summary(fit)[[1]]$`F value`[1])
			fvals<-fvals*signs
			indices<-c(indices,(which.max(fvals)))
			values<-c(values,(max(fvals)))
		}
		
	}
	if(return.vals){
		return(values)
	}else{
		return(indices)
	}
	
}