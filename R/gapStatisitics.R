#######################################################################################################################
# 
#  Gap statistics for selection of parameter r (number of underlying components)
#  
#  author: Pavlo Lutsik
# 
#######################################################################################################################

#
# select.r.gap
#
# Procedure for selecting the best r using the modified gap statistics 
#
# @details  The gap statistics for a given value of parameter r 
#		    is calculated as a difference between the "null" RMSE of 
#		    and the RMSE of the actual data factorization.
#			The "null" RMSE is obtained by averaging the factorization RMSE of the 
#			the randomly permuted input matrix
#			
#
# @author Pavlo Lutsik
#
select.r.gap<-function(D, method, minr=1, maxr=5, nperm=100, return.all=TRUE, plot=TRUE, ...){
		
	rmse.orig<-rep(NA, length(minr: maxr))
	rmse.perm<-rep(NA, length(minr: maxr))
	sd.rmse.perm<-rep(NA, length(minr: maxr))
	
	perms<-lapply(1:nperm, function(np){
		
		sample.int(ncol(D)*nrow(D), ncol(D)*nrow(D))
		
	})
	
	for(ix in 1:length(minr:maxr)){
		
		kk<-(minr:maxr)[ix]
		##	Original factorization
		if(method %in% ALGORITHMS[2:6]){
						
			fr.res<-factorize.alternate(D, k=kk, t.method=T_METHODS[method],...)
			rmse.orig[ix]<-fr.res$rmse
			
		}
		
		curr.rmses<-rep(NA, length(nperm))
		for(pp in 1:nperm){
			
			fr.res<-factorize.alternate(matrix(D[perms[[pp]]], ncol=ncol(D)), k=kk, t.method=T_METHODS[method],...)
			curr.rmses[pp]<-fr.res$rmse
	
		}	
		
		rmse.perm[ix]<-mean(curr.rmses)
		sd.rmse.perm[ix]<-sd(curr.rmses)
	
	}

	
	if(plot){
		
		layout(matrix(1:2, ncol=2))
		x<-minr:maxr
		plot(x,rmse.perm, col="red", pch=0, type="l", ylim=c(0, max(c(rmse.orig, rmse.perm))), main="RMSE")
		segments(x, rmse.perm-sd.rmse.perm,x, rmse.perm+sd.rmse.perm, col="red")
		epsilon = 0.02
		segments(x-epsilon,rmse.perm+sd.rmse.perm,x+epsilon,rmse.perm+sd.rmse.perm, col="red")
		segments(x-epsilon,rmse.perm-sd.rmse.perm,x+epsilon,rmse.perm-sd.rmse.perm, col="red")
		
		lines(rmse.orig, col="green", pch=1, type="o")
		
		
		plot(rmse.perm-rmse.orig, col="black", pch=1, type="o", main="Gap")
		
		
	}
	
	if(return.all){
	
		return(list(Rs=minr:maxr, RMSE=rmse.orig, RMSE.perm=rmse.perm, sd.RMSE.perm=sd.rmse.perm))
		
	}else{
		
		return((minr:maxr)[which.max(rmse.perm-rmse.orig)])		
	}
}

#######################################################################################################################
