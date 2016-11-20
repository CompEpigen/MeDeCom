#######################################################################################################################
#
#  factorizations.R
#  
#  R implementation of matrix factorization algorithms with various constraints and regularization
#
#######################################################################################################################
#  Constrained optimization for obtaining the mixing proportions
#######################################################################################################################
#'
#' factorize.regr
#' 
#' Get mixing proportions from the target data matrix and a matrix of latent factors
#' 
#' @param D 		m by n matrix with mixture data
#' @param Tt		either a m by k matrix of k true latent components
#' 					or a list of n such matrices, one per each column of D
#' @param A0		initialization for the mixture proportions matrix 
#' @param precision numerical tolerance of the optimization algorithm
#' 
#' @return 			a \code{list} with elements:
#' 					\describe{
#' 						\item{\code{A}}{matrix of mixing proportions}
#' 						\item{\code{T}}{Tt used}
#' 						\item{\code{rmse}}{RMSE of the regression model}
#' 					}
#' 
#' @export
#' 
factorize.regr<-function(D, Tt, A0=NULL, precision=1e-8){
	
	if(is.null(A0)){
		if(is.matrix(Tt)){
			A0<-matrix(1/ncol(Tt), nrow=ncol(Tt), ncol=1)
		}else{
			A0<-matrix(1/ncol(Tt[[1]]), nrow=ncol(Tt[[1]]), ncol=1)
		}
	}
	
	Alphas<-matrix(NA_real_, nrow=if(is.list(Tt)) ncol(Tt[[1]]) else ncol(Tt), ncol=ncol(D))
	
	#sapply(1:ncol(D), function(i){
	for(i in 1:ncol(Alphas)){
				if(is.list(Tt)){
					Tprof<-Tt[[i]]
				}else{
					Tprof<-Tt
				}
				
				G <- t(Tprof) %*% Tprof
				k <- ncol(Tprof)
				
				W <-t(Tprof) %*% D[,i]
				alpha<-RQuadSimplex(G, W, A0, precision)
				Alphas[,i]<-alpha[["A"]]
			}#)
	
	if(is.list(Tt)){
		rmse<-0
		for(ni in 1:length(Tt)){
			rmse<-rmse+sqrt(sum((D[,ni,drop=FALSE]-Tt[[ni]]%*%Alphas[,ni,drop=FALSE])^2)/ncol(D)/nrow(D))
		}
	}else{
		rmse<-sqrt(sum((D-Tt%*%Alphas)^2)/ncol(D)/nrow(D))
	}
	
	return(list("A"=Alphas, "T"=Tt, "rmse"=rmse))
}


#######################################################################################################################
#		T update methods
#######################################################################################################################
#
# Combinatorial update by selecting rows of T from all possible ones
#
updateT_integer<-function(G, W, TT, Poss, lambda){

	# add equal penalty
	first <- t(colSums(Poss * (G %*% Poss), 1)) + lambda * t(colSums(Poss))
	second <- -2*t(Poss)%*%W
		
	#[~,minix]=min(repmat(first,1,size(second,2))+second);
	minix<-apply(rep(t(first), ncol(second)) + second, 2, which.min)

	Tnew <- t(Poss[,minix,drop=FALSE]);
	return(Tnew)
	
# for j=1:d
#     %allobjs = transpose(sum(Poss .* (G * Poss), 1)) - 2 * Poss' * W(:,j);
#     %allobjs = first - 2*Poss'*W(:,j);
#     %[~, minix] = min(allobjs);
#     [~,minix]=min(first+second(:,j));
#		%     Tnew(j, :) = Poss(:, minix); 
#		% end
#		%                
#		end

}
#######################################################################################################################
# 
#   updateT_gini
#  
#   Use D.C. algorithm to solve problem
#
#   min_{T}     norm(D-TA, 'fro')^2 + lambdaT * gini(T), where gini(T) = sum_{j,k} T_jk (1 - T_jk)
#   sb.to       0 <= T_jk <= 1 
#
#   This problem is equivalent to the problem  
#    
#   tr(T * G * T') - 2 * tr(W * T') + lambdaT * gini(T), 
# 
#   where G = A * A',  W = D * A'
#
#   In each D.C. iteration, a convex quadratic problem of the form
#
#    tr(T * G * T') - 2 * tr(W * T') + lambda_T * tr(g * T')
#    sb.to   0 <= T_jk <= 1 
#
#    is solved, where g is the gradient at the current iterate.
#
# -- Input -- 
#
#
#   G, W    - matrices defined above
#   Tk      - initial solution  
#   lambdaT - regularization parameter
#   tol     - tolerance for the stopping condition for DC
#   f0      - objective evaluated at the initial solution
#                 
#
# -- output --
#  Tk1      - the solution
#  fk1      - objective evaluated at the solution:
#             tr(Tk1 * G * Tk1') - 2 * tr(W * Tk1') + lambdaT * gini(Tk1).
#
#  original MATLAB code by Martin Slawski
#
updateT_gini<-function(G, W, Tk, lambdaT, tol, f0, lower=0, upper=1){
	
	# DC Step
	fk <- f0;
	iter  <- 0;
	while(TRUE){

		iter <- iter + 1;
		
		# gradient of h at Tk
		gh <- 1 - 2 * Tk; 
		
		# solve DC step with SPG
        qp.res <- RQuadHC(G, W - 0.5 * lambdaT * gh, Tk, tol, lower, upper);
		Tk1<-qp.res[[1]]; Loss_temp<-as.numeric(qp.res[[2]])
	
		# subtract the linear part that we get from the gini penalty
		Loss_new <- Loss_temp - lambdaT * sum(Tk1 * gh);
		
		# difference at regularizer part
		Reg_new <- lambdaT * sum(Tk * (1 - Tk));
	
		fk1 <- Loss_new + Reg_new;
		fk1_fk <- fk1 - fk;
		red_f <- abs(fk1_fk / fk);
		# relative reduction, must be absolute value
		
		# check stopping criterion
		if(fk1_fk >= 0){
			break;
		}else{
			# update
			Tk <- Tk1;
			fk <- fk1;

			if(red_f < tol){
				break;
			}
		}
	}

	return(list(Tk1, fk1, iter))
}


#######################################################################################################################

updateT_empirical<-function(G, W, Poss, lambda){
	
	first <- t(colSums(Poss * (G %*% Poss), 1)) + lambda * t(colSums(Poss * (1-Poss)))
	
	second <- -2*t(Poss)%*%W

	minix<-apply(rep(t(first), ncol(second)) + second, 2, which.min)
	
	Tnew <- t(Poss[,minix,drop=FALSE])
	
	return(Tnew)
	
}
#######################################################################################################################

#updateT_resample<-function(D, A, TT, lambda, lower=0, upper=1, vsf=1, nsamp=100){
#	
#	G <- A %*% t(A); W <- A %*% t(D);
#	
#	Poss.l<-lapply(1:nrow(TT), function(ii){
#		
#		means<-TT[ii,,drop=FALSE]
#		var.scale <- sum((D[ii,,drop=FALSE]-means%*%A)^2)/length(means)
#				
#		Poss.mat<-matrix(rtruncnorm(nsamp*length(means),
#						a=rep(lower, length(means)),
#						b=rep(upper, length(means)),
#						means, 
#						sd=rep(sqrt(var.scale)*vsf, length(means))), 
#				ncol=nsamp)		
#			
#		Poss.mat[Poss.mat<lower]<-lower
#		Poss.mat[Poss.mat>upper]<-upper
#		
#		Poss.mat
#								
#	})
#
#	Poss.comb<-lapply(1:length(Poss.l), function(ii){
#				
#				first <- t(colSums(Poss.l[[ii]] * (G %*% Poss.l[[ii]]))) + lambda * t(colSums(Poss.l[[ii]] * (1-Poss.l[[ii]])))
#				
#				second <- -2*t(Poss.l[[ii]])%*%W[,ii,drop=FALSE]
#				
#				matrix(t(first)+second, ncol=1)
#			})
#	Poss.comb<-do.call("cbind", Poss.comb)
#	
#	minix<-apply(Poss.comb, 2, which.min)
#	
#	Tnew<-lapply(1:length(minix), function(mi) Poss.l[[mi]][,minix[mi],drop=FALSE])
#	
#	Tnew<-t(do.call("cbind", Tnew))
#	
#	return(Tnew)
#	
#}
#######################################################################################################################

updateT_multicore<-function(method="gini", G, W, Tk, lambdaT, tol, f0, lower=0, upper=1, ncores=2){
	
	if(method=="gini"){
		update.func<-updateT_gini
	}
	
	row.partition<-lapply(0:(ncores-1), function(chunck) (1:(ncol(W)%/%ncores))+chunck*(ncol(W)%/%ncores))
	
	opt_res<-mclapply(row.partition, function(row.ind) update.func(G, W[,row.ind], Tk[,row.ind], lambdaT, tol, f0, lower=0, upper=1), mc.cores=ncores)	
	#opt_res<-lapply(row.partition, function(row.ind) update.func(G, W[,row.ind], Tk[,row.ind], lambdaT, tol, f0, lower=0, upper=1))	
	
	Tk<-do.call("cbind", lapply(opt_res, el, 1))
	fk<-sum(sapply(opt_res, el, 2))
	avg_iter<-mean(sapply(opt_res, el, 3))
	
	return(list(Tk,fk,avg_iter))
}

#######################################################################################################################
#	 Wrappers for the T update methods 	 
#######################################################################################################################

getT.intfac<-function(D, A, V, lambda=0){
	
	G <- A %*% t(A); W <- A %*% t(D);
	Tupd<-updateT_integer(G, W, TT, combin(V,nrow(A)), lambda);
	Tupd
	
}

#######################################################################################################################

getT.empirical<-function(D, A, V, emp.dim=1000, lambda=0){
	
	G <- A %*% t(A); W <- A %*% t(D);
	
	Poss<-t(sapply(1:nrow(A), function(kk){
						
						sample(V,emp.dim, replace=TRUE)
						
					}))
	
	Tupd<-updateT_empirical(G, W, Poss, lambda);
	
	Tupd
	
}

#######################################################################################################################

#getT.resample<-function(D, A, TT, tol=1e-8, maxit=100, verbose=FALSE, ...){
#	
#	err<-norm(D - TT%*%A,'F')/ncol(D)/nrow(D)
#	nit<-0
#	while(err>tol && nit<maxit){
#		
#		if(verbose){
#			cat(sprintf("Number of iterations %g; Error %f\n", nit, err))
#		}
#		
#		nit<-nit+1
#		
#		Tnew<-updateT_resample(D, A, TT, ...)
#			
#		err<-norm(D - Tnew%*%A,'F')/ncol(D)/nrow(D)
#
#	}
#	
#	Tnew
#}

#######################################################################################################################

getT.hlasso<-function(D, A, lambda=0, T0=matrix(runif(nrow(D)*nrow(A)),nrow=nrow(D))){
	
	G <- A %*% t(A); W = A %*% t(D);
	mhcl<-RHLasso(G,W,t(T0),lambda)
	Tnew <- t(mhcl[[1]])
	
}

#######################################################################################################################

getT.gini<-function(D, A, lambda=0, lower=0, upper=1, T0=matrix(runif(nrow(D)*nrow(A)),nrow=nrow(D)),tol=1E-8){
	
	G <- A %*% t(A); W = A %*% t(D)
	f0 <- norm(D - T0 %*% A,'F')^2 + lambda * sum(sum(T0*(1-T0))) - norm(D,'F')^2
	res<-updateT_gini(G, W, t(T0), lambda, tol, f0, lower=lower, upper=upper)
		
	return(t(res[[1]]))
	
}
#######################################################################################################################

tof<-function(tr, G, wcol, lambda){
	
	#obj<-t(tr)%*%G%*%tr - 2 * t(wcol) %*% tr + lambda * t(tr) %*% (1-tr)
	obj<-t(tr)%*%G%*%tr - 2 * t(wcol) %*% tr + lambda * t(tr) %*% (1-tr)
	return(as.numeric(obj))
	
}


#######################################################################################################################

tof.grad<-function(tr, G, wcol, lambda){
	
	#obj<-t(tr)%*%G%*%tr - 2 * t(wcol) %*% tr + lambda * t(tr) %*% (1-tr)
	obj<-2 * G%*%tr - 2 * wcol - 2 * lambda * tr - lambda
	return(as.numeric(obj))
	
}
#######################################################################################################################
#
# Update T using R optimization tools
#
getT.optim<-function(D, A,T0=matrix(runif(nrow(D)*nrow(A)),nrow=nrow(D)), labmda=0){
	
	G <- A %*% t(A); W = A %*% t(D);
	#mhcl<-RHLasso(G,W,t(T0),lambda)
	Tnew<-T0
	for(i in 1:ncol(W)){
		result<-optim(par=T0[i,], fn=tof, gr=tof.grad, G, W[,i], lambda, 
				method="L-BFGS-B",
				lower=rep(0,ncol(Tnew)), upper=rep(1, ncol(Tnew)))
		Tnew[i,]<-result$par
	}
	Tnew
	
}

#######################################################################################################################
# Implemenation of the alternating optimization scheme
#######################################################################################################################
#
# onerun.alternate
#
# Worker routine for alternating optimization-based algorithms 
#
# original MATLAB code by Martin Slawski
#
# R port by Pavlo Lutsik
#
onerun.alternate<-function(
		D, 
		T0, 
		A0,
		Tfix=NULL,
		Tpartial=NULL,
		Tpartial.rows=NULL,
		Apartial=NULL,
		Apartial.cols=NULL,
		sample_var=NULL,
		t.method="quadPen", 
		lambda = 0,
		t.Poss = NULL, 
		normD = NULL, 
		itermax=100,
		qp.rangeT=c(0,1),
		#qp.Alower=rep(0,ncol(T0)),
		#qp.Aupper=rep(1,ncol(T0)),
		qp.Alower=NULL,
		qp.Aupper=NULL,
		emp.dim=500, 
		emp.resample=TRUE,
		emp.vsf=1,
		emp.borders=c(0,1),
		trace=FALSE,
		eps=1e-8,
		blocks=NULL,
		na.values=FALSE,
		verbose=TRUE,
		ncores=1){
	
	k<-ncol(T0)
	if(!is.null(Tfix)){
		kfix <- ncol(Tfix)
		fixix <- (k+1):(k+kfix);
		#qp.Alower <- c(qp.Alower, rep(0, kfix))
		#qp.Aupper <- c(qp.Aupper, rep(1, kfix))
	}else{
		kfix <- 0
	}
	
	
	if(!is.null(Tpartial)){
		if(is.null(Tpartial.rows)){
			Tpartial.rows<-1:nrow(Tpartial)			
		}
		target_rows<-setdiff(1:nrow(D), Tpartial.rows)
	}else{
		target_rows<-1:nrow(D)
	}
	
	if(!is.null(Apartial)){
		if(is.null(Apartial.cols)){
			Apartial.cols<-1:ncol(Apartial)			
		}
		target_cols<-setdiff(1:ncol(D), Apartial.cols)
	}else{
		target_cols<-1:ncol(D)
	}
	
	TT <- T0;
	
	if(!is.null(Tpartial)){
		TT[Tpartial.rows,]<-Tpartial
	}
	
	A <- A0;
	
	if(!is.null(Apartial)){
		A[,Apartial.cols]<-Apartial
	}
		
	
	A0 <-rbind(A0, matrix(0.0, nrow=kfix, ncol=ncol(A0)))
	
	## initialization
	if(is.null(normD)){
		if(!na.values){
			norm_val <- norm(D,'F')^2
		}else{
			norm_val<-sum(D[!is.na(D)]^2)
		}
	}else{
		if(is.null(Tfix)){
			norm_val <- normD
		}else{
			if(!na.values){
				norm_val <- norm(D, '2')
			}else{
				norm_val<-sum(D[!is.na(D)]^2)
			}
		}
	}
	
	if(!na.values){
		diff_norm<-norm(D - TT %*% A,'F')^2
	}else{
		diff<-D - TT %*% A
		diff_norm<-sum(diff[!is.na(diff)]^2)
	}
	
	if(t.method %in% c("Hlasso")){	
		
		f <- diff_norm + lambda * sum(sum(TT))
		
	}else if(t.method %in% c("quadPen", "empirical", "resample")){
		
		f <- diff_norm + lambda * sum(sum(TT*(1-TT)))
		
	}else{
		
		f <- diff_norm
	}
	
	if(verbose){
		cat(sprintf('** Initialization - Objective: %f\n', f));
	}
	
	## consider NA values
	
	if(na.values){
		nas.D<-is.na(D)
		nna.rows<-which(!apply(nas.D,1,any))
		D[nas.D]<-(TT %*% A)[nas.D]
	}
	
	
	## Alternating Optimization Scheme

	iter <- 0 
	Conv<-f
	Dt<-t(D)
	
	if(trace){
		As<-list()
		Ts<-list()
	}
	
	while(TRUE) {
		
		iter <- iter + 1;
		if(verbose){
			cat(sprintf('** Iter %d **\n', iter)) #iter starts from 2
		}
		
		####################################### UPDATE T #############################
		#
		# min_TT' norm(Dt - At * TT')^2 + lambda PSI(TT')
		# subject to lower >= each comp. of TT' >= upper
		#
		##############################################################################
		
		## if the missing value variant, start with updating A
		#if(iter>1 || !na.values){
		
			if(verbose){
				cat(sprintf("[Alternating run:] Updating T...\n"))
			}
			if(t.method == "integer"){
				
				if(is.null(t.Poss)){
					stop("A vector of possible values should be supplied for this method")
				}
				G <- A %*% t(A); W <- A %*% Dt;
				Tnew <- updateT_integer(G, W, TT, t.Poss, lambda);
				#ftemp = ftemp + norm_val;
			
			}else if(t.method =="empirical"){
							
				if(is.null(t.Poss)){
					stop("A vector of possible values should be supplied for this method")
				}
				
				if(emp.resample){
					
					V<-t.Poss
					t.Poss<-t(sapply(1:ncol(T0), function(kk){
						
							sample(V, emp.dim, replace=TRUE)
												
						}))
				}
				
				G <- A %*% t(A); W <- A %*% Dt;
				Tnew <- updateT_empirical(G, W, t.Poss, lambda);
				#ftemp = ftemp + norm_val;
			
			}else if(t.method=="Hlasso"){
				
				G <- A %*% t(A); W <- A %*% Dt;
				mhcl<-RHLasso(G,W,t(TT),lambda)
				Tnew <- t(mhcl[[1]]); ftemp <- mhcl[[2]]
				ftemp <- ftemp + norm_val;
				
				# if the regularization is too strong
				if(sum(Tnew!=0) == 0){
					if(verbose){
						cat(sprintf('\n [Alternating run:] T becomes zero matrix. Exit.\n'))
					}
					TT <- Tnew;
					f <- ftemp;
					break;
				}
				
			}else if(t.method=="optim"){
										
				Tnew<-getT.optim(D, A, TT, lambda)
				f<-norm(D - Tnew %*% A,'F')^2 #+ lambda * sum(sum(TT))
							
			}else if(t.method=="quadPen"){
				
				if(is.null(Tpartial)){
					Dt4T<-Dt
					Tstart<-t(TT)
				}else{
					Dt4T<-Dt[,target_rows]
					Tstart<-t(TT[target_rows,])
				}
							
				G <- A %*% t(A); W <- A %*% Dt4T
				##%%%mexHCLasso(G,W,T',lambda);
				
				if(ncores>1){
					res <- updateT_multicore("gini", G, W, Tstart, lambda, eps, f - norm_val, lower=qp.rangeT[1], upper=qp.rangeT[2], ncores=ncores);
				}else{
					res <- updateT_gini(G, W, Tstart, lambda, eps, f - norm_val, lower=qp.rangeT[1], upper=qp.rangeT[2]);
		 		}
				Trecov <- t(res[[1]]); ftemp<-res[[2]]
				
				if(!is.null(Tfix)){
					Tnew<-cbind(Trecov, Tfix)
				}else{
					Tnew<-Trecov
				}
				
				if(!is.null(Tpartial)){
					Ttmp<-matrix(NA_real_, nrow=nrow(T0), ncol=ncol(T0))
					Ttmp[target_rows,]<-Tnew
					Ttmp[Tpartial.rows,]<-Tpartial
					Tnew<-Ttmp
				}
								
			}else if(t.method=="resample"){
				
				Tnew<-updateT_resample(D, A, TT, lambda, lower=0, upper=1, vsf=1, nsamp=emp.dim)
							
			}
				
#		}else{
#			Trecov<-TT
#			if(!is.null(Tfix)){
#				Tnew<-cbind(Trecov, Tfix)
#			}else{
#				Tnew<-Trecov
#			}
#		}
		
		###################################### UPDATE A ##############################
		#
		# min_A norm(D - T * A)^2
		# subject to each clm of A on the simplex
		#
		##############################################################################
		
		if(verbose){
			cat(sprintf('[Alternating run:] Updating A...\n'))
		}
		
		if(is.null(blocks)){
			
			if(is.null(Apartial)){
				D4a<-D
			}else{
				D4a<-D[,target_cols]
			}
			
			if(!na.values){
				G <- t(Tnew) %*% Tnew; W <- t(Tnew) %*% D4a
			}else{
				G <- t(Tnew[nna.rows,]) %*% Tnew[nna.rows,]; W <- t(Tnew[nna.rows,]) %*% D4a[nna.rows,]
			}
			
			if(!is.null(qp.Alower) && !is.null(qp.Aupper)){
				mqs<-RQuadSimplexBox(G, W, A0, qp.Alower, qp.Aupper, eps)
			}else{
				mqs<-RQuadSimplex(G, W, A0, eps)
			}
			Anew <- mqs[[1]]; ftemp<-mqs[[2]]; A0<-Anew
			
			if(t.method %in% c("Hlasso", "integer")){
				fnew <- ftemp + normD + lambda * sum(sum(Trecov));
			}else if(t.method %in% c("quadPen", "empirical", "resample")){
				fnew <- ftemp + normD + lambda * sum(Trecov * (1 - Trecov));
			}else{
				fnew <- ftemp + normD
			}
			
		}else{
			
			Anew <- matrix(0, nrow(A), ncol(A));
			
			G0 <- t(Tnew[,c(blocks$c, blocks$s0)]) %*% Tnew[,c(blocks$c, blocks$s0)]; 
			W0 <- t(Tnew[,c(blocks$c, blocks$s0)]) %*% D[,blocks$pheno0]
			
			mqs<-RQuadSimplexBox(G0,W0,A[c(blocks$c,blocks$s0), blocks$pheno0], 
					qp.Alower[c(blocks$c,blocks$s0)], qp.Aupper[c(blocks$c,blocks$s0)], eps);
			Anew[c(blocks$c, blocks$s0),blocks$pheno0]<-mqs[[1]]; ftemp0<-mqs[[2]]
			
			G1 <- t(Tnew[,c(blocks$c, blocks$s1)]) %*% Tnew[,c(blocks$c, blocks$s1)]; 
			W1 <- t(Tnew[,c(blocks$c, blocks$s1)]) %*% D[,blocks$pheno1]
			
			mqs<-RQuadSimplexBox(G1,W1,A[c(blocks$c,blocks$s1), blocks$pheno1], 
					qp.Alower[c(blocks$c,blocks$s1)], qp.Aupper[c(blocks$c,blocks$s1)],  eps);
			Anew[c(blocks$c, blocks$s1),blocks$pheno1]<-mqs[[1]]; ftemp1<-mqs[[2]]
						
			if(t.method %in% c("Hlasso", "integer")){
				fnew <- ftemp0 + ftemp1 + normD + lambda * sum(sum(Trecov));
			}else if(t.method %in% c("quadPen", "empirical", "resample")){
				fnew <- ftemp0 + ftemp1 + normD + lambda * sum(Trecov * (1 - Trecov));
			}else{
				fnew <- ftemp0 + ftemp1 + normD
			}
		}
		
		if(!is.null(Apartial)){
			Atmp<-matrix(NA_real_, nrow=nrow(A0), ncol=ncol(A0))
			Atmp[,target_cols]<-Anew
			Atmp[,Apartial.cols]<-Apartial
			Anew<-Atmp
		}
		
		if(na.values){
			D[nas.D]<-(Tnew %*% Anew)[nas.D]
			Dt<-t(D)
		}
		
		
		# relative reduction, must be absolute value
		red_f <- (f - fnew) / f 
		
		if(verbose){
			cat(sprintf("Objective: %f, F reduc: %g\n", fnew, red_f))
		}
		
		red_f<-abs(red_f)
		
		#TT <- Tnew;
		TT <- Trecov;  
		if(!is.null(Tfix)){
			A <- Anew[1:k, ,drop=FALSE]
			Afix <- Anew[fixix, ,drop=FALSE]
			Dt <- t(D) - t(Afix) %*% t(Tfix)
			norm_val <-  sum(Dt^2)
		}else{
			A <- Anew
		}
		f <- fnew; Conv<-c(Conv,f)
		
		if(trace){
			Ts[[iter]]<-TT; As[[iter]]<-A
		}
		
		if(red_f < eps){
			if(verbose){
				cat(sprintf('[Alternating run:] alternating updates: objective value converges.\n'))
			}
			break
		}
		
		if(iter >= itermax){
			if(verbose){
				cat(sprintf('[Alternating run:] alternating updates: reach iteration limits.\n'))
			}
			break
		}
	}
	
	Fval<-f
	#Conv<-Conv[1:iter]
	if(trace){
		As<-As[2:length(As)]
		Ts<-Ts[2:length(Ts)]
		result<-list("T"=TT, "A"=A, "Fval"=Fval, "Conv"=Conv, "Ts"=Ts, "As"=As)
		if(!is.null(Tfix)){
			result$Afix<-Afix
		}
		return(result)
	}
	
	rmse<-sqrt(sum((D-TT%*%A)^2)/ncol(D)/nrow(D))
	result<-list("T" = TT, "A" = A, "Fval" = Fval, "Conv" = Conv, "rmse" = rmse)
	if(!is.null(Tfix)){
		result$Afix<-Afix
	}
	return(result)
}



#######################################################################################################################
# Implemenation of the alternating optimization scheme
#######################################################################################################################
#
# a wrapper for cppTAfact
#
# cppTAfact - alternating optimization framework to solve the following
# problem:
#		find T, A such that 0.5 * (||D - TA||_F)^2 + lambda ,
# where D is a mxn matrix with entries between 0 and 1,
# T is a mxr matrix with entries between 0 and 1,
# A is a rxn matrix with nonnegative values and columns summing up
# to 1.
#
# author: Nikita Vedeneev
#
# R port by Pavlo Lutsik
#

onerun.cppTAfact<-function(
		D, 
		T0, 
		A0,
		Tfix=NULL,
		Tpartial=NULL,
		Tpartial.rows=NULL,
		Apartial=NULL,
		Apartial.cols=NULL,
		sample_var=NULL,
		t.method="quadPen", 
		lambda = 0,
		t.Poss = NULL, 
		normD = NULL, 
		itermax=100,
		qp.rangeT=c(0,1),
		#qp.Alower=rep(0,ncol(T0)),
		#qp.Aupper=rep(1,ncol(T0)),
		qp.Alower=NULL,
		qp.Aupper=NULL,
		emp.dim=500, 
		emp.resample=TRUE,
		emp.vsf=1,
		emp.borders=c(0,1),
		trace=FALSE,
		eps=1e-8,
		blocks=NULL,
		na.values=FALSE,
		verbose=TRUE,
		ncores=1){

	res<-cppTAfact(
			t(D), #- a transposed D matrix,
			t(T0), #-Ttinit - a transposed init for T matrix,
			A0, # - an initial value for A matrix,
			lambda,# - regularizer parameter (0.0 by default),
			itermax, #itersMax, - a max number of alternations (1000 by default),
			eps, #tol - tolerance for alternations (1e-8 by default),
			10*eps, #tolA - tolerance for opt wrt A (1e-7 by default),
			10*eps #tolT - tolerance for opt wrt T (1e-7 by default)
	)
	### TODO: modify cppTAfact to output the list is identical to the output of onerun.alternate
	#
	#cppTAfact returns a named list where:
    #res$Tt - a transposed estimated of T matrix,
    #res$A - an estimate of A matrix,
    #res$niter - a total number of alternations
    #res$objF - objective value at res$Tt and res$A
	#
	result<-list("T" = t(res$Tt), "A" = res$A, "Fval" = res$objF, "Conv" = res$niter, "rmse"= res$rmse)
	return(result)
}

#######################################################################################################################
#'
#' factorize.alternate
#' 
#' Matrix factorization algorithms based on the alternating optimization scheme
#' 
#' @param D 			m by n input matrix with mixture data 
#' @param k 			number of latent components, \code{integer}
#' @param t.method 		method for updating the latent component matrix, one of 
#' 						\code{"integer", "empirical", "Hlasso"} or \code{"quadPen"}
#' @param Tfix			an optional matrix of a priori known fixed components
#' @param V				for \code{t.method} "integer" a small vector of possible values;
#' 						for \code{t.method} "empirical" a vector giving empirical distribution
#' 						of T values
#' 
#' @param init			type of initialization, either "random" (default) or "fixed".
#' 
#' @param opt			if \code{init} is "random" 
#' 						number of runs with independent initialization,
#' 						if \code{init} is "fixed"
#' 						starting values for T and A (see details)
#' 
#' @param emp.dim 		for \code{t.method} "empirical", 
#' 						number of randomly drawn samples for T row selection
#'  		 
#' @param emp.resample  for \code{t.method} "empirical",
#' 						a flag indicating whether resampling should
#' 						be done at each iteration
#' 
#' @param itermax		maximal number of iterations
#'
#' @param trace			a flag indicating whether to return the 
#' 						factorization results for each iteration 
#' 
#' @param eps			threshold for objective value change
#' 
#' @param ncores		number of CPU cores used for parallelization
#' 
#' @param pheno			a list with phenotypic information
#' 
#' @param verbose		flag specifying whether to show diagnostic
#' 						statements during the execution
#' 						
#' @details				In case \code{init} is "fixed" the starting values
#' 						for the m by k matrix of latent components 
#' 						and for the k by n matrix of mixing proportions 
#' 						should be specified as T and A elements of a list
#' 						supplied as \code{opt}
#' 
#' @return a \code{list} with the following elements:
#' 			\describe{
#'		 				\item{\code{T}}{matrix of latent components}
#'   					\item{\code{A}}{matrix of mixture proportions}
#' 						\item{\code{Fval}}{the final value of the objective function}
#' 						\item{\code{Conv}}{sequence of objective function values
#' 								attained after each iteration} 
#' 						\item{\code{rmse}}{RMSE of the factorization}
#' 			}
#' 
#' 
#' @author Martin Slawski
#' @author R port by Pavlo Lutsik
#' 
#' @export
#' 
factorize.alternate<-function(D, 
		k,
		method="MeDeCom.quadPen",
		t.method="quadPen",
		Tfix=NULL,
		Tpartial=NULL,
		Tpartial.rows=NULL,
		Apartial=NULL,
		Apartial.cols=NULL,
		sample_var=NULL,
		V=NULL,
		lambda = 0,
		init="random", 
		opt=5, 
		emp.dim=500,
		emp.resample=TRUE,
		emp.vsf=1,
		emp.borders=c(0,1),
		qp.rangeT=c(0,1),
		#qp.Alower=if(is.null(Tfix)) rep(0,k) else rep(0,k+ncol(Tfix)),
		#qp.Aupper=if(is.null(Tfix)) rep(1,k) else rep(1,k+ncol(Tfix)),
		qp.Alower=NULL,
		qp.Aupper=NULL,
		itermax=100, 
		trace=FALSE,
		eps=1e-8,
		ncores=1,
		pheno=NULL,
		na.values=FALSE,
		seed=NULL,
		verbosity=0L){
	
	if(!t.method %in% c("integer", "empirical", "resample", "Hlasso", "optim", "quadPen", "cppTAfact")){
		stop("supplied optimization method for T is not implemented")
	}
	
	n<-nrow(D);
	d<-ncol(D);
	
	if(!is.null(Tfix)){
		if(nrow(Tfix)!=n){
			stop("The supplied Tfix is not corresponding to the data matrix D")
		}
		fixedT<-TRUE
	}else{
		fixedT<-FALSE
	}
	
	if(!is.null(pheno)){
		
		d0 <- sum(pheno$labels == 0);
		d1 <- sum(pheno$labels == 1);
		
		if(d0 + d1 != d){ 
			stop('pheno.labels is not a 0-1 vector or number of labels does not agree with the columns of D') 
		}
		
		k0 <- pheno$k0;
		k1 <- pheno$k1;
				
		if(k0 > k1){
			stop('pheno.k0 has to be no larger than pheno.k1') 
		}
		
		if(k1 > 2 * k0){
			error('pheno.k1 should not be larger than 2*pheno.k0')  
		}
		
		c <- pheno$c
		
		if(pheno$c > 2 * k0 - 1){
			stop('pheno.c is too large') 
		}
		
		k <- k0 + k1 - c
		s0 <- k0 - c
		s1 <- k1 - c
		
		blocks<-list()
		
		if(c>0){
			blocks$c <- 1:c
		}else{
			blocks$c<- integer()
		}
		if(c+1<=c+s0){
			blocks$s0 <- (c+1):(c+s0)
		}else{
			blocks$s0 <- integer()
		}
		if(c+s0+1<=c+s0+s1){
			blocks$s1 <- (c+s0+1):(c+s0+s1)
		}else{
			blocks$s1 <- integer()
		}
		
		blocks$pheno0 <- which(pheno$labels == 0)
		blocks$pheno1 <- which(pheno$labels == 1)
		
	}else{
		blocks<-NULL
	}
	
	if(!na.values){
		normD <- norm(D,'F')^2; # temporary variable used in computation
	}else{
		normD <- sum(D[!is.na(D)]^2)
	}
	Poss<-NULL
	
	if(t.method=="integer"){
		
	  	Poss<-combin(V,k)
				
	}else if(t.method=="empirical"){

		if(emp.resample){ ## pass over an initial vector
			Poss<-V
		}else{  		  ## create a matrix of possibilities
			Poss<-t(sapply(1:k, function(kk){
								sample(V, emp.dim, replace=TRUE)
							}))
		}
	}
	
	if(init=="random"){
		numruns <- opt
		if(!is.null(seed)){
			set.seed(seed)
		}
	}else if(init=="fixed"){
		numruns <- 1
	}else{
		stop("No such initialization method")
	}
	
	T0s<-vector("list", numruns)
	A0s<-vector("list", numruns)
	Ts<-vector("list", numruns)
	As<-vector("list", numruns)
	if(fixedT){
		Afixs <- vector("list", numruns)
	}
	Fvals <- vector("numeric", numruns)
	Convs <- vector("list", numruns)
	
	#if(init=="random"){
	

	call_onerun<-function(run){
		
		if(verbosity>1L){
			cat(sprintf('[Alternating procedure:]----------- %d runs -----------\n', numruns));
			cat(sprintf('[Alternating procedure:]----------- Run %d ------------\n', run));
		}
		
		if(init=="random"){
		# generating random starting matrices
			
			if(t.method %in% c("integer")){
				T0 <- round(matrix(runif(n*k),nrow=n)) # enforce 0-1 constraints by rounding
			}else{
				T0 <- matrix(runif(n*k, min=qp.rangeT[1], max=qp.rangeT[2]),nrow=n)
			}
			
			if(is.null(blocks)){
				A0 <- randsplxmat(k,d)
				if(!is.null(qp.Alower) && !is.null(qp.Aupper)){
					for(kk in 1:d){
						A0[,kk] <- RProjSplxBox(A0[,kk,drop=FALSE], qp.Alower, qp.Aupper);
					}
				}
			}else{
				A0<-matrix(0, k, d)
				A0[c(blocks$c,blocks$s0), blocks$pheno0] <- randsplxmat(k0, d0);
				A0[c(blocks$c,blocks$s1), blocks$pheno1] <- randsplxmat(k1, d1);
			}
			
			if(trace){
				T0s[[run]]<<-T0
				A0s[[run]]<<-A0
			}
		
		}else if(init=="fixed"){
			T0 <- opt$T;  A0 <- opt$A;
			if(!is.null(qp.Alower) && !is.null(qp.Aupper)){
				for(kk in 1:d){
					A0[,kk] <- RProjSplxBox(A0[,kk,drop=FALSE], qp.Alower, qp.Aupper);
				}
			}
		}
		
		# solve the topic model
		if(method == "MeDeCom.cppTAfact"){
			onerun.function<-onerun.cppTAfact
		}else{
			onerun.function<-onerun.alternate
		}

		onerun <- onerun.function(
				D, T0, A0,
				Tfix = Tfix,
				Tpartial=NULL,
				Tpartial.rows=NULL,
				Apartial=NULL,
				Apartial.cols=NULL,
				sample_var=sample_var,
				t.method = t.method, 
				t.Poss = Poss, 
				normD=normD, 
				lambda=lambda,
				itermax=itermax, 
				qp.rangeT=qp.rangeT,
				qp.Alower=qp.Alower,
				qp.Aupper=qp.Aupper,
				emp.dim=emp.dim, 
				emp.resample=emp.resample, 
				emp.vsf=emp.vsf, 
				emp.borders=emp.borders,
				trace=trace,
				eps=eps,
				blocks=blocks,
				na.values=na.values,
				verbose=verbosity>1L,
				ncores=ncores);
		
#			Ts[[run]]<<-onerun[[1]]
#			As[[run]]<<-onerun[[2]]
#			if(fixedT){
#				Afixs[[run]]<<-onerun[[2]]
#			}
#			Fvals[[run]]<<-onerun[[3]]
#			Convs[[run]]<<-onerun[[4]]

		
			return(onerun)
					
		}
		
#		pcoordinates <- foreach(target = targets) %dopar%
#				tryCatch(RnBeads::rnb.execute.dreduction(rnb.set, target = target), error = function(e) { e$message } )
		
		#if(ncores>1){
		if(FALSE){
			#require(doMC)
			#registerDoMC(N_CORES)
			#cl<-makeCluster(N_CORES)
		    #result_list<-foreach(run = 1:numruns) %dopar% call_onerun(run)
			result_list<-mclapply(1:numruns, call_onerun, mc.cores=ncores)
		}else{
			result_list<-lapply(1:numruns, call_onerun)
		}

		dump<-sapply(1:length(result_list), function(run){
			onerun<-result_list[[run]]
			Ts[[run]]<<-onerun[[1]]
			As[[run]]<<-onerun[[2]]
			if(fixedT){
				Afixs[[run]]<<-onerun[[2]]
			}
			Fvals[[run]]<<-onerun[[3]]
			Convs[[run]]<<-onerun[[4]]
		})
		
		# pick the result with the smallest objective value
		if(numruns>1){
			idx <- which.min(Fvals)
		}else{
			idx <- 1L
		}
		
		Fval <- Fvals[[idx]];
		TT <- Ts[[idx]];
		A <- As[[idx]];
		Conv <- Convs[[idx]]
		
		if(!na.values){
			#rmse<-sqrt(sum((D-TT%*%A)^2)/ncol(D)/nrow(D))
			rmse <- Fval - lambda*sum(TT-TT^2)
		}else{
			diff_mat<-D-TT%*%A
			rmse<-sqrt(sum(diff_mat[!is.na(diff_mat)]^2)/ncol(D)/nrow(D))
		}
		
		if(trace){
			result<-list("T"=TT, "A"=A, "Fval"=Fval, "Conv"=Conv, "rmse"=rmse, "Ts"=Ts, "As"=As, "Fvals"=Fvals, "Convs"=Convs, "T0s"=T0s, "A0s"=A0s)
		}else{
			result<-list("T"=TT, "A"=A, "Fval"=Fval, "Conv"=Conv, "rmse"=rmse)
		}
		if(fixedT){
			result$Afix <- Afixs[[idx]]
			
		}
		return(result)
		
#	}else if(init=="fixed"){
#		
#		T0 <- opt$T;  A0 <- opt$A;
#		if(!is.null(qp.Alower) && !is.null(qp.Aupper)){
#			for(kk in 1:d){
#				A0[,kk] <- RProjSplxBox(A0[,kk,drop=FALSE], qp.Alower, qp.Aupper);
#			}
#		}
#		onerun <- onerun.alternate(
#				D, T0, A0,
#				Tfix = Tfix,
#				Tpartial=NULL,
#				Tpartial.rows=NULL,
#				Apartial=NULL,
#				Apartial.cols=NULL,
#				t.method = t.method, 
#				t.Poss = Poss, 
#				normD=normD, 
#				lambda=lambda, 
#				itermax=itermax,
#				qp.rangeT=qp.rangeT,
#				qp.Alower=qp.Alower,
#				qp.Aupper=qp.Aupper,
#				emp.dim=emp.dim, 
#				emp.resample=emp.resample, 
#				emp.vsf=emp.vsf,
#				emp.borders=emp.borders,
#				trace=trace, 
#				eps=eps, 
#				blocks=blocks,
#				na.values=na.values,
#				verbose=verbosity>1L);
#					
#		return(onerun)
#		
#	}else{
#		
#		stop("No such initialization method")
#	}
	
}
#######################################################################################################################
# Matrix factorization with binary components (from Slawski et al NIPS 2013)
#######################################################################################################################
#
#   factorize.exact
#
#   Exact factorization of the data matrices
#
# -- Input --
#
#   D      - input data matrix
#   V      - finite set of values the components are allowed to take
#   r      - number of components
#   opt    - various options, see Integerfac_findvert_top
#            Default options are used in case that this argument is omitted. 
#
#
# -- Output --
#
#   T       -  matrix of components
#
#   A       -  coefficient matrix
#
#   status  - a struct containing additional status information
#  
#   Original MATLAB code by Martin Slawski and Matthias Hein
#
#	R port by Pavlo Lutsik
#
factorize.exact<-function(D, r, V=c(0,1), opt=NULL, zeroprob = (0 %in% V)){
	
	
	if(is.null(opt)){
		opt <- opt.factorize.exact();# initialize with default options
	}	
	V<-sort(V)
	
	if(opt$affine) {
		zeroprob=FALSE
	}
	
	m=nrow(D); n=ncol(D);
	cardV <- length(V);
	myeps <- 1e-10;
	
	status<-list()
	
	if(opt$affine){
		meanD <- rowMeans(D);
		E <- D - meanD %*% ones(1,n);
		k <- r-1;
	}else{
		meanD <- zeros(m, 1); # to simplify the code
		E <- D;
		k <- r;
	}
	
#
	svds<-svd(E, nu=300);
	
	U<-svds$u; sigma<-svds$d
	
	if(sigma[k] < sqrt(myeps)){
		print('Data matrix does not match the specified number of components r: 
						one of the top singular values is (close) to zero')
	}
	ixx = 1:k; 
	TP<-matrix(NaN,nrow=m,ncol=r); # TP is the matrix for the found topics (vertices)
	All = 1:m;
	
	dT<-vector("numeric", r)
	
	if(opt$aggr=='no'){
		nsamples_basic = 1; nsamples_extra = 0;
	}else{
		nreplace <- round(opt.replace * k);
		nreplace <- max(1, nreplace); # at least one coordinate has to be replaced
		nsamples_basic <- floor((m - k)/nreplace);
		
		if(opt$nsamples == 0){
			nsamples_extra <- 0;
		}else{
			if(opt$nsamples > 0){
				nsamples_tot <- opt$nsamples;
				if(nsamples_tot > nsamples_basic){
					nsamples_extra = nsamples_tot - nsamples_basic;
				}else{
					nsamples_basic = min(opt$nsamples, nsamples_basic); nsamples_extra = 0;
				}
			}else{
				nsamples_extra = -(opt$nsamples);
			}
		}
	}
	nsamples = nsamples_basic + nsamples_extra;
	
	### some preparations outside the loop
	if(k>2){
		rest <- max(k-opt$chunksize,1);
		PossOver_0 <- combin(V, rest);
		Poss_0 <- combin(V, k-rest);
	}else{
		Poss_0 <- combin(V, k);
	}
	
	if(opt$nonnegative){
		if(opt$affine){
			A0 <- randsplxmat(r,n);
		}else{
			A0 <- rand(r,n);
		}
	}
	
	singleerr<-vector("numeric", nsamples)
	###
	
	for(i in 1:nsamples){
		if(i <= nsamples_basic){
			kxx <- licols(t(U[All,ixx,drop=F]))[["idx"]]
			
			jxx<-All[kxx];
			if(!opt$aggr=="no"){ 
				All<-All[-kxx[1:nreplace]]
			}
			
			if(opt$verbose) {
				cat(sprintf('Round: %d - Selecting: %s\n', i, paste(jxx,collapse=",") ));
			}else{ # choose extra coordinates by random selection
				rpm <- randperm(m);
				jxx <- rpm[1:k];
			}
			
			FD<-U[jxx,ixx]; 
			v <- mldivide(FD, eye(k));
			TT <- U[,ixx,drop=F]%*%v;
			Indices <- 1:m;  # we only have to test for the rows (Indices) which were not selected above (jxx)
			Indices <- Indices[-jxx];
			
			if(k>2){
				
				PossOver <- PossOver_0 - repmat(matrix(meanD[jxx[1:rest]], nrow=rest), 1, cardV^rest);
				Poss <-  Poss_0 - repmat(matrix(meanD[jxx[(rest+1):k]], nrow=length((rest+1):k)), 1, cardV^(k-rest));
				
				num <- nrow(PossOver); nump<-num+1;
				Poss <- rbind(Poss,ones(1,ncol(Poss))); # we append a one to the vector to integrate the fixed part
				
				TT1<-TT[Indices,1:num,drop=F]; 
				TT2<-TT[Indices,nump:ncol(TT),drop=F];   
				# thus we also select just the rows Indices from T (the rest is the identity matrix) and
				# divide it up into the fixed component (rest) and the chunk
				
				for(kk in 1:ncol(PossOver)){
					
					vec<-PossOver[,kk,drop=F];
					TPP1 <- TT1%*%vec; # these are the first components
					TPP <- cbind(TT2,TPP1)%*%Poss + repmat(matrix(meanD[Indices],nrow=length(Indices)),1,ncol(Poss));   
					# these are the first offset rows of the potential topics
					TTT <- abs(TPP - projectV(TPP, V));
					if(zeroprob && kk==1){
						TTT <- TTT[,2:ncol(TTT),drop=F]; # get rid of the all-zero-vector
					}
					sres<-sort(colSums(TTT^2),index.return=T);
					
					sdist<-sres$x; icc<-sres$ix;
					
					icc <- icc + (zeroprob && kk==1);
					if(kk==1){
						X1<-TPP[,icc[1:r],drop=F];
						X2<-rbind(repmat(vec+meanD[jxx[1:rest]],1,r),Poss[1:(nrow(Poss)-1),icc[1:r]]+
										repmat(matrix(meanD[jxx[(rest+1):k]],nrow=length((rest+1):k)),1,r));
						if(i==1){
							TP<-zeros(nrow(D),ncol(X1));
						}
						TP[Indices,((i-1)*r+1):(i*r)]<-X1;
						TP[jxx,((i-1)*r+1):(i*r)]<-X2;
						dT[1:r]<-sdist[1:r];
						tjx<-icc[1:r];
					}else{
						mixdT<-c(dT,sdist[1:r]); # join new and old distances
						sres<-sort(mixdT,index.return=T); smixdT<-sres$x; idd<-sres$ix
						thresh<-smixdT[r]; # this is the new maximal distance
						
						freeslots <- which(dT>thresh);
						tjx <- icc[1:r][sdist[1:r]<=thresh];
						if(length(tjx)>0){
							for(ii in 1:length(tjx)){
								
								X1<-TPP[,tjx[ii],drop=F];
								X2<-rbind(vec+meanD[jxx[1:rest]],Poss[1:(nrow(Poss)-1),tjx[ii],drop=F]+meanD[jxx[(rest+1):k]]);
								TP[Indices,(i-1)*r+freeslots[ii]]<-X1
								TP[jxx,(i-1)*r+freeslots[ii]]<-X2;
							}
							dT[freeslots]<-sdist[1:length(tjx)];
						}
					}
				}
				
			}else{
				Poss <- Poss_0 - repmat(matrix(meanD[jxx[1:k]],nrow=k),1,cardV^k);
				TPP <- TT%*%Poss + repmat(matrix(meanD, nrow=length(meanD)), 1, cardV^k);
				TTT<-abs(TPP - projectV(TPP, V));
				sres<-sort(colSums(TTT^2),index.return=T); icc<-sres$ix
				TP[,((i-1)*r+1):(i*r)]<-TPP[,icc[(1+zeroprob):(r+zeroprob)], drop=F];
				
			}
			
			if(opt$verbose) cat(sprintf('Single Fit - Iteration: %d\n',i))
			
			TP <- projectV(TP, V);
			
			T2 <- TP[,((i-1)*r+1):(i*r)];
			if(opt$nonnegative){
				G <- t(T2) %*% T2; W = t(T2) %*% D;
				if(opt$affine){
					mqr <- RQuadSimplex(G,W,A0, myeps);
					A2<-mqr[["A"]]
				}else{
					mqnn<-RHLasso(G,W,A0, myeps);
					A2<-mqnn[["A"]]
				}
			}else{
				if(opt$affine){
					A2 <- affls(T2, D);
				}else{
					A2 <- mldivide(T2,D); # least squares      
				}
			}
			singleerr[i]= sum(sum(abs(D-T2%*%A2)^2));
		}
	}	
	sres<-sort(singleerr,index.return=T);bestsingleerr<-sres$x;mix<-sres$ix
	bestsingleerr<-bestsingleerr[1];
	if(opt$verbose) cat(sprintf('Best single error: %f\n',bestsingleerr))
	
	if(opt$aggr %in% c('no', 'bestsingle')){
		T <- TP[,((mix[1]-1)*r+1):(mix[1]*r)]; 
		if(opt$nonnegative){
			G = t(T) %*% T; W = t(T) %*% D;
			if(opt$affine){
				mqr<- RQuadSimplex(G,W,A0,myeps); 
				A<-mqr[["A"]]
			}else{
				mqnn<-RHLasso(G,W,A0,myeps);
				A<-mqnn[["A"]]
			}	
		}else{
			if(opt$affine){
				A <- affls(T, D);
			}else{
				A <-mldivide(T,D);
			}
		}
		
		rmse<-sqrt(sum((D-T%*%A)^2)/ncol(D)/nrow(D))
		
		status$err <- sum(sum(abs(D-T%*%A)^2)); status$nsamples <- nsamples;
		
		return(list(T=T,A=A,status=status,rmse=rmse))
	}
	
	# merge the best candidate sets
	for(i in 1:min(nsamples, opt.naggsets)){
		TP2[,((i-1)*r+1):(i*r)]<-TP[,((mix(i)-1)*r+1):(mix(i)*r)];
	}
	TP<-TP2;
	
	TC[,1]<-TP[,1];
	counter<-2;
	for(i in 2:ncol(TP)){
		if(dist(t(TC),t(TP[,i]),method="euclidean")>1E-3*m){
			TC[,counter]<-TP[,i];
			counter<-counter+1;
		}
	}
	rm(TP);
	if(opt$verbose) cat(sprintf('Number of topics after duplicate checking: %d\n',ncol(TC)))
	
	#### TO BE MODIFIED
	
	if(opt$nonnegative){
		if(opt$affine){
			T2<-TC; G <- t(T2) %*% T2; A0 <- randsplxmat(ncol(T2),n); W <- t(T2) %*% D;
			mqr<- RQuadSimplex(G,W,A0, myeps); A2<-mqr[["A"]]
			err<-sum(sum(abs(D-T2%*%A2)^2));
			
			if(opt$verbose) cat(sprintf('With %dm topics - err: %f\n',ncol(TC),err))
			
			while(ncol(TC)>r){
				rm(err);
				for(i in 1:ncol(TC)){
					itt<-setdiff(1:ncol(TC),i);
					T2=TC[,itt,drop=F];
					if(rem[i,10]==0){
						if(opt$verbose) cat(sprintf('To do: %d - done: %d\n', ncol(TC),i))
						G <- t(T2) %*% T2; A0 <- randsplxmat(ncol(T2),n); W = t(T2) %*% D;
						mqr<- RQuadSimplex(G,W,A0, myeps); A2<-mqr[["A"]]
						err[i]<-sum(sum(abs(D-T2%*%A2)^2));
					}
				}
				if(ncol(TC)<4*r){
					disc<-which.min(err);
					if(opt$verbose){ cat(sprintf('Current size: %d - Minimal error: %f - discarding: %d\n',
										ncol(TC),min(err),disc))}
					itt<-setdiff(1:ncol(TC),disc);
					TC<-TC[,itt];
				}else{
					dec<-5;
					if(ncol(TC)>k+50){ dec=20; }
					if(ncol(TC)>k+100){ dec=50; }
					if(ncol(TC)>k+200){ dec=100; }
					if(ncol(TC)>k+2000){ dec=1000; }
					sres<-sort(err,index.return=T); serr<-sres$x; iuu<-sres$ix
					if(opt$verbose){ cat(sprintf('Minimal error: %f  Error %d : %f - discarding: %d\n',
										serr[1],dec,serr[dec], iuu[1:dec]))}
					itt<-setdiff(1:ncol(TC),iuu[1:dec]);
					TC<-TC[,itt,drop=F];
				}
			}
			
			# final fit
			T2<-TC; G <- t(T2) %*% T2; A0 <- randsplxmat(ncol(T2),n); W <- t(T2) %*% D;
			mqr<- RQuadSimplex(G,W,A0, myeps);  A2<-mqr[["A"]]
			
		}else{
			
			T2<-TC; G <- t(T2) %*% T2; A0 <- matrix(runif(ncol(T2)*n), ncol(T2), n); W <- t(T2) %*% D;
			mqnr<- mexQuadNN(G,W,A0, myeps); A2<-mqnr[["A"]]
			
			err<-sum(sum(abs(D-T2%*%A2)^2));
			
			if(opt$verbose) cat(sprintf('With %dm topics - err: %f\n',ncol(TC),err))
			
			while(ncol(TC)>r){
				rm(err);
				for(i in 1:ncol(TC)){
					itt<-setdiff(1:ncol(TC),i);
					T2=TC[,itt,drop=F];
					if(rem[i,10]==0){
						if(opt$verbose) cat(sprintf('To do: %d - done: %d\n', ncol(TC),i))
						G <- t(T2) %*% T2; A0 <- randsplxmat(ncol(T2),n); W = t(T2) %*% D;
						mqr<- RQuadSimplex(G,W,A0, myeps); A2<-mqr[["A"]]
						err[i]<-sum(sum(abs(D-T2%*%A2)^2));
					}
				}
				if(ncol(TC)<4*r){
					disc<-which.min(err);
					if(opt$verbose){ cat(sprintf('Current size: %d - Minimal error: %f - discarding: %d\n',
										ncol(TC),min(err),disc))}
					itt<-setdiff(1:ncol(TC),disc);
					TC<-TC[,itt];
				}else{
					dec<-5;
					if(ncol(TC)>k+50){ dec=20; }
					if(ncol(TC)>k+100){ dec=50; }
					if(ncol(TC)>k+200){ dec=100; }
					if(ncol(TC)>k+2000){ dec=1000; }
					sres<-sort(err,index.return=T); serr<-sres$x; iuu<-sres$ix
					if(opt$verbose){ cat(sprintf('Minimal error: %f  Error %d : %f - discarding: %d\n',
										serr[1],dec,serr[dec], iuu[1:dec]))}
					itt<-setdiff(1:ncol(TC),iuu[1:dec]);
					TC<-TC[,itt,drop=F];
				}
			}
			
			# final fit
			T2<-TC; G <- t(T2) %*% T2; A0 <- randsplxmat(ncol(T2),n); W <- t(T2) %*% D;
			# todo
			mqnr<- RHLasso(G,W,A0, myeps);  A2<-mqnr[["A"]]
		}
		
		
	}else{    
		
		if(opt$affine){
			T2<-TC;
			A2 <- affls(T2, D);
			err<-sum(sum(abs(D-T2%*%A2)^2));
			
			if(opt$verbose) cat(sprintf('With %dm topics - err: %f\n',ncol(TC),err))
			
			while(ncol(TC)>r){
				rm(err);
				for(i in 1:ncol(TC)){
					itt<-setdiff(1:ncol(TC),i);
					T2=TC[,itt,drop=F];
					if(rem[i,10]==0){
						if(opt$verbose) cat(sprintf('To do: %d - done: %d\n', ncol(TC),i))
						G <- t(T2) %*% T2; A0 <- randsplxmat(ncol(T2),n); W = t(T2) %*% D;
						mqr<- RQuadSimplex(G,W,A0, myeps); A2<-mqr[["A"]]
						err[i]<-sum(sum(abs(D-T2%*%A2)^2));
					}
				}
				if(ncol(TC)<4*r){
					disc<-which.min(err);
					if(opt$verbose){ cat(sprintf('Current size: %d - Minimal error: %f - discarding: %d\n',
										ncol(TC),min(err),disc))}
					itt<-setdiff(1:ncol(TC),disc);
					TC<-TC[,itt];
				}else{
					dec<-5;
					if(ncol(TC)>k+50){ dec=20; }
					if(ncol(TC)>k+100){ dec=50; }
					if(ncol(TC)>k+200){ dec=100; }
					if(ncol(TC)>k+2000){ dec=1000; }
					sres<-sort(err,index.return=T); serr<-sres$x; iuu<-sres$ix
					if(opt$verbose){ cat(sprintf('Minimal error: %f  Error %d : %f - discarding: %d\n',
										serr[1],dec,serr[dec], iuu[1:dec]))}
					itt<-setdiff(1:ncol(TC),iuu[1:dec]);
					TC<-TC[,itt,drop=F];
				}
			}
			
			# final fit
			T2 <- TC;
			A2 <- affls(T2, D);
			
		}else{    
			
			T2 <- TC;
			A2 <-mldivide(T2,D);
			
			err<-sum(sum(abs(D-T2%*%A2)^2));
			
			if(opt$verbose) cat(sprintf('With %dm topics - err: %f\n',ncol(TC),err))
			
			while(ncol(TC)>r){
				rm(err);
				for(i in 1:ncol(TC)){
					itt<-setdiff(1:ncol(TC),i);
					T2=TC[,itt,drop=F];
					if(rem[i,10]==0){
						if(opt$verbose) cat(sprintf('To do: %d - done: %d\n', ncol(TC),i))
						G <- t(T2) %*% T2; A0 <- randsplxmat(ncol(T2),n); W = t(T2) %*% D;
						mqr<- RQuadSimplex(G,W,A0, myeps); A2<-mqr[["A"]]
						err[i]<-sum(sum(abs(D-T2%*%A2)^2));
					}
				}
				if(ncol(TC)<4*r){
					disc<-which.min(err);
					if(opt$verbose){ cat(sprintf('Current size: %d - Minimal error: %f - discarding: %d\n',
										ncol(TC),min(err),disc))}
					itt<-setdiff(1:ncol(TC),disc);
					TC<-TC[,itt];
				}else{
					dec<-5;
					if(ncol(TC)>k+50){ dec=20; }
					if(ncol(TC)>k+100){ dec=50; }
					if(ncol(TC)>k+200){ dec=100; }
					if(ncol(TC)>k+2000){ dec=1000; }
					sres<-sort(err,index.return=T); serr<-sres$x; iuu<-sres$ix
					if(opt$verbose){ cat(sprintf('Minimal error: %f  Error %d : %f - discarding: %d\n',serr[1],
										dec,serr[dec], iuu[1:dec]))}
					itt<-setdiff(1:ncol(TC),iuu[1:dec]);
					TC<-TC[,itt,drop=F];
				}
			}
			# final fit
			T2 <- TC;
			A2 <- mldivide(T2,D);
		}
		
	}    
	
	err<-sum(sum(abs(D-T2%*%A2)^2));
	T<-T2; A<-A2;
	if(opt$verbose) cat(sprintf('With %d topics - err: %f\n',ncol(TC),err))
	status$nsamples <- nsamples;
	status$err <- err;
	
	return(list(T=T, A=A, status=status))
	
}
#######################################################################################################################
