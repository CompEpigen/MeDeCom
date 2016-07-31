#######################################################################################################################
#
# Routines for simulation of the DNA methylation profiles
#
# @author Pavlo Lutsik
#######################################################################################################################

#
# simulate.source.profile
#
# Generate the source DNA methylation profile
#
simulate.source.profile<-function(m, vals=c(0,1), probs=c(0.5,0.5)){

	sample(vals, m, replace=TRUE)
	
}
#
# simulate.ct.profile
#
# Generate the characteristic DNA methylation profiles of cell types
#
simulate.ct.profiles<-function(source, similarities){
	 
	CT.AVG<-sapply(similarities,function(sim){
				
				idx<-sample.int(length(source),floor((1-sim)*length(source)))
				CT.M<-source
				CT.M[idx]<-sample(c(0,1), floor((1-sim)*length(source)), replace=TRUE)
				CT.M
				
			})
	
	CT.AVG
  	
}

# simulate.ind.profiles
#
# Generate the characteristic DNA methylation profiles of cell types
#

simulate.ind.profiles<-function(
		ctm, 
		n,
		method="binom",
		vars=NULL, 
		chi.sq.df=0.11,
		success.prob=NULL,
		alpha=0.1,
		beta=10,
		vals=c(0,1)){

#	if(is.null(CHI.SQ.DF)){
#		CHI.SQ.DF=(20/n)
#	}
	
	
	if(method=="normal"){
		
		if(is.null(vars)){
			sds<-sqrt(rchisq(nrow(ctm), df=chi.sq.df)/n)
		}
		
		profiles<-lapply(1:ncol(ctm), function(iid){
					
					cell.sample<-sapply(1:m, function(cgi){
								
								cg.vector<-rnorm(n, mean=ctm[cgi,iid], sd=sds[cgi]*scale)
								
							})
					cell.sample<-MeDeCom:::projectV(cell.sample, vals)
					t(cell.sample)
					
				})
	}else if(method=="binom"){
		
		#success.prob<-scale*rchisq(nrow(ctm), df=chi.sq.df)/n
		
		#success.prob[success.prob>0.5]<-0.5
		
		profiles<-lapply(1:ncol(ctm), function(iid){
					
				if(is.null(success.prob))
					success.prob<-rbeta(nrow(ctm), shape1=alpha, shape2=beta)	
					
				ind.prof<-repmat(ctm[,iid,drop=F],m=n,n=1)
					
				change.v<-t(sapply(success.prob, function(sp) rbinom(n,1,prob=sp)))	
				
				ind.prof[change.v==1]<-sign(1-ind.prof[change.v==1])
								
				ind.prof
				})
		
	}
	
	
	profiles<-lapply(1:n, function(i) sapply(1:ncol(ctm), function(j) profiles[[j]][,i] ))
	
	
	profiles
	
}

# prepare.average.profile
#
# Prepare individual profiles from existing data 
#

prepare.average.profile<-function(pure.ct, 
		samp,
		o=0.45, 
		mean=F){
	
	n<-nrow(pure.ct)
	
	sds.samp<-apply(pure.ct[samp,,drop=F],1,sd)
	
	proj.profile.samp<-pure.ct[samp,,drop=F]
	if(mean){
		proj.profile.samp<-matrix(rowMeans(proj.profile.samp))
	}
	proj.profile.samp[proj.profile.samp<o]<-0
	proj.profile.samp[proj.profile.samp>1-o]<-1
	proj.profile.samp[proj.profile.samp>=o & proj.profile.samp<=1-o]<-0.5
	
	
	
	return(list(proj.profile.samp, sds.samp))
	
}

# simulate.populations
#
# Generate the cell populations
#
simulate.populations<-function(ind.profiles,
		sds.samp=NULL,
		vals=c(0,1),
		method="binom",
		success.prob=NULL,
		m=nrow(ind.profiles),
		nc=1000,
		alpha=0.1,
		beta=10){
		
		if(method=="normal")
		{
			profiles<-lapply(1:ncol(ind.profiles), function(iid){
				
				cell.sample<-sapply(1:m, function(cgi){
							
							cg.vector<-rnorm(nc, mean=ind.profiles[cgi,iid], sd=sds.samp[cgi]*scale)
							
						})
				cell.sample<-MeDeCom:::projectV(cell.sample, vals)
				t(cell.sample)
				
			})
	
		}else if(method=="binom"){
			
			#success.prob<-scale*rchisq(nrow(ind.profiles), df=chi.sq.df)/length(ind.profiles)
			
			if(is.null(success.prob))
				success.prob<-rbeta(nrow(ind.profiles), shape1=alpha, shape2=beta)
			
			#success.prob[success.prob>0.5]<-0.5
			
			profiles<-lapply(1:ncol(ind.profiles), function(iid){
						
						pop.prof<-repmat(ind.profiles[,iid,drop=F],m=nc,n=1)
						
						change.v<-t(sapply(success.prob, function(sp) rbinom(nc,1,prob=sp)))	
						
						pop.prof[change.v==1]<-sign(1-pop.prof[change.v==1])
						
						pop.prof
						
					})
			
		}
		profiles
}


# mix.populations
#
# Mix the cell populations according in proportions,
# given by the mixing matrix
#
mix.populations<-function(
		populations
		,mixing.matrix
		,cell.subset.size=NULL
		,noize=0){
	
	if(!is.list(populations))
		stop("Invalid value for populations")
	
	
	if(ncol(mixing.matrix)!=length(populations))
		stop("The second dimension of the mixing matrix has to the number of cell types")
		
	nt<-length(populations[[1]])
	
	results<-vector("list", ncol(mixing.matrix))
	
	if(length(populations[[1]])==nrow(mixing.matrix)){
		
		for(j in 1:ncol(mixing.matrix)){
			nc<-ncol(populations[[j]][[1]])
			
			if(is.null(cell.subset.size))
				css<-nc
			else
				css<-cell.subset.size
			
			subs.sizes<-sapply(mixing.matrix[,j], function(fr) floor(css*fr))
			cell.subsets<-lapply(1:nt,function(x) sample.int(ncol(populations[[j]][[nt]]), subs.sizes[x])) 
			
			
			mix<-list()
			for(i in 1:nt){
				mix[[i]]<-populations[[j]][[i]][,cell.subsets[[i]]]				
			}
			
			results[[j]]<-do.call("cbind", mix)
			rm(mix)
		}
			
		
	}else{
			warning("Not implemented yet")
			results<-NULL
		
	}
	
	return(results)
		
}

# introduce.imprinting
#
# Simulate genomic imprinting
#

introduce.imprinting<-function(populations,
		fraction=0.02,
		ixx=NULL){
	
	ixx<-sample.int(nrow(populations[[1]]),floor(fraction*nrow(populations[[1]])))
	

	for(i in 1:length(populations))
	{
		populations[[i]][ixx,]<-0.5
	}
	
	populations
	
}

# introduce.asm
#
# Simulate allele-specific methylation
#
introduce.asm<-function(populations,
		fraction=0.1,
		sites=NULL){
	
	
	for(i in 1:length(populations))
	{
		ixx<-sample.int(nrow(populations[[i]]),floor(fraction*nrow(populations[[i]])))
		snp.vals<-sample(c(0,0.5,1), length(ixx), replace=TRUE)
	
		populations[[i]][ixx,]<-snp.vals
		
	}

	populations
}

# introduce.effects
#
# Simulate true biological effects
#
introduce.effects<-function(populations,
		types,
		sites,
		signs,
		means, 
		sd.mean.frac=0.1){
	
	
	for(i in 1:length(populations)){
	
		for(j in types){
		
		effects<-rnorm(length(sites[[j]]), mean=means[j], sd=sd.mean.frac*means[j])		
		if(is.null(signs[[j]]))
			signs[[j]]<-sign(rnorm(length(sites[[j]])))
		
		new.values<-sapply(1:length(sites[[j]]), function(si){
					
				if(signs[[j]][si]==1){
					
					new.value<-sapply(populations[[i]][[j]][sites[[j]][si],], function(value){
						if(value %in% c(0,0.5)){
							if(runif(1)<effects[si]/(1-value)){
							    return(1)	
							}else{
								return(value)
							}
						}else{
							return(value)
						}
					})
										
				}else{
					new.value<-sapply(populations[[i]][[j]][sites[[j]][si],], function(value){
								if(value %in% c(0.5,1)){
									if(runif(1)<effects[si]/(value)){
										return(0)	
									}else{
										return(value)
									}
								}else{
									return(value)
								}
							})
					
				}		
					
				})
		
		populations[[i]][[j]][sites[[j]],]<-t(new.values)
		
		}
	}
	
	return(populations)	
}

summarize.patterns<-function(cell.pops){
	
	up_pops<-lapply(cell.pops, unique, MARGIN=2)
	
	ups<-do.call("cbind", up_pops)
	
	up<-unique(ups, MARGIN=2)
	
	#large.matrix<-do.call("cbind", cell.pops)
	
	fmat<-matrix(0, nrow=ncol(up), ncol=length(cell.pops))
	
	for(upi in 1:ncol(up)){
		
		for(popi in 1:length(cell.pops)){
			fmat[upi,popi]<-sum(colSums(abs(cell.pops[[popi]]-up[,upi]))==0)/ncol(cell.pops[[popi]])
		}
	}
	
	return(list(C=up, F=fmat))
}

# simulate.450k
#
# Simulate 450k microarray measurement 
#

simulate.450k<-function(
		cell.pops
		,ncs=ncol(cell.pops[[1]])
		,technical.effects=TRUE
		,totalInt=20000
		,totalInt.sd.stable=2000
		,totalInt.sd.ind=0
		,meth.channel.sd=300
		,umeth.channel.sd=300
		,meth.channel.bg.mean=1000
		,umeth.channel.bg.mean=1000
		,meth.channel.bg.sd=300
		,umeth.channel.bg.sd=300
){
	
	ti<-rnorm(nrow(cell.pops[[1]]), mean=totalInt, sd=totalInt.sd.stable)
	
	profiles<-sapply(cell.pops, function(ip){
				
				nc<-ncol(ip)
				
				cell.sample<-sample.int(nc, ncs)
				
				cell.sample.profile<-apply(ip,1,"[", cell.sample)
				
				profile<-rowMeans(t(cell.sample.profile))
				
				
				#			profile.sds<-beta.sd.model(profile)
				#			
				#			profile<-sapply(1:length(profile),function(j){
				#						rnorm(1,mean=profile[j],sd=beta.sd.model(profile[j]))
				#					})
				
				
				if(technical.effects){
					
					tii<-sapply(ti, rnorm, n=1, sd=totalInt.sd.ind)
					
					meth<-tii*profile
					umeth<-tii*(1-profile)
					
					meth<-sapply(meth, function(mv) rnorm(1,mv,meth.channel.sd)+rnorm(1,meth.channel.bg.mean,meth.channel.bg.sd))
					umeth<-sapply(umeth, function(uv) rnorm(1,uv,umeth.channel.sd)+rnorm(1,umeth.channel.bg.mean,umeth.channel.bg.sd))
					
					meth[meth<0]<-1
					umeth[umeth<0]<-1
					
					measured.profile<-meth/(meth+umeth)
					
					return(measured.profile)	
				}
				
				return(profile)
				
			})
	
	profiles
	
}

#
# get.confusion.table
#
# Get the confusion table from CpG site classification results 
#
get.confusion.table<-function(pvals, true.sites, sl=0.05, adjust=TRUE, adjust.method="hochberg"){

	if(adjust){
		pvals<-p.adjust(pvals, adjust.method)
	}
	
	## confusion table
	predict<-pvals<sl
	true<-rep(FALSE, length(pvals))
	true[true.sites]<-TRUE
	
	conf.table<-table(Outcome = predict, Truth = true)[2:1,2:1]
	
	sens<-conf.table[1,1]/(conf.table[1,1]+conf.table[2,1])
	spec<-conf.table[2,2]/(conf.table[1,2]+conf.table[2,2])
	
	return(list(conf.table=conf.table, sensitivity=sens, specificity=spec))
	
}