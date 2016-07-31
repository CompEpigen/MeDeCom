# Routines for selecting cell type-specific quantitative markers
# 
# Author: Pavlo Lutsik
###############################################################################


#
#
filter.by.mean.diff<-function(ref.data, min.diff=0.25, cell.types=colnames(ref.data)){
	
	if(length(cell.types)!=ncol(ref.data)){
		stop("invalid value for cell.types")
	}
	
	meandiff<-sapply(unique(cell.types), function(t){
				
				meandiff<-rowMeans(as.matrix(ref.data[,cell.types==t]))-rowMeans(ref.data[,cell.types!=t])
				meandiff
				
			})
	
	
	candidates<-lapply(unique(cell.types), function(t){
				
				mdiff<-meandiff[,t][abs(meandiff[,t])>min.diff]  
				names(mdiff[order(abs(mdiff), decreasing=T)])
				
			})
	
	names(candidates)<-unique(cell.types)
	
	candidates
		
}

#
#
filter.by.ttest<-function(ref.data, n.cand=10, preselected=NULL, cell.types=colnames(ref.data)){
	
	if(length(cell.types)!=ncol(ref.data)){
		stop("invalid value for cell.types")
	}
	
	
	if(!is.null(preselected)){
		if(!(is.list(preselected) && all(sapply(preselected, is.character)) && length(preselected)==length(unique(cell.types)) )){
			stop("invalid argument for preselected: expected list ")
		}
		
	}	
		
	ttest.stat<-lapply(unique(cell.types), function(t){ 
				apply(ref.data[preselected[[t]],],1,function(r) {
							if(length(r[cell.types==t])>1) {
								t.test(r[cell.types==t], r[cell.types!=t])$statistic 
							}else{ 
								t.test(mu=r[cell.types==t], r[cell.types!=t])$statistic
							}
						})
			})
	
	names(ttest.stat)<-unique(cell.types)
	
	candidates<-sapply(unique(cell.types),function(t){
				
				ord<-order(ttest.stat[[t]], decreasing=T)
				
				names(ttest.stat[[t]])[ord][1:n.cand]
				
			})
	
	colnames(candidates)<-unique(cell.types)
	
	candidates
}

#
#
filter.by.linearity<-function(candidates, data.set){
	
	cell.types<-colnames(candidates)
	
	#minAICcounts<-list()
	
	markers<-lapply(cell.types, function(ct){
				
				cand<-intersect(rownames(data.set),candidates[,ct])
				models<-paste("resp", cand, sep="~")
				models<-c(models, "resp~1")
				
				predictors<-t(data.frame(data.set)[cand,])
				
				design<-list()
				lapply(1:dim(predictors)[2], function(i){
							el(design, where=colnames(predictors)[i])<<-predictors[,i]; return(NULL)
						})
				
				design$resp<-t(data.set)
				logliks<-sapply(models, function(f) apply(residuals(lm(f, design)),2,logLik2))
				
				model.size<-sapply(strsplit(models,"+", fixed=T),  length)
				model.size[grep("~1", models)]<-0
				aics<-sapply(1:length(models), function(i) -2*logliks[,i]+2*model.size[i])
				
				bestModel<-models[apply(aics,1,which.min)]
				names(bestModel)<-rownames(aics)
				
				counts<-table(bestModel)
				#el(minAICcounts, where=ct)<<-counts
			
				counts["resp~1"]<-0
				counts<-counts[order(-counts)]
				print("Best models for:")
				print(ct)
				print(counts)
					
				return(strsplit(names(counts)[1], "~")[[1]][2])
				
			})
			
	markers
	
}


