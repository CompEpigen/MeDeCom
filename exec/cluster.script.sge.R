
if(any(grepl(".*cluster.script.gse.R", commandArgs()))){
	require(MeDeCom)
}

param.file<-commandArgs()[6]
params_list<-readRDS(param.file)

for(params in params_list){
	
	#print(params)
	
	load(file.path(params$DD,"data.set.RData"))
	D<-meth.data

	if(file.exists(file.path(params$DD, "trueA.RData"))){
		load(file.path(params$DD,"trueT.RData"))
		params$trueT<-trueT
		Tstar_present<-TRUE
	}else{
		Tstar_present<-FALSE
	}
	if(file.exists(file.path(params$DD, "trueA.RData"))){
		load(file.path(params$DD, "trueA.RData"))
		params$trueA<-trueA[,params$sample_subset]
		Astar_present<-TRUE
	}else{
		Astar_present<-FALSE
	}
	
	if(file.exists(file.path(params$DD, "start.RData"))){
		load(file.path(params$WD,"start.RData"))
	}
	
	load(file.path(params$DD, sprintf("cg_subset_%d.RDdata",params$cg_subset)))
	params$cg_subset<-cg_subset
	
	if(params$mode %in% c("initial_fine", "cv_fine")){
		if(file.exists(file.path(params$WD, params$lnr_result))){
			load(file.path(params$WD, params$lnr_result))
			old_result<-result
			rm(result)
		}
		
		if(file.exists(file.path(params$WD, params$clnr_result))){
			load(file.path(params$WD, params$clnr_result))
			params$startT<-result$That
			params$startA<-result$Ahat
			
		}
	}
	
	
	if(file.exists(file.path(params$DD, "cv_partitions.RDdata"))){
		load(file.path(params$DD, sprintf("cv_partitions.RDdata")))
	}
	#params$cg_subset<-cg_subsets[[params$cg_subset_id]]
	#params$sample_subset<-sample_subset
	#params$meth_matrix<-D
	
	if(params$mode %in% c("full","initial", "initial_fine")){
		incl_samples<-1:length(params$sample_subset)
	}else{
		fold_subset<-cv.partitions[params$FOLD,]
		incl_samples<-params$sample_subset[-fold_subset]
	}
	
	params$meth_matrix<-D[params$cg_subset,incl_samples,drop=FALSE]
	
	if(params$mode %in% c("full")){
		## average the CV results and use as the init
		cv_files<-params$cv_init_results
		print(cv_files)
		cv_result_list<-vector("list", length(cv_files))
		for(cv_file_idx in 1:length(cv_files)){
			cv.file<-file.path(params$WD, cv_files[[cv_file_idx]])
			if(file.exists(cv.file)){
				load.env<-new.env(parent=emptyenv())
				load(cv.file, envir=load.env)
				cv_result_list[[cv_file_idx]]<-get("result", envir=load.env)
			}
		}
		print(str(cv_result_list))
		cv_result_list<-cv_result_list[!sapply(cv_result_list, is.null)]
		inits<-MeDeCom:::summarizeCVinits(cv_result_list)
		params$startT<-inits$T
		if(!is.null(inits$A)){
			params$startA<-inits$A
		}else{
			params$startA<-MeDeCom:::factorize.regr(params$meth_matrix, params$startT)[["A"]]
		}
	}
	
	
	if(Tstar_present){
		if(!is.null(params$fixed_T_cols)){
			free_cols<-setdiff(1:ncol(trueT_ff), params$fixed_T_cols)
			params$fixedT<-trueT_ff[,params$fixed_T_cols, drop=FALSE]
			#trueT<-trueT_ff[,-fixed_T_cols, drop=FALSE]]
		}else{
			fixedT<-NULL
			free_cols<-1:ncol(trueT_ff)
		}
	}
	
	single_run_params<-intersect(names(as.list(args(MeDeCom:::singleRun))), names(params))
	result<-do.call("singleRun", params[single_run_params], envir=asNamespace("MeDeCom"))
	
	if(params$mode %in% c("full", "initial", "cv") || (result$Fval < old_result$Fval)){
		if(params$mode %in% c("full", "initial", "initial_fine")){
			trueT_prep<-trueA_prep<-NULL
			if(Tstar_present){
				trueT_prep<-trueT[params$cg_subset,free_cols,drop=FALSE]
			}
			if(Astar_present){
				trueA_prep<-trueA_ff[,incl_samples,drop=FALSE]
			}
			perf_result<-MeDeCom:::estimatePerformance(result,
					params$meth_matrix,
					trueT_prep, 
					trueA_prep)
		}else{
			perf_result<-MeDeCom:::estimateFoldError(
					result$That, 
					D[params$cg_subset,params$sample_subset[fold_subset],drop=FALSE],
					params$NFOLDS)			
		}
		for(elt in names(perf_result)){
			result[[elt]]<-perf_result[[elt]]
		}
		
		if(params$mode %in% c("initial_fine", "cv_fine")){
			print("found a better solution")
		}
		save(result, file=file.path(params$WD, params$lnr_result))
	}
}