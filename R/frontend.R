#######################################################################################################################
# 
#  frontend.R
#
#  Frontend functionality for DNA methylation deconvolution 
# 
#######################################################################################################################
#'
#' runMeDecom
#'
#' Perform a MeDeCom experiment
#'
#' @param data	 	DNA methylation dataset as a \code{numeric} matrix (methylation sites vs samples) or an ojbect of class \code{RnBeadSet}
#' @param Ks		values of parameter \deqn{k} to be tested, vector of type \code{integer}
#' @param lambdas	values of parameter \deqn{\lambda} to be tested, vector of type \code{numeric}
#' @param opt.method optimization method used. Currently supported values are \code{"MeDeCom.quadPen"} and \code{"MeDeCom.cppTAfact"}.
#' @param cg_subsets a \code{list} of \code{integer} vectors specifying row indices to include into the analysis
#' @param sample_subset samples to include into the analysis
#' @param startT    a \code{list} of length equal to \code{length(Ks)} or a \code{matrix} with \code{max(Ks)} columns
#' @param startA    a \code{list} of length equal to \code{length(Ks)} or a \code{matrix} with \code{max(Ks)} rows
#' @param trueT		a numeric matrix with as many rows as there are methylation sites in \code{data}
#' @param trueA		a numeric matrix with as many columns as there are methylation sites in \code{data}
#' @param fixed_T_cols  columsn of T which are known (to be implemented)
#' @param NINIT		number of random initializations
#' @param ITERMAX   maximal number of iterations of the alternating optimization scheme
#' @param NFOLDS    number of cross-validation folds
#' @param N_COMP_LAMBDA   the number of solutions to compare in the "smoothing" step
#' @param NCORES    number of cores to be used in the parallelized steps (at best a divisor of NINIT)
#' @param analysis.name a deliberate name of the analysis as a \code{character} singleton
#' @param use.ff    use \code{ff} package functionality for memory optimization
#' @param cluster.settings  a list with parameters for an HPC cluster
#' @param temp.dir  a temporary directory for the cluster-based analysis available on all nodes
#' @param cleanup	if \code{TRUE} the temporary directory will be removed on completion
#' @param verbosity verbosity level, \code{0} for the quiet execution 
#' @param time.stamps add timestamps to the diagnostic output
#' 
#' @details 
#' 
#' @return MeDeComSet object
#'
#' @author Pavlo Lutsik
#' @export
#' 
runMeDeCom<-function(
		data, 
		Ks, 
		lambdas,
		opt.method="MeDeCom.quadpen",
		cg_subsets=NULL,
		sample_subset=NULL,
		startT=NULL,
		startA=NULL,
		trueT=NULL,
		trueA=NULL,
		fixed_T_cols=NULL,
		NINIT=100, 
		ITERMAX=1000, 
		NFOLDS=10,
		N_COMP_LAMBDA=4,
		NCORES=1,
		analysis.name=NULL,
		use.ff=FALSE,
		cluster.settings=NULL,
		temp.dir=NULL,
		cleanup=TRUE,
		verbosity=1L,
		time.stamps=FALSE
){
	ts<-function(){
		if(time.stamps){
			paste0(format(Sys.time(),"%Y-%m-%d %H:%M:%S"), ", ")
		}else{
			""
		}
	}
	if(verbosity>0){
		cat(sprintf("[%sMain:] checking inputs\n",ts()))
	}
	if(inherits(data, "RnBSet")){
		D<-meth(data)
		pheno<-pheno(data)
	}else if(inherits(data, "list")){
		D<-data[[1]]
		pheno<-data[[2]]
	}else if(inherits(data, "matrix")){
		D<-data
	}else{
		stop("unsupported object supplied for data")
	}
	cluster_run<-!is.null(cluster.settings)
	
	#### data preparation
	if(verbosity>0){
		cat(sprintf("[%sMain:] preparing data\n",ts()))
	}
	LAMBDA_GRID<-lambdas
	if(length(LAMBDA_GRID)<4L){
		N_COMP_LAMBDA=0
	}
	if(is.null(cg_subsets)){
		cg_subsets=list("All CpGs"=1:nrow(D))
	}
	if(is.null(sample_subset)){
		sample_subset=1:ncol(D)
	}
	if(is.null(analysis.name)){
		analysis.name<-sprintf("MeDeCom_run_%s", format(Sys.time(),"%Y%m%d_%H_%M_%S"))
	}
	if(is.null(temp.dir)){
		temp.dir<-file.path(tempdir(), "MeDeCom")
		if(!file.exists(temp.dir)){
			dir.create(temp.dir)
		}
	}
	
	#### marker selection
#	if(!is.null(top.var)){
#		if(!(is.integer(top.var) || length(top.var)!=1 || top.var<2 || top.var>nrow(D))){
#			stop("invalid value for top.var supplied")
#		}
#		vars<-apply(D, var)
#	}
	
	#### final data matrix
	if(use.ff){
		D_ff<-as.ff(D)
	}else{
		D_ff<-D
	}
	
	Tstar_present<-!is.null(trueT)
	if(Tstar_present){
		if(use.ff){
			trueT_ff<-as.ff(trueT)
		}else{
			trueT_ff<-trueT
		}
	}
	
	Astar_present<-!is.null(trueA)
	if(Astar_present){
		trueA_ff<-trueA[,sample_subset]
	}else if(Tstar_present){
		regr.est<-MeDeCom:::factorize.regr(D[,sample_subset], trueT)
		trueA_ff<-regr.est$A
		rm(regr.est)
	}
	
	### cross-validation preparation
	cv.partitions<-do.call("rbind", lapply(0:(NFOLDS-1), function(fld) (1:(length(sample_subset)%/%NFOLDS))+fld*(length(sample_subset)%/%NFOLDS)))
	
	#out <- cv.partitions[FOLD,];
	#inn <- setdiff(1:length(sampl_subset), out);
	#D_in  <- D[,inn,drop=FALSE]
	#D_out <- D[,out,drop=FALSE]
	
	if(cluster_run){
	#### cluster preparation
		
		wd<-file.path(temp.dir, analysis.name)
		if(!file.exists(wd)){
			dir.create(wd)
		}
		WD<-wd
		DD<-file.path(WD,"data")
		if(!file.exists(DD)){
			dir.create(DD)
		}
	}else{
		WD<-DD<-temp.dir
	}
	#### analysis run
	if(verbosity>0){
		cat(sprintf("[%sMain:] preparing jobs\n",ts()))
	}
#	run_param_list<-vector("list",
#					length(cg_subsets)*
#					length(Ks)*
#					(NFOLDS+1)*
#					(
#						length(LAMBDA_GRID)+
#						2*(
#							(length(LAMBDA_GRID)-N_COMP_LAMBDA)*N_COMP_LAMBDA + 
#							(N_COMP_LAMBDA)*(N_COMP_LAMBDA-1)/2
#						)
#					)
#				)
	lambda_smoothing<-N_COMP_LAMBDA>0
	if(lambda_smoothing){
		run_count<-length(cg_subsets)*
				length(Ks)*
				(
					(NFOLDS+1)*length(LAMBDA_GRID)+
					(NFOLDS)*
					(
						2*(
							(length(LAMBDA_GRID)-N_COMP_LAMBDA)*N_COMP_LAMBDA + 
							(N_COMP_LAMBDA)*(N_COMP_LAMBDA-1)/2
							)
						)
					)
	}else{
		run_count<-length(cg_subsets)*
				length(Ks)*
				(
					(NFOLDS+1)*length(LAMBDA_GRID))
	}
	run_param_list<-vector("list", run_count)
				
	result_index<-array(dim=c(
					length(cg_subsets),
					(NFOLDS+1),
					length(Ks),
					length(LAMBDA_GRID)))
	result_list<-vector("list", 
					length(cg_subsets)*
					length(Ks)*
					(NFOLDS+1)*
					length(LAMBDA_GRID))
	param_ctr<-0
	res_ctr<-0
	#run_param_list<-list()
	cg_subset_ids<-1:length(cg_subsets)
	#for(run in c("initial", "cv"))
	for(run in c("cv", "full"))
	{
		#for(K in Ks){
		for(gr in cg_subset_ids)
		{
			if(run=="cv"){
				folds=1:NFOLDS
			}else{
				folds=0
			}
			for(fold in folds)
			{	
				for(K in Ks)
				{
					for(lnr in 1:length(LAMBDA_GRID))
					{
						param_ctr<-param_ctr+1
						res_ctr<-res_ctr+1
						param_list<-list(
								cg_subset_id=gr,
								sample_subset=sample_subset,
								K=K, 
								lambda=LAMBDA_GRID[lnr],
								lambda_nr=lnr, 
								lambda_cnr=lnr,
								NINIT=NINIT, 
								NFOLDS=NFOLDS, 
								FOLD=fold, 
								NFOLDS=NFOLDS,
								mode=run, 
								direction="topbottom", 
								ITERMAX=ITERMAX,
								NCORES=NCORES,
								WD=WD, 
								DD=DD,
								METHOD=opt.method)
						if(run=="full" && cluster_run){
							cv_results_idx<-result_index[
									match(param_list$cg_subset_id, cg_subset_ids),
									match(1:NFOLDS, c(0,1:NFOLDS)),
									match(param_list$K, Ks),
									param_list$lambda_nr]
							param_list$cv_init_results<-result_list[cv_results_idx]
						}
						run_param_list[[param_ctr]]<-param_list
						if(run=="full"){
							jname<-sprintf("init_%s_%d_%d_%d", analysis.name, gr, K, lnr)
						}else{
							jname<-sprintf("cv_%s_%d_%d_%d_%d", analysis.name, gr, fold, K, lnr)
						}

						attr(run_param_list[[param_ctr]], "jname")<-jname
						if(run=="full"){
							if(lambda_smoothing){
								deps<-sprintf("cvfine_%s_%d_*_%d_bt_%d*", analysis.name, gr, K, lnr)
							}else{
								deps<-sprintf("cv_%s_%d_*_%d_%d*", analysis.name, gr, K, lnr)
							}
							attr(run_param_list[[param_ctr]], "depends_on")<-deps
						}
						if(cluster_run){
							result_file<-sprintf("result_%s_GROUP_%d_K_%d_FOLD_%d_LAMBDA_%g_NINIT_%d.RData", 
									run,
									gr,
									K,
									fold,
									LAMBDA_GRID[lnr],
									NINIT)
							result_list[[res_ctr]]<-result_file
							run_param_list[[param_ctr]]$lnr_result<-result_file
						}else{
							result_list[[res_ctr]]<-list()
						}
						result_index[
								#match(run, c("initial", "cv")),
								match(gr, cg_subset_ids),
								match(fold, c(0, 1:NFOLDS)),
								match(K, Ks),
								lnr]<-res_ctr
					}
					if(run=="cv" && lambda_smoothing){
						for (direction in c("topbottom", "bottomtop")){
						if(direction=="topbottom"){
							lindexes<-(length(LAMBDA_GRID)-1):1
						}else{
							lindexes<-(2:length(LAMBDA_GRID))
						}
						for(lnr in lindexes){
								if(direction=="topbottom"){
									comparison_lambdas<-(lnr+1):min(length(LAMBDA_GRID),(lnr+N_COMP_LAMBDA))
									token="tb"
								}else{
									comparison_lambdas<-max(1L,(lnr-N_COMP_LAMBDA)):(lnr-1)	
									token="bt"
								}
								for(comp_lambda in comparison_lambdas){
									param_ctr<-param_ctr+1
									res_idx<-result_index[
											#match(run, c("initial", "cv")),
											match(gr, cg_subset_ids),
											match(fold, c(0, 1:NFOLDS)),
											match(K, Ks),
											lnr]
									comp_res_idx<-result_index[
											#match(run, c("initial", "cv")),
											match(gr, cg_subset_ids),
											match(fold, c(0, 1:NFOLDS)),
											match(K, Ks),
											comp_lambda]
									run_param_list[[param_ctr]]<-list(
											cg_subset_id=gr,
											sample_subset=sample_subset,
											K=K,
											lambda=LAMBDA_GRID[lnr],
											lambda_nr=lnr, 
											NINIT=NINIT, 
											NFOLDS=NFOLDS,
											FOLD=fold,
											mode=sprintf("%s_fine",	run),
											direction=direction,
											lambda_cnr=comp_lambda,
											ITERMAX=ITERMAX,
											NCORES=NCORES,
											WD=WD,
											DD=DD,
											METHOD=opt.method,
											lnr_result=result_list[[res_idx]],
											clnr_result=result_list[[comp_res_idx]]
											)
									
									#if(run=="initial"){
									if(run=="full"){
										jname<-sprintf("initfine_%s_%d_%d_%s_%d_%d", analysis.name, gr, K, token, lnr, comp_lambda)
										deps<-sprintf("init_%s_%d_%d*", analysis.name, gr, K)
									}else{
										jname<-sprintf("cvfine_%s_%d_%d_%d_%s_%d_%d", analysis.name, gr, fold, K, token, lnr, comp_lambda)
										deps<-sprintf("cv_%s_%d_%d_%d*", analysis.name, gr, fold, K)
									}
									if(lnr>1){
										deps<-c(deps, attr(run_param_list[[param_ctr-1]], "jname"))
									}
									attr(run_param_list[[param_ctr]], "jname")<-jname
									attr(run_param_list[[param_ctr]], "depends_on")<-deps
								}
							}
						}
					}
				}
			}
		}
	}
	if(verbosity>0){
		cat(sprintf("[%sMain:] %d factorization runs in total\n", ts(), length(run_param_list)))
		if(cluster_run){
			cat(sprintf("[%sMain:] starting the jobs\n",ts()))
		}
	}
	if(cluster_run){
		
		#submitAllClusterJobs(param_list)
		meth.data<-D
		save(meth.data, file=file.path(DD, "meth.data.RData"))
		
		save(meth.data, file=file.path(DD,"data.set.RData"))
		#params$meth_matrix<-meth.data
		
		save(trueT, file=file.path(DD,"trueT.RData"))
		if(!is.null(trueA)){
			save(trueA, file.path(DD, "trueA.RData"))
		}
		
		if(!is.null(startT) && !is.null(startA)){
			start<-list(startT, startA)
			save(start, file=file.path(DD,"start.RData"))
		}
		
		for(cgsi in cg_subset_ids){
			cg_subset<-cg_subsets[[cgsi]]
			save(cg_subset, file=file.path(DD, sprintf("cg_subset_%d.RDdata",cgsi)))
		}
		
		save(sample_subset, file=file.path(DD, sprintf("sample_subset.RDdata")))
		save(cv.partitions, file=file.path(DD, sprintf("cv_partitions.RDdata")))
		
		for(idx in 1:length(run_param_list)){
			id<-attr(run_param_list[[idx]], "jname")
			deps<-attr(run_param_list[[idx]], "depends_on")
			submitClusterJob(id, deps, run_param_list[idx], WD, 
					RDIR=cluster.settings$R_bin_dir, hosts=cluster.settings$host_pattern, ram_limit=cluster.settings$mem_limit)
		}
		
		waitForClusterJobs(analysis.name, verbose=verbosity>0L)
		
		for(idx in 1:length(result_list)){
			if(!is.null(result_list[[idx]]) && file.exists(file.path(WD, result_list[[idx]]))){
				load_env<-new.env(parent=emptyenv())
				load(file.path(WD, result_list[[idx]]), envir=load_env)
				result_list[[idx]]<-get("result", envir=load_env)
			}
		}
		
		if(cleanup){
			unlink(WD, recursive=TRUE)
		}
	}else{
		final_results<-vector("list", length(run_param_list))
		
		for(idx in 1:length(run_param_list)){
			
			if(verbosity>0L && idx %% 100==0){
				cat(sprintf("[%sMain:] %d runs complete\n",ts(), idx))
			}
			
			params<-run_param_list[[idx]]
			params$cg_subset<-cg_subsets[[params$cg_subset_id]]
			params$sample_subset<-sample_subset
			#params$meth_matrix<-D
		
			if(params$mode %in% c("full", "initial", "initial_fine")){
				incl_samples<-1:length(sample_subset)
			}else{
				fold_subset<-cv.partitions[params$FOLD,]
				incl_samples<-params$sample_subset[-fold_subset]
			}
			
			params$meth_matrix<-D_ff[params$cg_subset,incl_samples,drop=FALSE]
			#params$trueT<-trueT
			if(Tstar_present){
				params$trueT<-trueT_ff[params$cg_subset,]
				colnames(params$trueT)<-rownames(params$trueT)<-NULL
				params$trueA<-trueA_ff[,params$sample_subset]
				colnames(params$trueA)<-rownames(params$trueA)<-NULL
			}
			
			res_idx<-result_index[
					#match(run, c("initial", "cv")),
					match(params$cg_subset_id, cg_subset_ids),
					match(params$FOLD, c(0, 1:NFOLDS)),
					match(params$K, Ks),
					params$lambda_nr]
			comp_res_idx<-result_index[
					#match(run, c("initial", "cv")),
					match(params$cg_subset_id, cg_subset_ids),
					match(params$FOLD, c(0, 1:NFOLDS)),
					match(params$K, Ks),
					params$lambda_cnr]
			params$lnr_result<-result_list[[res_idx]]
					
			if(params$mode %in% c("initial_fine", "cv_fine")){
				params$startT<-result_list[[comp_res_idx]]$T
				params$startA<-result_list[[comp_res_idx]]$A
			}
			
			if(params$mode %in% c("full")){
				## average the CV results and use as the init
				cv_results_idx<-result_index[
						match(params$cg_subset_id, cg_subset_ids),
						match(1:NFOLDS, c(0,1:NFOLDS)),
						match(params$K, Ks),
						params$lambda_cnr
						]
				inits<-summarizeCVinits(result_list[cv_results_idx])
				params$startT<-inits$T
				if(!is.null(inits$A)){
					params$startA<-inits$A
				}else{
					params$startA<-factorize.regr(params$meth_matrix, params$startT)[["A"]]
				}
			}
			
			if(Tstar_present && !is.null(fixed_T_cols)){
				free_cols<-setdiff(1:ncol(trueT_ff), fixed_T_cols)
				params$fixedT<-trueT_ff[,fixed_T_cols, drop=FALSE]
				#trueT<-trueT_ff[,-fixed_T_cols, drop=FALSE]]
				fixed_T_cols<-integer()
			}else{
				fixedT<-NULL
				free_cols<-1:K
			}
			
			####################### START FACTORIZATION RUN
			single_run_params<-intersect(names(as.list(args(singleRun))), names(params))
			result<-do.call("singleRun", params[single_run_params])
			####################### END FACTORIZATION RUN
			
			if(params$mode %in% c("full", "cv") || (result$Fval < result_list[[comp_res_idx]]$Fval)){
				if(params$mode %in% c("full", "initial_fine")){
						trueT_prep<-trueA_prep<-NULL
						if(Tstar_present){
							trueT_prep<-trueT_ff[params$cg_subset,-fixed_T_cols,drop=FALSE]
						}
						if(Astar_present){
							trueA_prep<-trueA_ff[,incl_samples,drop=FALSE]
						}
						perf_result<-estimatePerformance(result, 
								trueT_prep, 
								trueA_prep)
					}else{
						perf_result<-estimateFoldError(
								result$That, 
								D_ff[params$cg_subset,params$sample_subset[fold_subset],drop=FALSE],
								NFOLDS)			
					}
				for(elt in names(perf_result)){
					result[[elt]]<-perf_result[[elt]]
				}
				
				result_list[[res_idx]]<-result
				if(params$mode %in% c("initial_fine", "cv_fine") && verbosity>1L){
					cat("found a better solution\n")
				}
			}
		}
	}
	
	if(verbosity>0L){
		cat(sprintf("[%sMain:] finished all jobs. Creating the object\n",ts()))
	}
	
	#### object creation
	
	#system(sprintf("%s/Rscript %s/collect.results.R %s %d %d %d", RDIR, SRCDIR, WD, K, NINIT, NFOLDS))
	result_object<-collectResults(result_list, cg_subset_ids, Ks, lambdas, NFOLDS, result_index)
	
	dataset_info<-list(m=nrow(D), n=ncol(D))
	
	params_to_save<-c("NFOLDS","N_COMP_LAMBDA","NINIT","ITERMAX")
	for(param in params_to_save){
		result_object$parameters[[param]]<-get(param)
	}
	
	MeDeComSet(result_object$parameters, result_object$outputs, dataset_info)
}
#######################################################################################################################
submitClusterJob<-function(job_name, dependencies, params, WD, RDIR="/usr/bin", hosts="*", ram_limit="5G"){
	
	src_file<-system.file("exec/cluster.script.sge.R", package="MeDeCom")
	param_file<-file.path(WD, sprintf("%s_params.RDS", job_name))
	saveRDS(params, param_file)
	
	if(!is.null(dependencies)){
		deps<-paste(dependencies, collapse=",")
		deps_string<-sprintf("-hold_jid '%s'", deps)
	}else{
		deps_string<-NULL
	}
	qsub_string<-sprintf("qsub -cwd -j y -o %s.log -b y -V -N %s -l h='%s' -l mem_free=%s", 
										file.path(WD,job_name),	job_name, hosts,  ram_limit)
	if(!is.null(deps_string)){
		qsub_string<-paste(qsub_string, deps_string)
	}
	script_string<-sprintf("%s/Rscript %s %s", RDIR, src_file, param_file)
	
	job_cmd<-paste(qsub_string, script_string)
	res<-system(job_cmd, intern=TRUE)
}
#######################################################################################################################
waitForClusterJobs<-function(analysis_id, lookup_int=10, verbose=TRUE){
	
	repeat{
		lookup_cmd<-sprintf("qstat -r | grep \"Full jobname\" | grep -e %s", analysis_id)
		suppressWarnings({
			running_jobs<-system(lookup_cmd, intern=TRUE)
		})
		if(length(running_jobs)>0){
			cat(sprintf("[SGE jobs:] %d remaining\n", length(running_jobs)))
			Sys.sleep(lookup_int)
		}else{
			break;
		}
	}
}
#######################################################################################################################
#
# A single factorization run with fixed parameters
#
#######################################################################################################################
singleRun<-function(
		meth_matrix,
		K,
		lambda,
		startA=NULL,
		startT=NULL,
		fixedT=NULL,
		Alower=NULL,
		Aupper=NULL,
		fixed_T_cols=integer(),
		NINIT,
		ITERMAX,
		NCORES=1L,
		METHOD="MeDeCom.quadPen",
		verbosity=1L
){
	D<-meth_matrix

		ALG_LAMBDA=lambda#LAMBDA_GRID[lambda_nr];

		if(!is.null(startT) && !is.null(startA)){
			ALG_INT="fixed"
			ALG_OPT=list(T=startT, A=startA)
		}else{
			ALG_INT="random"
			ALG_OPT=NINIT;
		}
		
		fr<-factorize.alternate(	
				D,
				k=K,
				method=METHOD,
				t.method=MeDeCom:::T_METHODS[METHOD],
				lambda=ALG_LAMBDA,
				init=ALG_INT,
				opt=ALG_OPT,
				itermax=ITERMAX,
				qp.Alower=Alower,
				qp.Aupper=Aupper,
				Tfix=fixedT,
				ncores=NCORES
				);
		
	fr$cve<-NA
	names(fr)<-c("That", "Ahat", "Fval", "Conv", "RMSE", "cve")
	return(fr[c(1:2,6,3,5)])
}
#######################################################################################################################
#
# Summarize the results of the CV runs for initialization
#
#######################################################################################################################
summarizeCVinits<-function(result_list){
	nfolds<-length(result_list)
	#Ts<-vector("list", length(result_list))
	#As<-vector("list", length(result_list))
	fold_errs<-vector("numeric", nfolds)
	for(idx in 1:nfolds){
		fold_errs[idx]<-result_list[[idx]]$cve
	}
	best_result<-which.min(fold_errs)
	refT<-result_list[[best_result]]$That
	AvgT<-refT
	incl_folds<-1
	for(idx in (1:nfolds)[-best_result]){
		perm<-matchLMCs(result_list[[idx]]$That, refT)
		## add result to the average
		if(!is.null(perm)){
			incl_folds<-incl_folds+1
			AvgT<-AvgT+result_list[[idx]]$That[,perm,drop=FALSE]
		}
	}
	AvgT<-AvgT/incl_folds
	### TODO: also summarize As
	AvgA<-NULL
	return(list(T=AvgT, A=AvgA))
}
#######################################################################################################################
#
# Estimate performance of a factorization run
#
#######################################################################################################################
estimatePerformance<-function(fr, trueT=NULL, trueA=NULL){
	
				#if(mode %in% c("initial", "initial_fine") ){
				#### Performance estimation part
				if(!is.null(trueT)){
					perm<-MeDeCom:::matchLMCs(fr$T, trueT, check=FALSE)
					print(perm)
					rmseT<-MeDeCom:::RMSE_T(fr$T, trueT, perm)
					
					if(length(unique(perm))>0){
						if(length(unique(perm))<ncol(trueT)){
							
							trueT_adapted<-trueT[,unique(perm),drop=FALSE]
							
							regr.est_adapted<-MeDeCom:::factorize.regr(D, trueT_adapted)
							trueA_adapted<-regr.est_adapted$A
							rm(regr.est_adapted)
							
							if(is.null(dim(trueA_adapted))){
								trueA_adapted<-matrix(trueA_adapted, nrow=1)
							}
							
							perm2<-MeDeCom:::matchLMCs(fr$T, trueT_adapted, check=FALSE)
							maeA<-MeDeCom:::MAE_A(fr$A, trueA_adapted, perm2)
							
						}else{
							if(!is.null(trueA)){
								maeA<-MeDeCom:::MAE_A(fr$A, trueA, perm)
							}else{
								maeA<-NA
							}
						}
					}else{
						maeA<-NA
					}
				}else{
					rmseT<-maeA<-NA
				}
				
				list(rmseT=rmseT, maeA=maeA, trueA=trueA)
}
#######################################################################################################################
#
# Estimate the cross-validation error
#
#######################################################################################################################
estimateFoldError<-function(Tin, D_out, NFOLDS){
	
	K<-ncol(Tin)
	
	G <- t(Tin) %*% Tin; W <- t(Tin) %*% D_out;
	
	mqsr <- MeDeCom:::RQuadSimplex(G, W, MeDeCom:::randsplxmat(K, ncol(D_out)), 1E-8);
	Anew <- mqsr$A; ftemp<-mqsr$Loss
	
	Dhat <- Tin %*% Anew;
	
	fold.error<-norm(D_out - Dhat, 'F')^2/(ncol(D_out)/NFOLDS);
	
	list(cve=fold.error)
}
#######################################################################################################################
#
# Collect the results of a MeDeCom run
# 
# author Pavlo Lutsik
# author Nina Baumgartner
#  
#######################################################################################################################
collectResults<-function(result_list, cg_subsets, Ks, lambdas, NFOLDS, result_index){
	
	library(methods)
	
	#### Prepare containers
	
	elts<-c("T", "A", "cve", "Fval", "rmse", "rmseT", "maeA")# ,"dist2C")
	
	elt_types<-c(
			"T"=quote(list()), 
			"A"=quote(list()), 
			"cve"=quote(NA_real_), 
			"Fval"=quote(NA_real_), 
			"rmse"=quote(NA_real_), 
			"rmseT"=quote(NA_real_), 
			"maeA"=quote(NA_real_)#, 
			#"dist2C"=quote(NA_real_)
	)
	
	all_results<-list()
	
	all_results$outputs<-list()
	
	## distance to center calculation, require input data 
	if("dist2C" %in% elts){
		source(file.path(result_dir, "analysis_settings.RDump"))
		print(file.path(GLOBAL_DATA_DIR, sprintf("%s_%s_%s", DATASET, DATA_SUBSET, NORMALIZATION), "data.set.RData"))
		load(file.path(GLOBAL_DATA_DIR, sprintf("%s_%s_%s", DATASET, DATA_SUBSET, NORMALIZATION), "data.set.RData"))
		if("SAMPLE_SUBSET" %in% ls()){
			sample_subset<-SAMPLE_SUBSET
		}else{
			sample_subset<-1:ncol(meth.data)
		}
	}

	for(gr in cg_subsets){
		
		results<-list()
		
		folds<-1:NFOLDS
		
		for(elt in 1:length(elts)){
			
			results[[elt]]<-matrix(
					eval(elt_types[[elt]]),
					nrow=length(Ks), ncol=length(lambdas),
					dimnames=list(
							paste("K", Ks, sep="_"), 
							paste("lambda", lambdas, sep="_"))
					)
			### start filling up the matrix
			for(lli in lambdas){
				#finds out the ll_index of lli in lambdas
				ll_index <- NULL
				for (i in 1:length(lambdas)){
					if(lambdas[i]== lli){
						ll_index<- as.numeric(i)
					}
				}			
				for(K in Ks){
					# finds out the index_ks of K in Ks
					K_index<- NULL
					for (i in 1:length(Ks)){
						if(Ks[i]== K){
							K_index<- as.numeric(i)
						}
					}			
					
					if(elts[elt]=="dist2C"){
						ind<-readRDS(file.path(result_dir, sprintf("cg_group_%d.RDS", gr)))
						data_center<-rowMeans(meth.data[ind,sample_subset])
					}		
					if(elts[elt]=="cve"){
						cv.errs<-sapply(folds, function(fold){	
										res_idx<-result_index[
												#match(run, c("initial", "cv")),
												gr,
												match(fold, c(0,folds)),
												K_index,
												ll_index]
										res<-result_list[[res_idx]]
										return(el(res, where=3))
									
								})
						
						results[[elt]][K_index, ll_index]<-mean(cv.errs)							
					}else if (elts[elt]=="dist2C"){
						res<-NULL
							res_idx<-result_index[
									gr,
									1,
									K_index,
									ll_index]
							res<-result_list[[res_idx]]
						if(!is.null(res)){
							That<-el(res, where="That")
							
							dists<-apply(That, 2, function(vr) sum((vr-data_center)^2))
							results[[elt]][[K_index, ll_index]]<-mean(dists)
						}
					}else{
						res<-NULL
							res_idx<-result_index[
									gr,
									1,
									K_index,
									ll_index]
							res<-result_list[[res_idx]]
						if(!is.null(res)){
							if(is.list(eval(elt_types[[elt]]))){
								results[[elt]][[K_index, ll_index]]<-el(res, where=elt)
							}else{
								results[[elt]][K_index, ll_index]<-el(res, where=elt)
							}
						}
					}
					
				}
			}
		}
		
		names(results)<-elts
		
		all_results$outputs[[as.character(gr)]]<-results
	}
	
	help<- NULL
	for ( i in cg_subsets){
		help <- append(help, cg_subsets[[i]])
	}
	all_results$parameters$cg_subsets<-help
	all_results$parameters$Ks<-Ks
	all_results$parameters$lambdas<-lambdas
	
	####
	# add parameters for description in shiny
	####
	
	#parameters we already know
	#all_results$parameters$LAMBDA_GRID<-sort(lambdas)
	#all_results$parameters$Ks<-sort(Ks)
	
	#parameters from analysis_settings.RDump
	if(FALSE){
		fn<-file.path(result_dir, "analysis_settings.RDump")
		if (file.exists(fn)){
			#read file analysis_settings.RDump
			analysis_settings<- source(fn)
			
			#analysis name
			all_results$parameters$ANALYSIS <- ANALYSIS
			#groups
			help<- NULL
			if (length(groups )==1){	
				help <- append(help, groups[[1]][i])
			}else{
				for ( i in groups){
					print("das ist group list")
					print (groups[i])
					help <- append(help, groups[[i]])
				}
			}
			all_results$parameters$groups<-help
			#normalization                                                  
			all_results$parameters$NORMALIZATION <- NORMALIZATION    
			#maximal number of iterations                                   
			all_results$parameters$ITERMAX<-ITERMAX                     
			
			#marker_selection                                                                            
			all_results$parameters$MARKER_SELECTION<-MARKER_SELECTION           
			#Folds                                                                                       
			all_results$parameters$NFOLDS<- NFOLDS           
			#addition information                                                                        
			all_results$parameters$ANALYSIS_TOKEN<-ANALYSIS_TOKEN 
			#number of random initialization                                                    
			all_results$parameters$NINIT<-NINIT
			#Data set                                                               
			all_results$parameters$DATASET<-DATASET 
			#Data subset                                                               
			all_results$parameters$DATA_SUBSET<- DATA_SUBSET 
			#original group_list                                                       
			help<- NULL
			# necessary , otherwise it is a nested list [[]]
			for ( i in 1:length(groups)){
				help <- append(help, groups[[i]])
			}
			all_results$parameters$ORIGINAL_groups<- help
			
		}else{
			stop("Analysis_settings.RDump don't exist!!!")
		}
	}
	
	if(FALSE){
		print(all_results)
		# save the final R object
		saveRDS(all_results, file=sprintf("collected.results.RDS"))
	}
	
	return(all_results)
}
#######################################################################################################################