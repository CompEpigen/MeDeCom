#######################################################################################################################
# plotting.R
#
# Plotting functions
#
# author: Pavlo Lutsik
#######################################################################################################################


#######################################################################################################################
###   PLOTTING UTILITIES
#######################################################################################################################
openPSfile<-function(filename, w, h, ...){
  
  setEPS()
  postscript(file=filename, horizontal=F,
             onefile=F,
             width=w, height=h,
             family=c("/home/lutsik/R/fonts/arial.afm",
                      "/home/lutsik/R/fonts/arialbd.afm",
                      "/home/lutsik/R/fonts/ariali.afm",
                      "/home/lutsik/R/fonts/arialbi.afm"),
             pointsize=12, paper="special", ...)
  
  
}
#######################################################################################################################
openPDFfile<-function(filename, w, h, ...){
	
#	setEPS()
	pdf(file=filename,
			onefile=F,
			width=w, height=h,
			family=c("/home/lutsik/R/fonts/arial.afm",
					"/home/lutsik/R/fonts/arialbd.afm",
					"/home/lutsik/R/fonts/ariali.afm",
					"/home/lutsik/R/fonts/arialbi.afm"),
			pointsize=12, paper="special", ...)
	
	
}
#######################################################################################################################
drawMarkerBarplot<-function(data, cg)
{
  high<-grep("_H", colnames(data))
  low<-grep("_L", colnames(data))
  resp<-as.integer(rbind(1:(dim(data)[2]/2),1:(dim(data)[2]/2)))
  resp[low]<-(-resp[low])
  if(cg %in% rownames(data))
    barplot(data[cg,order(2*abs(resp)+sign(resp))],beside=T, col=c("blue","red"), las=2, ylim=c(0,1), main=cg)
  
}
#######################################################################################################################
pcaVarBarplot<-function(pca){
  
  heights<-as.matrix((pca$sdev^2)/sum(pca$sdev^2))
  rownames(heights)<-paste("PC", 1:length(heights), sep="")
  barplot(heights, beside=F, cex.names=0.75, las=2, ylab="fraction of variance",border=T, legend.text=F, col=grey.colors(34, gamma=5))
  
}
#######################################################################################################################
pcaProjPlot<-function(pca,comps=c(1,2), labs=T, ...){
  
 
  plot(as.formula(paste(paste("PC", comps[2:1], sep=""), collapse="~")),
       pca$x[,comps],
       ...)
  
  if(labs) text(pca$x[,comps], labels=rownames(pca$x), pos=4, offset=0.5, ...)
  
  
}
#######################################################################################################################
densityScatterPlot<-function(data, columns, legend.only=F, ...){
  
  colors<-densCols(data[,columns])
#   par(mai=c(0.5,0.5,0.05,0.05), mar=c(3,3,1,1))
   par(mar=c(2,2,1,1), mgp=c(1,0,0))
#    par(oma=c(1,1,1,1))
  
#   si<-sample.int(dim(data)[1], min(dim(data)[1],50000), replace=F)
  
  
  if(!legend.only) {
    bins<-matrix(nrow=100, ncol=100)
    
    sapply(1:dim(data)[1], function(row) {
      x=data[row, columns[1]]
      y=data[row, columns[2]]
      bin.x=ceiling((1-x)*100)
      bin.y=ceiling((1-y)*100)
      if(is.na(bins[bin.x, bin.y])) bins[bin.x, bin.y]<<-row
    })
    si<-na.omit(as.numeric(bins))
    print(paste("Plotting ",length(si), " points"))
    
    
    plot(data[si,columns], col=colors[si], pch=".", cex.axis=0.8, cex.lab=0.8, tcl=0.2)
#   smoothScatter(data[,columns])
  text(0.2, 0.9, paste("r=", sprintf("%1.3f", cor(data[,columns[1]], data[,columns[2]]))),cex=0.75)
  }else{
    counts<-binXYdata2D(data[,columns[1]],data[,columns[2]], nbins=128)
  
    bins<-matrix(nrow=128, ncol=128)
    
    sapply(1:dim(data)[1], function(row) {
      x=data[row, columns[1]]
      y=data[row, columns[2]]
      bin.x=ceiling((1-x)*100)
      bin.y=ceiling((1-y)*100)
      if(is.na(bins[bin.x, bin.y])) bins[bin.x, bin.y]<<-row
    })
    si<-na.omit(as.numeric(bins))
    par(mai=c(0,0,0,0))
    plot.new()
    image.plot(legend.only=T, col=c("#FFFFFF", sort(unique(colors), decreasing=T)), zlim=c(0,max(counts)), legend.shrink = 0.5,
               graphics.reset=T, ...)
    
  }
}
#######################################################################################################################
data.lineplot<-function(data, quant.trait=NULL, 
		platform="probes450", 
		annot=NULL, 
		col="cgs",
		pts=col,
		lty=col,
		summary=FALSE, 
		xlabel="quantitative trait", 
		ylabel="DNA methylation", 
		legend=TRUE, 
		title=NULL,
		plot=TRUE){
	
	if(!require(ggplot2)){
		stop("This functionality depends on the ggplot2 package. Please install it and repeat again.")
	}
	
	if(!is.matrix(data) || !is.numeric(data)){
		stop("invalid argument for data")
	}
	
	if(is.null(quant.trait)){
		quant.trait<-c(0, (1:(ncol(data)-2))/(ncol(data)-1), 1)
	}
	
	if(is.null(annot)){
		if(!is.null(platform)){
			annot<-rnb.annotation2data.frame(rnb.get.annotation(platform), add.names=TRUE)[rownames(data),]
		}else{
			annot<-NULL
		}
	}else{
		annot<-annot[rownames(data),]
	}
			
	line.data<-data.frame(data=as.numeric(t(data)), 
			cgs=as.character(sapply(rownames(data), rep, times=ncol(data))),
			quant.trait=as.numeric(t(sapply(sort(quant.trait), rep, times=nrow(data))))
			)
	if(length(col)==nrow(data)){
		line.data$Group<-as.character(sapply(col, rep, times=ncol(data)))
	}
	if(length(pts)==nrow(data)){
		line.data$Group2=as.character(sapply(pts, rep, times=ncol(data)))
	}
	
	
	if(!is.null(annot)){
		for(n in colnames(annot)){
			line.data[[n]]=as.character(sapply(annot[[n]], rep, times=ncol(data)))
		}
	}
	
	lp<-ggplot(line.data, aes(x=quant.trait, y=data, group=cgs)) + xlab(xlabel) + ylab(ylabel)
	
	if(is.character(col) && length(col)==1){
		lp<-lp + geom_line(aes_string(color=col, linetype=col), size=0.5)
	}else{
		lp<-lp + geom_line(aes_string(color="Group", linetype="Group2"), size=0.5)
	}
	
	if(is.character(pts) && length(pts)==1){
		lp<-lp + geom_point(aes_string(color=col, shape=pts))
	}else{
		lp<-lp + geom_point(aes_string(color="Group", shape="Group2"))
	}
	
	
	if(summary) lp<-lp + stat_summary(aes(group=NA), fun.y='median', geom='line', size=2)
	
	if(!legend) lp<-lp + theme(legend.position = "none")
	
	if(!is.null(title)) lp<-lp + ggtitle(title)
	
	if(plot)
		print(lp) else lp
}
#######################################################################################################################
# Heatmap with an x-axis on top
# 
# Author: Pavlo Lutsik 
#######################################################################################################################
heatmap.mod<-function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
		distfun = dist, hclustfun = hclust, reorderfun = function(d, 
				w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
				"Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
		margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
				1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
		labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
		verbose = getOption("verbose"), 
		xaxis=1, xaxt="nn", ylab.side=4,...) 
{
	scale <- if (symm && missing(scale)) 
				"none"
			else match.arg(scale)
	if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
		stop("'x' must be a numeric matrix")
	nr <- di[1L]
	nc <- di[2L]
	if (nr <= 1 || nc <= 1) 
		stop("'x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2L) 
		stop("'margins' must be a numeric vector of length 2")
	doRdend <- !identical(Rowv, NA)
	doCdend <- !identical(Colv, NA)
	if (!doRdend && identical(Colv, "Rowv")) 
		doCdend <- FALSE
	if (is.null(Rowv)) 
		Rowv <- rowMeans(x, na.rm = na.rm)
	if (is.null(Colv)) 
		Colv <- colMeans(x, na.rm = na.rm)
	if (doRdend) {
		if (inherits(Rowv, "dendrogram")) 
			ddr <- Rowv
		else {
			hcr <- hclustfun(distfun(x))
			ddr <- as.dendrogram(hcr)
			if (!is.logical(Rowv) || Rowv) 
				ddr <- reorderfun(ddr, Rowv)
		}
		if (nr != length(rowInd <- order.dendrogram(ddr))) 
			stop("row dendrogram ordering gave index of wrong length")
	}
	else rowInd <- 1L:nr
	if (doCdend) {
		if (inherits(Colv, "dendrogram")) 
			ddc <- Colv
		else if (identical(Colv, "Rowv")) {
			if (nr != nc) 
				stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
			ddc <- ddr
		}
		else {
			hcc <- hclustfun(distfun(if (symm) 
										x
									else t(x)))
			ddc <- as.dendrogram(hcc)
			if (!is.logical(Colv) || Colv) 
				ddc <- reorderfun(ddc, Colv)
		}
		if (nc != length(colInd <- order.dendrogram(ddc))) 
			stop("column dendrogram ordering gave index of wrong length")
	}
	else colInd <- 1L:nc
	x <- x[rowInd, colInd]
	labRow <- if (is.null(labRow)) 
				if (is.null(rownames(x))) 
					(1L:nr)[rowInd]
				else rownames(x)
			else labRow[rowInd]
	labCol <- if (is.null(labCol)) 
				if (is.null(colnames(x))) 
					(1L:nc)[colInd]
				else colnames(x)
			else labCol[colInd]
	if (scale == "row") {
		x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
		sx <- apply(x, 1L, sd, na.rm = na.rm)
		x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
	}
	else if (scale == "column") {
		x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
		sx <- apply(x, 2L, sd, na.rm = na.rm)
		x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
	}
	lmat <- rbind(c(NA, 3), 2:1)
	lwid <- c(if (doRdend) 1 else 0.05, 4)
	lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
			4)
	if (!missing(ColSideColors)) {
		if (!is.character(ColSideColors) || length(ColSideColors) != 
				nc) 
			stop("'ColSideColors' must be a character vector of length ncol(x)")
		lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
		lhei <- c(lhei[1L], 0.2, lhei[2L])
	}
	if (!missing(RowSideColors)) {
		if (!is.character(RowSideColors) || length(RowSideColors) != 
				nr) 
			stop("'RowSideColors' must be a character vector of length nrow(x)")
		lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
						1), lmat[, 2] + 1)
		lwid <- c(lwid[1L], 0.2, lwid[2L])
	}
	lmat[is.na(lmat)] <- 0
	if (verbose) {
		cat("layout: widths = ", lwid, ", heights = ", lhei, 
				"; lmat=\n")
		print(lmat)
	}
	on.exit(dev.flush())
	op <- par(no.readonly = TRUE)
	on.exit(par(op), add = TRUE)
	layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
	if (!missing(RowSideColors)) {
		par(mar = c(margins[1L], 0, 0, 0.5))
		image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
	}
	if (!missing(ColSideColors)) {
		par(mar = c(0.5, 0, 0, margins[2L]))
		image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
	}
	if (xaxis==3){
		par(mar = c(margins[1L],margins[2L], margins[1L], 0))
	}else{
		par(mar = c(margins[1L], 0, 0, margins[2L]))  
	}
	
	
	if (!symm || scale != "none") 
		x <- t(x)
	if (revC) {
		iy <- nr:1
		if (doRdend) 
			ddr <- rev(ddr)
		x <- x[, iy]
	}
	else iy <- 1L:nr
	image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
	if(xaxt!="n") axis(xaxis, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
				cex.axis = cexCol)
	if (!is.null(xlab)) 
		mtext(xlab, side = 1, line = margins[1L] - 1.25)
	axis(2, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
			cex.axis = cexRow)
	
	if (!is.null(ylab)) 
		mtext(ylab, side = ylab.side, line = margins[2L] - 1.25)
	if (!missing(add.expr)) 
		eval(substitute(add.expr))
	par(mar = c(margins[1L], 0, 0, 0))
	if (doRdend) 
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	else frame()
	if(xaxis==3) par(mar = c(0, margins[2L], if (!is.null(main)) 1 else 0, 0))
	else par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
	if (doCdend) 
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	else if (!is.null(main)) 
		frame()
	if (!is.null(main)) {
		par(xpd = NA)
		title(main, cex.main = 1.5 * op[["cex.main"]])
	}
	invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
							doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}
#######################################################################################################################
###   VISUALIZING THE RESULTS OF DECOMPOSITION
#######################################################################################################################
plot.errs<-function(
		applyAllList, 
		data, 
		trueT=NULL, 
		trueA=NULL,
		plots=c("RMSE",
				if(is.null(trueT)) character(0) else "RMSE.T",
				if(is.null(trueA)) character(0) else "MAE.A"
		)){
	
	rmse.D<-sapply(applyAllList, function(fr){
				
				sqrt(sum((data-fr$T%*%fr$A)^2)/ncol(data)/nrow(data))
				
			})
	
	rmse.T<-NULL
	if(!is.null(trueT)){
		
		rmse.T<-sapply(applyAllList, function(fr){
					
					perm<-matchLMCs(fr$T, trueT)
					if(!is.null(perm)){
						sqrt(sum((trueT-fr$T[,perm])^2)/ncol(data)/nrow(data))
					}else{
						NA
					}
				})
	}
	
	mae.A<-NULL
	if(!is.null(trueA)){
		
		mae.A<-sapply(applyAllList, function(fr){
					if(!is.null(trueT)){
						perm<-matchLMCs(fr$T, trueT)
					}else{
						perm<-1:nrow(fr$A)
					}
					if(!is.null(perm)){
						sum(abs(trueA-fr$A[perm,]))/ncol(trueA)/nrow(trueA)
					}else{
						NA
					}
					
				})
	}
	
	layout(matrix(1:length(plots),nrow=1, ncol=length(plots)))
	
	if("RMSE" %in% plots){	
		barplot(rmse.D, las=2, ylab="RMSE")
	}
	
	if(!is.null(rmse.T) && "RMSE.T" %in% plots){
		barplot(rmse.T, las=2, ylab="RMSE, T")
	}
	
	if(!is.null(mae.A) && "MAE.A" %in% plots){
		barplot(mae.A, las=2, ylab="MAE, A")
	}
	
	return(list(rmseD=rmse.D, rmseT=rmse.T, maeA=mae.A))
	
}
#######################################################################################################################
plot.cves<-function(
		methods=ALGORITHMS[-c(1,5,6)],
		ymax=1,
		lambdas=NULL,
		lambdamax=1,
		log.scale=TRUE, 
		leg.pos="bottomleft", 
		...
){
	
	if(is.null(lambdas)){
		lambdas<-RELGRID*lambdamax
	}
	
	#plot(NA, ylim=c(0,ymax), xlim=c(0,max(lambdas)))
	
	cves<-list()
	for(method in methods){
		cves[[method]]<-MeDeCom:::test.lambdagrid(...,method=method, lambdagrid=lambdas)
		
	}
	
	method<-names(cves[1])
	if(log.scale){
		plot(log10(lambdas), cves[[1]]$cve, 
				pch=ALGORITHM.PCH[method],
				col=ALGORITHM.COLS[method],
				lty=2, type="o",
				ylim=c(0,ymax), xlim=c(min(log10(lambdas)),max(log10(lambdas))))
	}else{
		plot(lambdas, cves[[1]]$cve,
				pch=ALGORITHM.PCH[method],
				col=ALGORITHM.COLS[method], 
				lty=2, type="o",
				ylim=c(0,ymax), xlim=c(0,max(lambdas)))
	}
	dump<-lapply(methods[-1], function(method){
				if(log.scale){
					lines(log10(lambdas), cves[[method]]$cve, 
							pch=ALGORITHM.PCH[method],
							col=ALGORITHM.COLS[method], lty=2, type="o")
				}else{
					lines(lambdas, cves[[method]]$cve,
							pch=ALGORITHM.PCH[method],
							col=ALGORITHM.COLS[method], lty=2, type="o")
				}
				
			})
	
	legend(leg.pos, legend=methods,
			col=MeDeCom:::ALGORITHM.COLS[methods],
			pch=MeDeCom:::ALGORITHM.PCH[methods],
			lty=2)

	return(invisible(list(CVEs=cves, lambdas=lambdas)))
}
#######################################################################################################################
components.number.plot<-function(
		data, 
		trueT=NULL, 
		maxk=2, 
		V=c(0,1),
		methods=names(ALGORITHM.COLS),
		opt.exact=opt.factorize.exact(affine=FALSE), 
		max.rmse=1,
		leg.pos="topright",
		main="",
		data.complement=NULL,
		lambda=0,
		remove.regression=FALSE,
		ktest=NULL,
		...
){
	
	if(is.null(ktest)){
		if(missing(data)){
			stop("Neither the data nor the results of the previous computation have been supplied")
		}
		ktest<-sapply(1:maxk, function(kk){
					all<-MeDeCom:::applyAll(data, trueT, k=kk, v=V,
							opt.exact=opt.exact, methods=methods,
							lambda=lambda,...)
					all
				})
	}
	ktest<-ktest[methods,,drop=F]
	
	plot(NA, ylim=c(0,max.rmse), xlim=c(1,maxk), xlab="r", ylab="RMSE", main=main)
	meths<-rownames(ktest)
	if("regression" %in% meths && remove.regression)
		meths<-setdiff(meths, "regression")
	dump<-sapply(meths, function(meth){
				errs<-sapply(ktest[meth,], el, where="rmse")
				lines(1:maxk, errs, type="o", lty=2,
						pch=MeDeCom:::ALGORITHM.PCH[meth],
						col=MeDeCom:::ALGORITHM.COLS[meth])
			})
	legend(leg.pos, legend=meths,
			col=MeDeCom:::ALGORITHM.COLS[meths],
			pch=MeDeCom:::ALGORITHM.PCH[meths],lty=2)
	
	
	if(!is.null(data.complement)){
		
		plot(NA, ylim=c(0,max.rmse), xlim=c(1,maxk), xlab="r", ylab="RMSE", main=main)
		
		dump<-sapply(meths, function(meth){
					errs<-sapply(ktest[meth,], function(result){
								
								if(meth %in% names(ALGORITHM.COLS)[c(1,2)]){
									Tcompl=getT.hlasso(data.complement, result$A)
								}else if(meth %in% names(ALGORITHM.COLS)[c(4)]){
									Tcompl=getT.empirical(data.complement, result$A,
											V=as.numeric(trueT))						
								}else{
									Tcompl=getT.intfac(data.complement, result$A, 
											V=V)	
								}
								err<-sqrt((sum((data.complement-Tcompl%*%result$A)^2))/
												ncol(result$A)/nrow(Tcompl))
								err
							})
					
					lines(1:maxk, errs, type="o", lty=2, pch=15,
							col=MeDeCom:::ALGORITHM.COLS[meth])
				})
		legend(leg.pos, legend=meths,
				col=MeDeCom:::ALGORITHM.COLS[meths],
				pch=MeDeCom:::ALGORITHM.PCH[meths],
				lty=2)
	}
	return(invisible(ktest))
}
#######################################################################################################################
plot.K.selection<-function(
		MeDeComSet, 
		statistic="cve", 
		Ks=integer(), 
		lambdas=numeric(),
		cg_subset=1,
		D=NULL, 
		cg_subsets=NULL, 
		sample_subset=NULL,
		KvsRMSElambdaLegend=TRUE,
		normalizedCVE=FALSE,
		addPlotTitle=FALSE
){

	all_results<-MeDeComSet
	#gr<-as.integer(input$cg_group)
	gr<-cg_subset
	
	if(length(Ks)==0){
		Ks<-MeDeComSet@parameters$Ks
	}
	
	if(length(lambdas)<1){
		lambdas<-MeDeComSet@parameters$lambdas
	}
	
  if(statistic=="rmse" && !is.null(D)){
  	meth.data<-D
  
  	#if("SAMPLE_SUBSET"%in% names(getRuns()[[input$analysisrun]])){
  	#	sample_subset<-getRuns()[[input$analysisrun]][["SAMPLE_SUBSET"]]
  	if(is.null(sample.subset)){
  		sample_subset<-1:ncol(meth.data)
  	}
  
  	if(is.null(cg_subsets)){
  		cg_subsets<-list(c(1:nrow(D)))
  	}
  
  #		ind<-readRDS(sprintf("%s/cg_group_%d.RDS",
  #						getRuns()[[input$analysisrun]][["run.dir"]],
  #						#dataset()$groups[gr]
  #						gr
  #				))
  
  	startRMSE<-sqrt(mean((meth.data[ind,sample_subset]-rowMeans(meth.data[ind,sample_subset]))^2))
  }else{
  	startRMSE<-NA
  }
  
	plotTitle<-""

	if(addPlotTitle){
		plotTitle<-sprintf("GG group %d", cg_subset)
	}

	allRMSE<-getStatistics(MeDeComSet, Ks, lambdas, cg_subset, statistic=statistic)
	if(!is.null(dim(allRMSE))){
		allRMSE<-matrix(allRMSE, nrow=length(Ks), ncol=length(lambdas))
	}
	if(statistic=="CVE" && normalizedCVE){
		allRMSE<-sqrt(allRMSE/length(getCGsubset()))
	}
	
	if(statistic=="RMSE" && !is.na(startRMSE)){
		ymaxv=max(na.omit(c(as.numeric(allRMSE), startRMSE)))
		yminv=min(na.omit(c(as.numeric(allRMSE), startRMSE)))
	}else if(statistic %in% c("cve","rmse") &&  normalizedCVE){
		ymaxv=1.0
		yminv=0.0
	}else{
		ymaxv=max(na.omit(as.numeric(allRMSE)))
		yminv=min(na.omit(as.numeric(allRMSE)))
	}
	
	if(KvsRMSElambdaLegend){
		layout(matrix(1:2, ncol=2),widths=c(0.75, 0.25))
		par(mar=c(4,4,2,2))
	}
	if(all(is.na(allRMSE))){
	  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
	  text(x=0.5,y=0.5, paste("Statistic:", statistic,"not available for this result"))
	  return(NULL)
	}
	matplot(Ks, 
			allRMSE,
			type="o", 
			pch=15,
			xaxt = "n",
			ylim=c(yminv-0.1*(ymaxv-yminv), 
					ymaxv+0.1*(ymaxv-yminv)), 
			col=1:length(all_results@parameters$lambdas),
			lty=1:length(all_results@parameters$lambdas),
			xlab="k",
			ylab=names(PERFORMANCE_MEASURES)[which(PERFORMANCE_MEASURES==statistic)],
			main=plotTitle
	)
	axis(1, at=Ks, labels=as.character(Ks))
	
	if(statistic=="rmse"){
		abline(h=startRMSE)
	}
	if(KvsRMSElambdaLegend){
		plot.new()
		legend("right", legend=as.character(all_results@parameters$lambdas), 
				#col=lambdaCols[1:length(all_results$lambdas)], 
				col=1:length(all_results@parameters$lambdas),
				lty=1:length(all_results@parameters$lambdas),
				pch=15,
				xpd=TRUE,
				title="lambda")
	}
}
#######################################################################################################################
plot.lambda.selection<-function(
		MeDeComSet,
		K,
		cg_subset=1,
		minLambda=0,
		maxLambda=Inf,
		scale="native",
		includeRMSE_T=FALSE,
		includeMAE_A=FALSE,
		includeDist2C_T=FALSE
		){
	
	results<-MeDeComSet
	llab <- "lambda"
	llsubs<-results@parameters$lambdas >= as.numeric(minLambda) & results@parameters$lambdas <= as.numeric(maxLambda)
	
	elts<-PERFORMANCE_MEASURES
	
	if(!includeRMSE_T){
		elts<-elts[-grep("rmseT", elts)]
	}
	if(!includeMAE_A){
		elts<-elts[-grep("maeA", elts)]
	}
	if(!includeDist2C_T){
		elts<-elts[-grep("dist2C", elts)]
	}
	
	layout(matrix(1:length(elts), ncol=1))
	par(oma=c(5,3,5,3))
	par(mar=c(3,8,3,0))
	lmbd<-results@parameters$lambdas[llsubs]
	
	for(elt in 1:length(elts)){
		
		vals <- as.numeric(getStatistics(MeDeComSet, K, lmbd, cg_subset, statistic=elts[elt]))
		offs<-(max(vals[!is.na(vals)])-min(vals[!is.na(vals)]))/5
		
		if(!all(is.na(vals))){
			if(scale=="log"){
				plotting_subset<-!is.na(vals) & lmbd>0
			}else{
				plotting_subset<-!is.na(vals)
			}
			plot_order<-order(lmbd[plotting_subset], decreasing=FALSE)
			if(scale=="log"){
				plot(lmbd[plotting_subset][plot_order], 
						vals[plotting_subset][plot_order],
						type="o", pch=25, 
						axes=F,
						log="x",
						cex=2, ps=5, xpd=TRUE,
						xlab=NA, ylab="", col="blue", bg="blue")
			}else{
				plot(lmbd[plotting_subset][plot_order], 
						vals[plotting_subset][plot_order],
						type="o", pch=25, 
						axes=F,
						cex=2, ps=5, xpd=TRUE,
						xlab=NA, ylab="", col="blue", bg="blue")
			}
			if(elt==length(elts)){
				if(scale=="native"){
					par(xaxs="r")
				}
				axis(side = 1, cex.axis=1.5, xpd=TRUE, line=3)
			}
			ticks<-seq(max(0,min(vals[!is.na(vals)])-offs), max(vals[!is.na(vals)])+offs, by=offs)
			axs<-axis(side = 2, 
					las = 1, 
					cex.axis=1.5, 
					xpd=TRUE,
					at=ticks,
					labels=sprintf("%#1.5g", ticks)
			)
			mtext(names(elts)[elt],
					side=2, 
					line=9, 
					cex=1.2)
			abline(h=axs, lty=3)
			abline(v=lmbd[!is.na(vals)], lty=3, 
					xpd=TRUE
			)
			if(scale=="log" && any(lmbd[!is.na(vals)]==0)){
				abline(h=vals[!is.na(vals)][which(lmbd[!is.na(vals)]==0)], lty=1)
				text(0.9*max(lmbd[!is.na(vals) & lmbd>0]), 
						vals[!is.na(vals)][which(lmbd[!is.na(vals)]==0)] + 0.1*(max(vals[!is.na(vals)])-min(vals[!is.na(vals)])), 
						expression(lambda == paste(0)), cex=1.5)
			}
		}
	}
	mtext(
			expression(lambda), 
			side=1, outer=TRUE, 
			line=4,
			cex=1.2, 
	)
	par(mfrow=c(1,1))
}
#######################################################################################################################
#'
#' plotParameters
#' 
#' Parameter selection plots 
#' 
#' @param MeDeComSet	MeDeCom object
#' @param cg_subset		integer index of the CpG subset (if no subsets were used)
#' @param Ks			values of parameter k to use (by default all k values available in the \code{MeDeComSet} are used)
#' @param lambdas		values of parameter lambda to use (by default all lambda values available in the \code{MeDeComSet} are used)
#' @param statistic		if multiple k values are supplied, the statistic which is plotted (defaults to cross-validation error)
#' @param minLambda		minimal lambda value
#' @param maxLambda		maximal lambda value
#' @param lambdaScale character indicating if native scale or logarithmic scale should be employed for plotting lambda
#' @param ...  further paramters passed to \code{plot.lambda.selection} or \code{plot.K.selection}
#' 
#' @export
plotParameters<-function(
		MeDeComSet,
		cg_subset=1,
		Ks=integer(),
		lambdas=integer(),
		statistic="cve",
		minLambda=0,
		maxLambda=Inf,
		lambdaScale="native",
		...
){
	check_inputs(MeDeComSet, cg_subset, Ks, lambdas)
	
	if(length(Ks)==1){
		plot.lambda.selection(MeDeComSet,
						Ks,
						cg_subset,
						minLambda,
						maxLambda,
						scale=lambdaScale,
						...)
	}else{
		plot.K.selection(MeDeComSet,
				statistic,
				Ks,
				lambdas,
				cg_subset,
				...)
	}
}
#######################################################################################################################
component.scatter<-function(
		TT, 
		Tref, 
		lmc=1, 
		match=NA, 
		Ahat=NULL, 
		Aref=NULL,
		highlight=0, 
		method="pearson", 
		ncols=2, 
		nrows=2, 
		point.cat=character(), 
		cg.feature=NULL,
		smooth=FALSE
){
	
	layout(t(matrix(1:(ncols*nrows), ncol=nrows, nrow=ncols)))
	ind<-1:nrow(TT)
	if(!smooth && !is.null(cg.feature)){
		if(is.factor(cg.feature) || is.character(cg.feature)){
			cats<-unique(cg.feature)
			catv<-as.integer(cg.feature)
			catv[is.na(cg.feature)]<-0L
			col.legend<-c("grey", rainbow(length(cats)))
			point.cols<-col.legend[catv+1L]
			text.legend<-as.character(cats)
		}else if(is.numeric(cg.feature)){
			palette<-grey.colors(start=0.9, end=0.1, min(20, length(cg.feature)))
			col.legend<-palette
			point.cols<-palette[cut(cg.feature,min(20, length(cg.feature)),labels=FALSE)]
			point.subs<-which(cg.feature>input$quantFilterMin & cg.feature<input$quantFilterMax)
			text.legend<-levels(cut(cg.feature,min(20, length(cg.feature))))
		}else{
			point.cols<-rep("black", length(ind))
		}
	}else{
		point.cols<-rep("black", length(ind))
	}

	if(!is.na(match) && length(match)!=ncol(TT)){
		stop("Wrong value for match: if not NA the length should be equal to ncol(TT)")
	}

	if(ncol(TT)>ncols*nrows){
		stop("too little plotting windows")
	}
	
	if(is.na(match)){
		that_cols<-lmc
		tref_cols<-1:ncol(Tref)
	}else{
		if(length(lmc)>0){
			that_cols<-lapply(1:ncol(Tref), function(trefc) which(match==trefc))
			tref_cols<-1:ncol(Tref)
		}else{
			that_cols<-1:ncol(That)
			tref_cols<-match[lmc]
		}
	}
	
	for(lmci in that_cols){	
		for(i in tref_cols){
			if(length(lmci)>0){
				if(method=="anova"){
					if(is.null(Ahat)){
						comp.data<-jitter(rowMeans(TT[,lmci,drop=FALSE]))
					}else{
						comp.data<-jitter(rowMeans(TT[,lmci,drop=FALSE] %*% sweep(Ahat[lmci,,drop=FALSE], 2, colSums(Ahat[lmci,,drop=FALSE])+10^(-8), "/"))) 
					}
					fit<-aov(ct~comp1,data=data.frame(ct=Tref[,i], comp1=rowMeans(TT[,lmci])))
					header<-sprintf("ANOVA F-value: %2.1f",
							round(summary(fit)[[1]]$`F value`[1]))
				}else if(method=="pearson"){
					if(is.null(Ahat)){
						comp.data<-rowMeans(TT[,lmci,drop=FALSE])
					}else{
						comp.data<-rowMeans(TT[,lmci,drop=FALSE] %*% sweep(Ahat[lmci,,drop=FALSE], 2, colSums(Ahat[lmci,,drop=FALSE])+10^(-8), "/")) 
					}
					header<-sprintf("Pearson correlation: %1.3f",
							cor(Tref[,i],rowMeans(TT[,lmci,drop=FALSE]),
									use="pairwise.complete.obs"))
				}
				
				if(length(unique(Tref[,i]))<10){
					ref.data<-jitter(Tref[,i])
				}else{
					ref.data<-Tref[,i]
				}
				
				xl=sprintf("average of LMC %s", paste(lmci, collapse=", "))
				yl=colnames(Tref)[i]
				if(smooth && method!="anova"){
					smoothScatter(comp.data, ref.data,
							colramp = colorRampPalette(c("white", "blue","yellow",  "red", "black")),
							xlab=xl,
							ylab=yl,
							main=header)
				}else{
					plot(comp.data, ref.data, cex=0.5, 
							xlab=xl,
							ylab=yl,
							col=point.cols,
							main=header)
				}
				if(i %in% highlight){
					box(which="figure", lwd=2, col="red")
				}
			}
		}
	}
	
	if(!smooth && !is.null(cg.feature)){
		plot(1, type="n", axes=FALSE, xlab="", ylab="")
		legend("center", text.legend,  col=col.legend, pch=20)
	}
}
components.scatter<-component.scatter
#######################################################################################################################
component.heatmap<-function(
		TT, 
		Tref=NULL, 
		method="pearson", 
		color.int=1000, 
		centered=FALSE, 
		...
){
	
	if(is.null(Tref)){
		Tref<-TT
	}
	
	if(!method %in% c('pearson', 'spearman', 'euclidean', 'angular')){
		stop("Not implemented yet")
	}
	
	if(centered){
		th<-sweep(TT, 1, rowMeans(TT), "-")
		if(!is.null(Tref)){
			tr<-sweep(Tref, 1, rowMeans(Tref), "-")
		}else{
			tr<-th
		}
	}else{
		th<-TT
		if(!is.null(Tref)){
			tr<-Tref
		}else{
			tr<-th
		}
	}
	
	if(method %in% c('pearson', 'spearman')){
		data<-cor(th, tr, use='pairwise.complete.obs', method=method)
	}else if (method %in% c("euclidean","angular")){
		data<-matrix(ncol=ncol(tr), nrow=ncol(th))
		for(ri in 1:ncol(th)){
			for(ci in 1:ncol(tr)){
				if(method=="euclidean"){
					data[ri,ci]<-sqrt(sum((tr[,ci]-th[,ri])^2))
				}else if(method=="angular"){
					data[ri,ci]<-sum(tr[,ci]*th[,ri])/sqrt(sum(tr[,ci]^2))/sqrt(sum(th[,ri]^2))
				}
			}
		}
	}
	
	notes<-matrix(NA_character_, nrow(data), ncol(data))
#	for(cn in 1:ncol(data)){
#		notes[which.max(data[,cn]),cn]<-sprintf("%1.3f", max(data[,cn]))
#	}
	
	for(cn in 1:nrow(data)){
		notes[cn,which.max(data[cn,])]<-sprintf("%1.3f", max(data[cn,]))
	}
	
	for(cn in 1:ncol(data)){
		if(is.na(notes[which.max(data[,cn]),cn])){
			notes[which.max(data[,cn]),cn]<-"*"
		}else{
			notes[which.max(data[,cn]),cn]<-paste0(notes[which.max(data[,cn]),cn], "*")
		}
	}
	
	#colfunc<- function(color.int) brewer.pal(11,"RdBu")[color.int]
	
	if(method %in% c("pearson", "spearman","angular")){
		colfunc <- colorRampPalette(c("brown", "red", "pink","white", "grey", "blue", "darkblue"))
		col_vector<-colfunc(color.int)[floor(color.int*(min(data)/2+0.5)):ceiling(color.int*(max(data)/2+0.5))]
	}else{
		col_vector<-heat.colors(max(100,color.int%/%10))
	}
	key_title<-sprintf("%s correlation coef.", method) 
	heatmap.2(t(data), trace="none", scale="none",
			col=col_vector, 
			dendrogram="none",
			Colv=NA, Rowv=NA, labRow=colnames(Tref), ylab="",
			density.info="none", cellnote=t(notes), key.xlab=key_title, 
			key.title=key_title, 
			...)	
}
#######################################################################################################################
components.heatmap<-component.heatmap
#######################################################################################################################
component.mds<-function(
		That,
		Tref=NULL,
		D=NULL,
		sample.characteristic=NULL,
		dist.method=c("euclidean"),
		center=FALSE,
		data.ch=NULL,
		title="MDS"
){

	if(is.null(colnames(T))){
		colnames(That)<-paste("cmp", 1:ncol(That), sep="_")		
	}

	mdd<-list(That)
	
	if(!is.null(Tref)){
		if(is.null(colnames(Tref))){
			colnames(Tref)<-paste("ref", 1:ncol(Tref), sep="_")		
		}
		mdd[[2]]<-Tref
	}
	if(!is.null(D)){
		mdd[[3]]<-D
	}
	
	d<-get.distance.matrix(mdd, measure=dist.method, centered=center)

	fit <- cmdscale(d, k=2, eig=FALSE)
	
	x <- fit[,1]
	y <- fit[,2]
	
	cols<-rep("red", ncol(That))
	if(!is.null(Tref)){
		cols<-c(cols, rep("blue", ncol(Tref)))
	}
	if(!is.null(D)){
		if(is.null(sample.characteristic)){
			data.cols<-rep("grey", ncol(D))
		}else{
			if(is.numeric(sample.characteristic)){
				palette<-grey.colors(min(20, length(sample.characteristic)))
				data.cols<-palette[cut(sample.characteristic,min(20, length(sample.characteristic)))]
			}else{
				if(is.character(sample.characteristic)){
					sample.characteristic<-as.factor(sample.characteristic)
				}
				palette<-rainbow(length(levels(sample.characteristic)))
				data.cols<-palette[as.integer(sample.characteristic)]
			}
		}
		cols<-c(cols, data.cols)
	}
	
	plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
			main=title, col=cols, pch=20)
	
	text(x[1:ncol(That)], y[1:ncol(That)],
			labels=sprintf("LMC%d", 1:ncol(That)),
			pos=3-(sign(x[1:ncol(That)])))
	
	if(!is.null(Tref) && !is.null(colnames(Tref))){
		text(x[(ncol(That)+1):(ncol(That)+ncol(Tref)+1)], 
				y[(ncol(That)+1):(ncol(That)+ncol(Tref)+1)],
				labels=colnames(Tref))
	}
	
	if(!is.null(D) && !is.null(colnames(D))){
		
		start_D<-ncol(That)+ifelse(is.null(Tref), 0, ncol(Tref))+1
		end_D<-ncol(That)+ifelse(is.null(Tref), 0, ncol(Tref))+ncol(D)+1
		text(x[start_D:end_D], 
				y[start_D:end_D],
				labels=colnames(D))
	}

}
#######################################################################################################################
component.dendrogram<-function(
		That, 
		Tref=NULL, 
		dist.measure="correlation", 
		centered=TRUE, 
		label.cols=NULL
){
	if(!is.null(Tref)){
		if(is.null(colnames(Tref))){
			colnames(Tref)<-paste("Profile", 1:ncol(Tref), sep="")
		}
		d<-get.distance.matrix(list(That, Tref), measure=dist.measure, centered=centered)
	}else{
		d<-get.distance.matrix(list(That), measure=dist.measure, centered=centered)
	}
	hcl_obj<-hclust(d, method="average")
	dendr<-as.dendrogram(hcl_obj)
	
	#coloring
	if(!is.null(label.cols)){
		colLab <- function(n) {
			if (is.leaf(n)) {
				a <- attributes(n)
				labCol <- label.cols[as.integer(grepl("LMC", a$label))+1]
				attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
			}
			n
		}
		dendr<-dendrapply(dendr, colLab)
	}

	par(mar=c(7.1,4.1,4.1,2.1))
	plot(dendr, pch=NA)
	if(dist.measure=="correlation"){
		yl="1-r"
	}else{
		yl=paste(dist.measure, "distance")
	}
	title(ylab=yl)
}
#######################################################################################################################
#'
#' plotLMCs
#' 
#' A wrapper for various plotting methods for the visualization of LMCs
#' 
#' @param MeDeComSet     an object with MeDeCom results
#' @param type		     plot type, a \code{character} of length 1 (see Details)
#' @param K			     value of parameter k to use
#' @param lambda	     value of parameter lambda to use
#' @param cg_subset	     which CpG subset to use
#' @param lmc		     which LMC to use for visualization
#' @param Tref		     a matrix with reference methylomes
#' @param distance	     distance measure to use
#' @param center  		 if \code{TRUE} the LMC and reference methylome matrices will be row-centered
#' @param n.most.var	 is not \code{NA} a respective number of CpGs with the highest standard deviation will be plotted
#' @param D				 input data matrix used to derive the LMCs
#' 
#' @details
#' Available plot types include:
#' \describe{
#' \item{\bold{\code{dendrogram}}}{
#'        Dendrogram visualizing a joint hierarchical clustering of LMCs and, if available, the reference methylomes.}
#' \item{\bold{\code{heatmap}}}{
#'        Heatmap visualizing a distance between LMCs and the reference methylomes.}
#' \item{\bold{\code{mds}}}{
#'        Joint multidimensional scaling of LMCs and the reference methylomes.}
#' \item{\bold{\code{scatterplot}}}{
#'        Multi-panel scatterplot of LMCs and reference methylomes.}
#' \item{\bold{\code{extremality}}}{
#'        Barplot visualizing the value of the regularizer term for each LMC.}
#' \item{\bold{\code{distance to center}}}{
#'        Barplot visualizing a distance to the data center for each LMC. Input data matrix used to derive the LMCs should be 
#' 		  supplied as argument \code{D}.}
#' }
#' 
#' 
#' @export
plotLMCs<-function(
		MeDeComSet,
		type,
		K=NA,
		lambda=NA,
		cg_subset=1,
		lmc=NA,
		Tref=NULL,
		distance="correlation",
		center=FALSE,
		n.most.var=NA,
		D=NULL,
		sample.characteristic=NULL,
		scatter.matching=FALSE,
		scatter.smooth=TRUE,
		scatter.cg.feature=NULL,
		locus.parameters=list()
){
	plot.types<-c("dendrogram","MDS","heatmap","similarity graph", "scatterplot", "extremality", "distance to center", "locus plot")
	
	if(missing(type) || !type %in% plot.types){
		stop(sprintf("Please specify the plot type, one of \"%s\"", paste(plot.types, collapse="\", \"")))
	}
	
	That<-getLMCs(MeDeComSet,
			K,
			lambda,
			cg_subset)
	Ahat<-getProportions(MeDeComSet,
			K,
			lambda,
			cg_subset)
	
	if(is.null(colnames(That))){
		colnames(That)<-paste("LMC",1:ncol(That), sep="")
	}
	if(!is.null(Tref) && is.null(colnames(Tref))){
		colnames(Tref)<-paste("Profile", 1:ncol(Tref), sep="")
	}
	
	if(!is.na(n.most.var)){
		sds<-apply(That, 1, sd)
		top.cgs<-sds[1:n.most.var]
		That<-That[top.cgs,,drop=FALSE]
		if(!is.null(Tref)){
			Tref<-Tref[top.cgs,,drop=FALSE]
		}
		if(!is.null(D)){
			D<-D[top.cgs,,drop=FALSE]
		}
	}
	
	if(type=="dendrogram"){
		
		component.dendrogram(That, Tref, dist.measure=distance, centered=center, label.cols=c("red", "blue"))
	
	}else if(type=="heatmap"){
					
		if(distance=="correlation"){
			dist_method<-"pearson"
		}else{
			dist_method<-distance
		}
		component.heatmap(That, Tref, margins=c(5,7), method=dist_method, cexRow=1, centered=center)
					
	}else if(type=="MDS"){
		
		component.mds(That, Tref, D, dist.method=distance, center=center,
				data.ch=sample.characteristic, title="MDS")
		
	}else if(type=="scatterplot"){
		
		if(is.null(Tref)){
			stop("Tref should be supplied for this plot")
		}
		if(scatter.matching){
			perm<-MeDeCom:::corrmatch(That, Tref)
		}else{
			perm<-NA
		}
		
		ncols=min(3,ncol(Tref))
		nrows=(ncol(Tref) %/% min(3,as.integer(K)+1)) + as.integer(ncol(Tref) %% min(3,as.integer(K)+1)>0)						
		component.scatter(
				That,#[point.subs,,drop=FALSE],
				Tref,#[point.subs,,drop=FALSE], 
				lmc=lmc,
				match=perm,
				ncols=ncols,
				nrows=nrows,
				Ahat=Ahat,
				cg.feature=scatter.cg.feature,
				smooth=scatter.smooth
		)
		
	}else if(type == "matching plot"){
		
		lls<-sort(all_results$lambdas)
		Kvals<-all_results$Ks
		recovery_matrix<-matrix(0, ncol=length(Kvals),nrow=length(lls))
		match_freq_matrix<-matrix(0, ncol=length(Kvals),nrow=length(lls))
		
		blue.cols<-colorRampPalette(c("white", "blue"))
		
		for(k in Kvals){
			for(ll in sort(lls)){
				
				That<-all_results@outputs[[as.character(k)]]$T[[gr,which(all_results$lambdas==ll)]]
				
				mtch.forward<-MeDeCom:::matchLMCs(That, Tref, check=FALSE)
				mtch.reverse<-MeDeCom:::matchLMCs(Tref, That, check=FALSE)
				num.mutual<-sum(1:ncol(That)==mtch.reverse[mtch.forward])
				match.freq<-num.mutual/(min(ncol(That), ncol(Tref)))
				recovery_matrix[which(lls==ll),which(Kvals==k)]<-num.mutual
				match_freq_matrix[which(lls==ll),which(Kvals==k)]<-match.freq
				#as.numeric(length(mtch)==length(unique(mtch)) && length(mtch)==ncol(Tref))
			}}
		N_COL_BINS<-10

		image(x=1:length(Kvals),y=1:length(lls),
				t(match_freq_matrix),
				#col=gray.colors(N_COL_BINS-1, start=1, end=0, gamma=1.5),	
				col=blue.cols(N_COL_BINS-1),
				breaks=(1:N_COL_BINS)*(1/N_COL_BINS), 
				axes=FALSE,xlab="k", ylab="lambda")
		for(ri in 1:nrow(recovery_matrix)){
			for(ci in 1:ncol(recovery_matrix)){
				text(x=ci,y=ri,t(recovery_matrix)[ci,ri])
			}
		}
		axis(1, at = 1:ncol(recovery_matrix), labels=as.character(Kvals), tick=FALSE)
		axis(2, at = 1:nrow(recovery_matrix), labels=sprintf("%g", lls),tick=FALSE, las=2)
		
	}else if(type == "similarity graph"){
		d<-get.distance.matrix(lit(That, Tref), measure=distance, centered=centered)
		if(distance=="euclidean"){
			d <- d/max(d)				
		}else if(distance=="correlation"){
			diag(d)<-0
		}
		
		require(igraph)
		G <- graph.adjacency(d, mode="upper", weighted=TRUE,)
		G <-delete.edges(G, which(E(G)$weight <= as.numeric(input$minGraphSimilarity)))
		#G <- remove.edges(G, V(G)[ degree(G)==0 ])
		lo<-layout.fruchterman.reingold(G, weights=rep(1, ecount(G)))
		#lo<-layout.spring(G, repulse=TRUE)
		plot(G, layout=lo,
				vertex.size=3, edge.arrow.size=0.2,
				vertex.label=colnames(mdd), rescale=FALSE,
				xlim=range(lo[,1]), ylim=range(lo[,2]), vertex.label.dist=0.2,
				vertex.label.color="black")
		
		#E(G)$weight <- E(G)$weight - 1
		#G <- graph.adjacency(adjmat, mode="undirected", weighted=TRUE)
		#G <- remove.vertices(G, V(G)[ degree(G)==0 ])
	
	}else if(type=="distance to center"){
		
		center<-rowMeans(D)
		distances<-colSums(sweep(That, 1, center, "-")^2)
		barplot(distances, las=2, names.arg=as.character(1:ncol(That)), xlab="r", ylab="||t_r - mean(D)||")
		
	}else if(type=="extremality"){
		
		extremality<-colMeans(That*(1-That))
		barplot(extremality, las=2, names.arg=as.character(1:ncol(That)), xlab="r", ylab="t_r(1-t_r)")
		
	}else if(type=="locus plot"){
		
		locus.parameters[["That"]]<-That
		locus.parameters[["cgs"]]<-MeDeComSet@parameters$cg_subset_lists[match(cg_subset, MeDeComSet@parameters$cg_subset_lists)]
		locus.parameters[["D"]]<-D
		
		do.call("locus_plot", locus.parameters)
	}
}
#######################################################################################################################
##
## locus_plot
##
## Plot results for a selected locus
##
## param That a matrix of LMCs
## param ann CpG annotation
## param cgs indices of CpGs data for which are present in That with respect to ann
## param locus.chr chromosome
## param locus.start start coordinate
## param locus.end end coordinate
## param flank.start number of basepairs to extend the locus upstream
## param flank.end number of basepairs to extend the locus downstream
## param comp.cols color code for LMCs
## param legend.pos location of the legend, in accordance with \link{legend}
## param Tstar matrix of reference profiles
## param D matrix of input methylation data used to produce That
## param plot.genes if \cs{TRUE} a track with gene locations will be plotted
## param ann.genes gene annotation necessary for the gene plotting
##
locus_plot<-function(
		That,
		ann,
		cgs,
		locus.chr
		,locus.start
		,locus.end
		,locus.name=sprintf("%s:%d-%d", locus.chr, locus.start,locus.end)
		,locus.forward=FALSE
		,flank.start=1000
		,flank.end=1000
		,comp.cols=rainbow(ncol(That))
		,legend.pos="topleft"
		,Tstar=NULL
		,D=NULL
		,plot.genes=FALSE
		,ann.genes=NULL
){
	
	#require(RnBeads)
	#ann<-rnb.annotation2data.frame(rnb.get.annotation("probes450"))
	
	ann.cgs<-ann[cgs,]
	cgs.locus<-which(ann.cgs$Chromosome==locus.chr & ann.cgs$Start>locus.start-flank.start & ann.cgs$End<locus.end+flank.end)
	ann.cgs.locus<-ann.cgs[cgs.locus, ]
	
	if(plot.genes && !is.null(ann.genes)){
		#ann.genes<-rnb.annotation2data.frame(rnb.get.annotation("genes"))
		
		genes.locus.start<-which(ann.genes$Chromosome==locus.chr & ann.genes$Start>locus.start-flank.start & ann.genes$Start<locus.end+flank.end)
		genes.locus.end<-which(ann.genes$Chromosome==locus.chr & ann.genes$End>locus.start-flank.start & ann.genes$End<locus.end+flank.end)
		genes.locus.cover<-which(ann.genes$Chromosome==locus.chr & ann.genes$Start<locus.start-flank.start & ann.genes$End>locus.end+flank.end)
		ann.genes.locus<-ann.genes[unique(c(genes.locus.start, genes.locus.end, genes.locus.cover)), ]

	}
	
	plot(NA,NA, 
			type="n",lty=3,pch=15, col="white", 
			ylim=c(-0.3-plot.genes*0.2,1), 
			xlim=c(locus.start-flank.start, locus.end+flank.end), 
			ylab="methylation level", xlab=locus.chr,
			yaxt="n")
	axis(2, at=c(0,0.25,0.5,0.8,1))
	abline(h=0)
	for(cmp in 1:ncol(That)){
		lines(ann.cgs.locus$Start, That[cgs.locus,cmp], type="o",lty=2, pch=15, cex=1, col=comp.cols[cmp])
	}
	if(!is.null(Tstar)){
		prof.cols<-rainbow(ncol(That)+ncol(Tstar))[(ncol(That)+1):(ncol(That)+ncol(Tstar))]
		for(prof in 1:ncol(Tstar)){
			lines(ann.cgs.locus$Start, Tstar[cgs.locus,prof], type="o", lty=4, pch=19, cex=1, col=prof.cols[prof])
		}
	}else{
		prof.cols<-character()
	}
	if(!is.null(D)){
		for(prof in 1:ncol(D)){
			lines(ann.cgs.locus$Start, D[cgs.locus,prof], type="o", lty=3, pch=21, cex=0.5, col="lightgrey")
		}
	}
	arrows(if(locus.forward) locus.start else locus.end, 
			-0.15, 
			if(locus.forward) locus.end else locus.start, 
			-0.15, length = 0.1)
	text(locus.start+(locus.end-locus.start)/2, -0.2, locus.name)
	
	if(plot.genes && nrow(ann.genes.locus)>0){
		for(row in 1:nrow(ann.genes.locus)){
			print(nrow(ann.genes.locus))
			gene.forward<-ann.genes.locus[row,"Strand"]=="+"
			line_y_pos<-(-1)*(0.25)-(((row-1) %% min(4, nrow(ann.genes.locus))))*(0.25/min(4,nrow(ann.genes.locus)))
			print(line_y_pos)
			arrows(if(gene.forward) ann.genes.locus[row,"Start"] else ann.genes.locus[row,"End"], 
					line_y_pos, 
					if(gene.forward) ann.genes.locus[row,"End"] else ann.genes.locus[row,"Start"], 
					line_y_pos, length = 0.1)
			lab_pos<-ann.genes.locus[row,"Start"]+(ann.genes.locus[row,"End"]-ann.genes.locus[row,"Start"])/2
			if(lab_pos<locus.start-flank.start || lab_pos>locus.end+flank.end) 
				if(ann.genes.locus[row,"Start"]>locus.start-flank.start )				
					lab_pos<-(locus.end+flank.end+ann.genes.locus[row,"Start"])/2
				else if(ann.genes.locus[row,"End"]<locus.end)
					lab_pos<-(locus.start-flank.start+ann.genes.locus[row,"End"])/2
				else lab_pos<-((locus.start-flank.start)+(locus.end+flank.end))/2
			if(!is.na(ann.genes.locus[row,"symbol"]))
				text(lab_pos, line_y_pos-0.05, ann.genes.locus[row,"symbol"])
			else
				text(lab_pos, line_y_pos-0.05, ann.genes.locus[row,"ID"])
		}
	}
	legend(legend.pos, pch=c(rep(15,ncol(That)),if(!is.null(Tstar)) rep(20,ncol(Tstar)) ), 
			col=c(comp.cols,prof.cols), lty=c(rep(2, ncol(That)), if(!is.null(Tstar)) rep(4, ncol(Tstar))),
			legend=if(!is.null(colnames(That))) c(colnames(That), colnames(Tstar)) else c(sprintf("LMC %d", 1:ncol(That)), colnames(Tstar)))
	return(invisible(list(locus.data=That[cgs.locus,], locus.ann=ann.cgs.locus)))
}
#######################################################################################################################
proportion.lineplot<-function(
		MeDeComSet,
		K,
		lambda,
		lmc,
		cg_subset=1,
		Aref=NULL,
		Tref=NULL,
		ref.profile=NA, 
		reorder="increasing", 
		legend.pos="topleft", 
		add=T, 
		assignment.method="pearson", 
		...){
	
	ltys<-legs<-cols<-pchs<-c()
	
	applyAllList<-list()
	
	if(!is.null(Aref)){
		if(is.null(rownames(Aref))){
			rownames(Aref)<-sprintf("Ref. profile %s", 1:nrow(Aref))
		}
		yl<-paste(rownames(Aref)[ref.profile], collapse=" + ")
		applyAllList$truth<-list(A=Aref)
	}
	
	if(!is.null(MeDeComSet@outputs$Aprime)){
		#all$regression<-MeDeCom:::factorize.regr(meth.data[ind,sample_subset],trueT[ind,])
		applyAllList$regression<-list(A=MeDeComSet@outputs$Aprime)
		if(is.null(Aref)){
			yl<-paste(rownames(applyAllList$regression$A)[ref.profiles], collapse=" + ")
		}
	}
	
	if(!is.na(lmc)){
		applyAllList$MeDeCom<-list()
		applyAllList$MeDeCom$A<-getProportions(MeDeComSet, K, lambda, cg_subset)
		if(is.null(Aref) && is.null(applyAllList$regression)){
			yl<-paste(rownames(applyAllList$MeDeCom$A)[lmc], collapse=" + ")
		}
	}
	
	if(!is.na(reorder)){
		if(!is.null(Aref)){
			sample_order<-order(applyAllList$truth$A[ref.profile,], decreasing=reorder=="decreasing")
		}else if(!is.null(applyAllList$regression)){
			sample_order<-order(applyAllList$regression$A[ref.profile,], decreasing=reorder=="decreasing")
		}else{
			sample_order<-order(applyAllList$MeDeCom$A[lmc,], decreasing=reorder=="decreasing")
		}
	}else{
		sample_order<-1:ncol(applyAllList$MeDeCom$A)
	}
	nsamp<-max(sapply(applyAllList, function(res) if(!is.null(res)) ncol(res$A) else 0))
	
	xl="samples"
	
	plot(NA, ylim=c(0,1), xlim=c(1,nsamp), xlab=xl, 
			ylab=sprintf("Proportion %s", yl),
			yaxt="n",
			...)
	axis(2, las=2)			
	for(mm in names(applyAllList)){
			if(!is.null(applyAllList[[mm]])){
				if(mm %in% c("truth","regression")) {
					peer<-ref.profile
				}else if(is.na(lmc) && !is.null(Tref)){
					#peer<-which.max(apply(applyAllList[[mm]]$T, 2, cor, y=trueT[,component], use="pairwise.complete.obs"))
					if(mm %in% c("IntFac", "VertexSearch")){
						assignment.method<-"anova"
					}else{
						assignment.method<-"spearman"
					}
					peer<-which(corrmatch(applyAllList[[mm]]$T, Tref, assignment.method)==ref.profile)
					if(length(peer)>1){
						peer<-peer[which.max(corrmatch(applyAllList[[mm]]$T[,peer], Tref, assignment.method, return.vals=TRUE))]
					}
				}else{
					peer<-lmc
				}
				legs<-c(legs,mm)
				alg.col<-ALGORITHM.COLS[mm]
				alg.pch<-ALGORITHM.PCH[mm]
				pchs<-c(pchs, alg.pch)
				cols<-c(cols, alg.col)
				if(mm=="truth"){
					ltys<-c(ltys, 1)
				}else{
					ltys<-c(ltys, 2)
				}
				if(length(peer)==1){
					lines(applyAllList[[mm]]$A[peer,sample_order], 
							col=alg.col, type="o", pch=alg.pch, lty=2)
				}else{
					lines(colSums(applyAllList[[mm]]$A[peer,sample_order]), 
							col=alg.col, type="o", pch=alg.pch, lty=2)
				}
			}
	}
	if(!is.na(legend.pos)){
		legend(legend.pos, legend=legs, col=cols, pch=pchs, lty=ltys)
	}
}
#######################################################################################################################
proportion.scatter<-function(Ahat, Aref, lmc, ref.profile){
	
	ref_data<-colSums(Aref[ref.profile,,drop=FALSE])
	est_data<-colSums(Ahat[lmc,,drop=FALSE])
	
	if(is.null(rownames(Aref))){
		rownames(Aref)<-sprintf("Ref. profile %s", 1:nrow(Aref))
	}
	
	plot(ref_data, est_data, xlim=c(0,1), ylim=c(0,1), 
			xlab=paste(rownames(Aref)[ref.profile], collapse=" + "), ylab=paste(sprintf("LMC%d", lmc), collapse=" + "))
	abline(h=(0:10)/10, lwd=0.5, col="grey")
	abline(v=(0:10)/10, lwd=0.5, col="grey")
	
	fitData<-list()
	fitData$proportion1<-ref_data
	fitData$proportion2<-est_data
	fit<-lm(proportion2~proportion1, fitData)
	#abline(fit)
	rp <- vector('expression',1)
	rp[1] <- substitute(expression(italic(r) == MYVALUE), 
			#list(MYVALUE = format(summary(fit)$r.squared,dig=3)))[2]
			list(MYVALUE = format(cor(fitData$proportion1,fitData$proportion2),dig=3)))[2]
	legend("topleft", legend=rp, bty="n")
	
}
#######################################################################################################################
proportion.heatmap<-function(
		Ahat, 
		sample.characteristic=NULL,
		clusterCols=FALSE,
		clusterRows=FALSE
		){
	
	if(!is.null(sample.characteristic)){
		data.ch<-sample.characteristic
		
		if(is.numeric(data.ch)){
			palette<-grey.colors(start=0.9, end=0.1, min(20, length(data.ch)))
			data.cols<-palette[cut(data.ch,min(20, length(data.ch)),labels=FALSE)]
			data.cols[is.na(data.cols)]<-"#ffffff"
		}else{
			if(is.character(data.ch)){
				data.ch<-factor(data.ch, levels=unique(data.ch))
			}
			palette<-rainbow(length(levels(data.ch)))
			data.cols<-palette[as.integer(data.ch)]
		}
	}else{
		data.cols<-NA
	}
	
	N_COL_BINS<-50
	hm_args<-list(x=Ahat,
			scale="none", trace="none", density.info="none", 
			keysize=1, key.xlab="proportion",  key.title="proportion",
			Colv=clusterCols, Rowv=clusterRows, 
			dendrogram=c("none", "column", "row", "both")[1+clusterCols+2*clusterRows],
			cexRow=2,
			#margins=if(!all(is.na(data.cols))) c(7,12) else c(5,5),
			margins=c(7,12),
			#col=c("white",gray.colors(N_COL_BINS-1, start=0.9, end=0.3), "black"),
			col=gray.colors(N_COL_BINS-1, start=1, end=0, gamma=1.5),
			#breaks=c(0.0001,(1:N_COL_BINS)*(1/N_COL_BINS),0.9999),
			breaks=(1:N_COL_BINS)*(1/N_COL_BINS))
	if(!all(is.na(data.cols))) {
		hm_args$ColSideColors=data.cols 
	}
	do.call(gplots::heatmap.2, hm_args)
	
	if(!all(is.na(data.cols))){
		
		if(is.factor(data.ch)){
			var_levels<-levels(data.ch)
		}else if(is.numeric(data.ch)){
			var_levels<-levels(cut(data.ch,min(20, length(data.ch))))
		}
		nc<-1+length(var_levels)%/%10
		xcoord<-0.85-0.05*nc
		ycoord<-1.15-0.05*nc #/(length(data.ch)%%10+1)
		legend(x=xcoord, y=ycoord, var_levels, col=palette, pch=15,
				ncol=nc, 
				xpd=TRUE, horiz=FALSE, bty="n", title=attr(data.ch, "name"))
	}
}
#######################################################################################################################
proportion.feature.corr<-function(Ahat, lmc, data.ch, includeRegressionLine=FALSE){
	
	if(is.null(rownames(Ahat))){
		rownames(Ahat)<-sprintf("LMC%d", 1:nrow(Ahat))
	}
	yl<-paste(rownames(Ahat)[lmc], collapse=" + ")
	if(is.numeric(data.ch)){
		plot(data.ch, Ahat[lmc,], 
				xlab="", 
				ylab=yl, las=2)
		if(includeRegressionLine){
			fitData<-list()
			fitData$feature<-data.ch
			fitData$proportion<-Ahat[lmc,]
			fit<-lm(proportion~feature, fitData)
			abline(fit)
			rp <- vector('expression',1)
			rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
					list(MYVALUE = format(summary(fit)$r.squared,dig=3)))[2]
			legend("topright",legend=rp, bty="n")
		}
	}else if(is.factor(data.ch) || is.character(data.ch)){
		if(is.character(data.ch)){
			data.ch<-as.factor(data.ch)
		}
		boxplot(prop~pheno, data.frame(prop=Ahat[lmc,], pheno=data.ch),
				xlab=attr(data.ch, "name"), ylab=yl, las=2)			
	}
}			
#######################################################################################################################
#'
#' plotProportions
#' 
#' A wrapper for various plotting methods for the visualization of mixing proportions
#' 
#' @param MeDeComSet     an object with MeDeCom results
#' @param type		     plot type, a \code{character} of length 1 (see Details)
#' @param K			     value of parameter k to use
#' @param lambda	     value of parameter lambda to use
#' @param cg_subset	     which CpG subset to use
#' @param lmc		     which LMC to use for visualization
#' @param Aref		     a matrix with reference methylomes
#'  
#' @details
#' Available plot types include:
#' \describe{
#' \item{\bold{\code{heatmap}}}{
#'        Lineplot of proportions recovered by MeDeCom and reference proportions.}
#' \item{\bold{\code{barplot}}}{
#'        Stacked barplot of proportions recovered by MeDeCom.}
#' \item{\bold{\code{lineplot}}}{
#'        Lineplot of proportions recovered by MeDeCom and, if available, reference proportions.}
#' \item{\bold{\code{scatterplot}}}{
#'        Lineplot of proportions recovered by MeDeCom and reference proportions.}
#' }
#' 
#' 
#' @export
plotProportions<-function(
		MeDeComSet,
		type,
		K,
		lambda,
		cg_subset=1,
		lmc=NA,
		Aref=NULL,
		ref.profile=NA, 
		assignment.method="pearson",
		sample.characteristic=NULL,
		heatmap.clusterCols=FALSE,
		heatmap.clusterRows=FALSE
){
	
	Ahat<-getProportions(MeDeComSet, K, lambda, cg_subset)
	
	if(type=="heatmap"){
		
		proportion.heatmap(Ahat, 
				sample.characteristic,
				clusterCols=heatmap.clusterCols,
				clusterRows=heatmap.clusterRows)
		
	}else if(type=="barplot"){
		
		barplot(Ahat, las=2, xlab="Samples", ylab="Proportion")
		
	}else if(type=="lineplot"){
		
		proportion.lineplot(
				MeDeComSet,
				K,
				lambda,
				cg_subset,
				lmc=lmc,
				Aref=Aref,
				ref.profile=ref.profile, 
				reorder=TRUE, 
				assignment.method=assignment.method)
		
	}else if(type=="scatterplot"){
		
		if(is.null(Aref) && is.null(Tref)){
			stop("one of Aref or Tref should be supplied for this plot")
		}else if(is.na(lmc) || is.na(ref.profile)){
			stop("one of lmc or ref.profile should be supplied for this plot")
		}
		
		proportion.scatter(Ahat, Aref, lmc, ref.profile)

	}else if(type=="sample characteristics"){
		
		if(is.null(sample.characteristic)){
			stop("sample.characteristic should be supplied for this plot")
		}
		
		proportion.feature.corr(Ahat, lmc, sample.characteristic)
	}
}
#######################################################################################################################