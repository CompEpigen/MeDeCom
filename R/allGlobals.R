
#######################################################################################################################
# GLOBALS
#######################################################################################################################

ALGORITHMS<-c(
		"truth",
		"regression",
		"houseman2012",
		"houseman2016",
		"MeDeCom",
		"MeDeCom.quadPen",
		"MeDeCom.cppTAfact", 
		"HLasso", 
		"IntFac", 
		"IntEmpirical",
		"Resample",
		"VertexSearch")

T_METHODS<-c(
		NA,
		NA,
		NA,
		NA,
		"quadPen",
		"quadPen",
		"cppTAfact",
		"Hlasso", 
		"integer", 
		"empirical",
		"resample",
		NA)

names(T_METHODS)=ALGORITHMS

ALGORITHM.COLS=c(
		"truth"="black",
		"regression"="blue",
		"houseman2012"="deepskyblue",
		"houseman2016"="skyblue",
		"MeDeCom"="red",
		"MeDeCom.quadPen"="red",
		"MeDeCom.cppTAfact"="tomato", 
		"HLasso"="orange", 
		"IntFac"="green", 
		"IntEmpirical"="plum",
		"Resample"="magenta",
		"VertexSearch"="brown")

ALGORITHM.PCH=c(15,0,0,0,2,2,2,3,4,5,6,1)

names(ALGORITHM.PCH)<-ALGORITHMS

# a basis for the lambda parameter grid
RELGRID = c(
		1E-10, 1E-5, 1E-4, 5E-4, 1E-3, 2E-3, 5E-3, 
		1E-2, 2E-2, 5E-2, 1E-1, 2E-1, 5E-1, 1)

PERFORMANCE_MEASURES<-c("Objective"="Fval", "RMSE"="rmse", "CV error"="cve", "RMSE, T"="rmseT", "MDC, T"="dist2C", "MAE, A"="maeA")


#######################################################################################################################
