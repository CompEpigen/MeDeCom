/******************************************************
 * f = mexProjSplxBox(x, l, u)
 *
 * Input: x - a vector in Rn.
 * Output: f - Projection of x onto {y: sum(y) == 1, l <= y <= u}
 *
 *********************************************************/

#include <stdio.h>
#include <cstddef>
#include "dynblas.h"
#include <omp.h>
#include <limits>
#include <Rcpp.h>
#include <float.h>
using namespace Rcpp;

inline double dabs(double x){
    if ( x < 0 )
        x = -x;
    return x;
}

void Proj(double* f, double* x, double* l, double* u, double* y, double* z, double* p, double* q, int k) {
    /* 
     * projection uses Dykstra's algorithm
     * 
     */
   
     double tol1 = 0;
     double tol2 = 0;	
     double tolcur = 0;	
     double tolterm = 1E-12; 
     double mean = 0;
     double sum = 0;	 
     double delta = 0;
     double overk = 1/(double)k;
     /* int maxiters = 10000; */
     /* int itercur = 0; */	

	for(int i = 0; i < k; i++){
	 	p[i] = 0;
       	}
    

	for(int i = 0; i < k; i++){
	 	q[i] = 0;
       	}

      
     for(int i = 0; i < k; i ++){
       f[i] = x[i];		
     }

		

     do{/* break; */
       /* itercur = itercur + 1; */
       /* reset before each new round */
       mean = 0;
       tol1 = 0;
       tol2 = 0;

       /* copy old iterate */
       for(int i = 0; i < k; i++){
	 y[i] = f[i];
       }
	

       /* projection 1: box constraints, primal */

       for(int i = 0; i < k; i++) {
	 
	 z[i] = f[i] + p[i];

	 if(z[i] < l[i]){
	   z[i] = l[i];
	 }

	 if(z[i] > u[i]){
	   z[i] = u[i];
	 }	
       }

       /* projection 1: box constraints, dual */ 

       for(int i = 0; i < k; i++) {
	 p[i] = f[i] + p[i] - z[i];
       }

       /* projection 2: sum constraint, primal */

       for(int i = 0; i < k; i++) {
	 f[i] = z[i] + q[i]; 
	 mean = mean + f[i] * overk; 
       }



       for(int i = 0; i < k; i++) {
	 f[i] = f[i] - mean + overk;
       }

       /* projection 2: sum constraint, primal */

       for(int i = 0; i < k; i++) {
	 q[i] = z[i] + q[i] - f[i];
       }

       /* tolerance checking */

       for(int i = 0; i < k; i++) {
	 delta = dabs(f[i] - y[i]);
	 if(delta > tol1){
	   tol1 = delta;
	 }
       }

       for(int i = 0; i < k; i++) {
	 delta = l[i] - f[i];
	 if(delta > tol2){
	   tol2 = delta;
	 }

	 delta  = f[i] - u[i];
	 if(delta > tol2){
	   tol2 = delta;
	 }
       }

       sum = 0;	
       for(int i=0; i < k; i++){		
	sum = sum + f[i];
       }

       if(dabs(sum - 1) > tol2){
	tol2 = dabs(sum - 1);	 	
       }	

       if(tol1 < tol2){
	 tolcur = tol2;
       }
       else{
	   tolcur = tol1;	
       }

     /* printf ("tols: %1.6f \n", tol2); */

     } while(tolcur > tolterm);/* && itercur < maxiters); */     
}


//[[Rcpp::export]]
 NumericMatrix RProjSplxBox(NumericMatrix Xinp, NumericVector linp, NumericVector uinp) {
    
	Rcpp::NumericMatrix Xi(clone(Xinp));
	Rcpp::NumericVector li(clone(linp));
	Rcpp::NumericVector ui(clone(uinp));

    /*Get pointer to input data*/

	double* x =  Xi.begin();
	double* l =  li.begin();
	double* u =  ui.begin();

    //int cols = Xi.ncol();
    int rows = Xi.nrow();
	
    double* f = NULL;      /* output  */
    
    /* create the output vector */
    Rcpp::NumericMatrix newX((int)rows,1);
    f = newX.begin();

    /* int* ix = (int*) malloc(rows * sizeof(int)); */
    double* y = (double*) malloc(rows * sizeof(double));
    double* z = (double*) malloc(rows * sizeof(double)); 
    double* p = (double*) malloc(rows * sizeof(double));
    double* q = (double*) malloc(rows * sizeof(double));
      
    Proj(f, x, l, u, y, z, p, q, rows); 
    
    /* free(ix); */
    free(z);
    free(y);
    free(p);
    free(q);
    
    return(newX);
    
}
