#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#define EPSILON DBL_EPSILON

SEXP ZeroIn(SEXP x, SEXP y, SEXP flow, SEXP fup, SEXP fcall, SEXP st,
               SEXP stol, SEXP smaxit, SEXP env)
{
   SEXP sexp_lvec, sexp_uvec, sexp_mvec, ans;
   int maxit = INTEGER(smaxit)[0];
   int T = INTEGER(st)[0];
   double tol = REAL(stol)[0];
   int done[T];
   double fa[T], fc[T], a[T], c[T];
   int count = 0;
   double prev_step, new_step, tol_act;
   double p, q, b, fb;
   
   // Initialize fa, fb and a,b

   for(int k=0; k < T; k++) {
      fa[k] = REAL(flow)[k];
      //fb[k] = REAL(fup)[k];
      fc[k] = fa[k];
      a[k] = REAL(x)[k];
      //b[k] = REAL(y)[k];
      c[k] = a[k];
   }
   
   PROTECT(ans = allocVector(REALSXP, T + 1));

   for(int j=0; j < maxit; j++) {
       defineVar(install("x"), y, env);
       PROTECT(sexp_mvec = eval(fcall, env));
       for(int i=0; i < T; i++) {
           // on first iteration, set done[i] = 0, for each i;
           if(j==0) {
              done[i] = 0;
           }
           if(done[i] == 0) {
               fb = REAL(sexp_mvec)[i];
               b = REAL(y)[i];

               if( (fb > 0 && fc[i] > 0) || (fb < 0 && fc[i] < 0) ) {
	                /* Adjust c for it to have a sign opposite to that of b */
	                c[i] = a[i];  fc[i] = fa[i];
	           }
               prev_step = b - a[i];
               
               if( fabs(fc[i]) < fabs(fb) )  {				
	               a[i] = b;  b = c[i];  c[i] = a[i];	/* best approximation		*/
	               fa[i] = fb;  fb = fc[i];  fc[i] = fa[i];
	           }
	   
	           tol_act = 2*EPSILON*fabs(b) + tol/2;
	           new_step = (c[i] - b)/2;
	           
     	       if( fabs(new_step) <= tol_act || fb == (double)0 ) {
	                 done[i] = 1;
	                 count++;  
	           }
	           else {
	                if( fabs(prev_step) >= tol_act  && fabs(fa[i]) > fabs(fb) ) {	
	                    register double t1,cb,t2;
	                    cb = c[i] - b;
	                    if( a[i]==c[i] ) {		
		                    t1 = fb/fa[i];		
		                    p = cb*t1;
		                    q = 1.0 - t1;
	                    }
	                    else {			/* Quadric inverse interpolation*/
                            q = fa[i]/fc[i];  t1 = fb/fc[i];	 t2 = fb/fa[i];
		                    p = t2 * ( cb*q*(q-t1) - (b-a[i])*(t1-1.0) );
		                    q = (q-1.0) * (t1-1.0) * (t2-1.0);
	                    }
	                    if( p>(double)0 ) {
		                     q = -q;			/* opposite sign; make p positive */
	                    }
	                    else  {		
		                     p = -p;			/* q				*/
                        }
	                    if( p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2) )	{
		                     new_step = p/q;
		                }			
	                }
	                if( fabs(new_step) < tol_act) {
	                    if( new_step > (double)0 )	/* than tolerance		*/
		                    new_step = tol_act;
	                    else
		                    new_step = -tol_act;
	                }
	                a[i] = b;	fa[i] = fb;			/* Save the previous approx. */
	                b += new_step;	

	                REAL(y)[i] = b;
	                REAL(ans)[i] = b;           
	           }
           }
       }
       UNPROTECT(1);
       if(count==T) {
          // all roots have been found if count=T
          // Rprintf("[%d]", j+1); 
         // break;
          REAL(ans)[T] = j + 1;
          UNPROTECT(1);
          return(ans);
       }
   }
   REAL(ans)[T] = -1;
   UNPROTECT(1);
   return(ans);
}

