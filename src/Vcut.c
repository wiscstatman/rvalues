#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP Vcut(SEXP Vmat, SEXP lamfun, SEXP nunits, SEXP ngrid, SEXP Agrid) {
  /* Vmat - nunits x ngrid matrix */
  /* lam_fun - vector of length ngrid */
  SEXP ans;
  int nrow, ncol, ii, jj;
  double tst, g0, gdelt, slop;
  
  nrow = INTEGER(nunits)[0];
  ncol = INTEGER(ngrid)[0];
  
  PROTECT(ans = allocVector(REALSXP, nrow));
  
  for (ii=0; ii < nrow; ii++) {
     for(jj=0; jj < ncol; jj++)  {
         tst = REAL(Vmat)[ii + nrow*jj] - REAL(lamfun)[jj];
         // maybe change the tst>= 0 condition or
         // compare(Vmat) with (lamfun) directly
         // Interpolation:
         if(tst >= 0) {
            //tstp1 = REAL(Vmat)[ii + nrow*jj - 1] - REAL(lamfun)[jj - 1];
            // check if "jj=0"
            if(jj > 0) {
                g0 = REAL(Vmat)[ii + nrow*jj - 1] - REAL(lamfun)[jj - 1]; 
                gdelt = REAL(Vmat)[ii + nrow*jj] - REAL(lamfun)[jj] - g0;
                slop = (REAL(Agrid)[jj] - REAL(Agrid)[jj - 1])/gdelt;
                REAL(ans)[ii] = REAL(Agrid)[jj - 1] - g0*slop;
            }
            else {
                REAL(ans)[ii] = REAL(Agrid)[0];
            }
            //INTEGER(ans)[ii] = jj + 1;
            break;
         }
         if(jj == ncol - 1) {
            //INTEGER(ans)[ii] = jj + 1;
            REAL(ans)[ii] = REAL(Agrid)[jj + 1];
            // If no root exists across columns return max 
         }
     }
  }
  
  UNPROTECT(1); 
  return(ans);
} 
