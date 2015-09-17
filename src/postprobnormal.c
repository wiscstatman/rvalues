#include <R.h>
#include <Rmath.h>

void postprobnormal(double *XX, double *std_err, double *support, double *mix_prop,
              int *nn, int *TT, double *num, double *post, double *loglik) {  

  int n=*nn, m=*TT, i, j, min_ind=0, count=0;
  double x, ss, ff, denom_sum, min_val=0.0;

  *loglik = -(double)n * 0.91893853320467274178; /* n/2 times log(2pi) */
 
  for (i=0; i<n; i++){
    x=XX[i];
    ss = std_err[i];
    for (j=0; j<m; j++) {
      ff = (x-support[j])/ss;
      ff = (ff*ff)/2.0;
      num[j] = -ff;

      if (count==0 || ff < min_val) {
         if(mix_prop[j] > 0) {        
             min_val = ff;
             min_ind = j;
             count += 1;
         }
      }
    }
    denom_sum = 1.0;
    for (j=0; j<m; j++) {
      if (j==min_ind) 
        num[j] = 1.0;
      else {
        ff = num[j] + min_val - log(mix_prop[min_ind]);
        num[j] = mix_prop[j]*exp(ff);
        denom_sum += num[j];
      }
    }
    for (j=0; j<m; j++) {
      post[i + n*j] = num[j]/denom_sum;
    }
    
    *loglik += log(denom_sum) - min_val + log(mix_prop[min_ind]) - log(ss);
  }
}

