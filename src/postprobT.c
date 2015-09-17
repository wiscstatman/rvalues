#include <R.h>
#include <Rmath.h>

void postprobT(double *XX, double *std_err, double *df, double *support, double *mix_prop,
              int *nn, int *TT, double *num, double *post, double *loglik) {  

  int n=*nn, m=*TT, i, j, min_ind=0;
  double x, ss, pw, deg, ff, denom_sum, min_val=0.0;

  *loglik = -(double)n * lgamma(0.5); // lgamma(1/2) = log(sqrt(pi))
 
  for (i=0; i<n; i++){
    x=XX[i];
    ss = std_err[i];
    deg = df[i];
    pw = (deg + 1)/2.0;
    for (j=0; j<m; j++) {
      ff = (x-support[j])/ss;
      ff = (ff*ff)/deg;
      ff = pw*log(1 + ff);
      num[j] = log(mix_prop[j]) - ff;

      if (j==0 || ff < min_val) {        
         min_val = ff;
         min_ind = j;
      }
    }
    denom_sum = 1.0;
    for (j=0; j<m; j++) {
      if (j==min_ind) 
        num[j] = 1.0;
      else {
        ff = num[j] + min_val - log(mix_prop[min_ind]);
        num[j] = exp(ff);
        denom_sum += num[j];
      }
    }
    for (j=0; j<m; j++) {
      post[i + n*j] = num[j]/denom_sum;
    }
    
    *loglik += log(denom_sum) - min_val + log(mix_prop[min_ind]) - log(ss) + lgamma(pw) - lgamma(pw - 0.5) - log(deg)/2.0;
  }
}
