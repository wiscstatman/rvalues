#include <R.h>
#include <Rmath.h>

void postprobbinomial(double *XX, double *ntrials, double *support, double *mix_prop,
              int *nn, int *TT, double *num, double *post_prob, double *loglik) {  

  int n=*nn, m=*TT, i, j, max_ind=0;
  double x, mm, ff, denom_sum, max_val=0.0;

  *loglik = 0.0;
 
  for (i=0; i<n; i++){
    x=XX[i];
    mm = ntrials[i];
    for (j=0; j<m; j++) {
      ff = x*log(support[j]) + (mm - x)*log(1 - support[j]);
      num[j] = ff + log(mix_prop[j]);
      
      if (j==0 || ff > max_val) {        
          max_val = ff;
          max_ind = j;
      }
    }
    denom_sum = 1.0;
    for (j=0; j<m; j++) {
      if (j==max_ind) 
          num[j] = 1.0;
      else {
          ff = num[j] - max_val - log(mix_prop[max_ind]);
          num[j] = exp(ff);
          denom_sum += num[j];
      }
    }
    for (j=0; j<m; j++) {
       post_prob[i + n*j] = num[j]/denom_sum;
    }
    *loglik += log(denom_sum) + lgamma(mm+1) - lgamma(mm - x + 1) - lgamma(x + 1) + max_val + log(mix_prop[max_ind]); 
  }
}
