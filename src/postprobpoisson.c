#include <R.h>
#include <Rmath.h>

// Function to compute a matrix of posterior probabilities in the Poisson
// mixture model. The (i,j) entry of the output matrix will denote
// the posterior probability that the class label is i given that 
// the support parameter is j.

void postprobpoisson(double *XX, double *eta, double *support, double *mix_prop,
              int *nn, int *TT, double *num, double *post_prob, double *loglik) {  

  int n=*nn, m=*TT, i, j, max_ind=0;
  double x, et, ff, denom_sum, max=0.0;

  *loglik = 0.0;
 
  for (i=0; i<n; i++){
    x=XX[i];
    et = eta[i];
    for (j=0; j<m; j++) {
      ff = x*log(support[j]) - et*support[j];
     
      if (j==0 || ff > max) {        
          max = ff;
          max_ind = j;
      }
    }
    denom_sum = 1.0;
    for (j=0; j<m; j++) {
      if (j==max_ind) 
          num[j] = 1.0;
      else {
          ff = log(mix_prop[j]) - log(mix_prop[max_ind]) + et*support[max_ind] - et*support[j];
          ff = ff + x*log(support[j]) - x*log(support[max_ind]);
          num[j] = exp(ff);
          denom_sum += num[j];
      }
    }
    for (j=0; j<m; j++) {
      post_prob[i + n*j] = num[j]/denom_sum;
    }
    *loglik += log(denom_sum) + log(mix_prop[max_ind]) - et*support[max_ind] - lgamma(x + 1) + x*log(et);
  }
}

