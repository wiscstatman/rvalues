\name{rvaluesMCMC}
\alias{rvaluesMCMC}
\title{
R-values from MCMC output.
}
\description{
Returns r-values from an array of MCMC output.
}
\usage{
rvaluesMCMC(output, qtheta, alpha.grid = NULL, ngrid = NULL, smooth = "none")
}
\arguments{
  \item{output}{
       a matrix contatining mcmc ouput. The ith row should
       represent a sample from the posterior of the ith parameter
       of interest. 
   }
   \item{qtheta}{
       either a function which returns the quantiles (for upper tail probs.) of 
       theta or a vector of theta-quantiles. 
   }
   \item{alpha.grid}{
       grid of values in (0,1); used for the discrete approximation
       approach for computing r-values.
   }
   \item{ngrid}{
       number of grid points for alpha.grid; only relevant when \code{alpha.grid = NULL}
   }
   \item{smooth}{
       either \code{smooth="none"} or \code{smooth} takes 
               a value between 0 and 10; this determines the level of smoothing applied to the
               estimate of \eqn{\lambda(\alpha)}; if \code{smooth} is given a number, the
               number is used as the \code{bass} argument in \code{\link{supsmu}}.
   }
}
\value{
An object of class "rvals" which is a list containing at least the following components:
\item{main}{a data frame containing the r-values, the r-value rankings along with
the rankings from several other common procedures}
\item{aux}{a list containing other extraneous information}
\item{rvalues}{a vector of r-values}
}
\references{
Henderson, N.C. and Newton, M.A. (2016). \emph{Making the cut: improved ranking and selection
for large-scale inference}. J. Royal Statist. Soc. B., 78(4), 781-804.
\url{https://doi.org/10.1111/rssb.12131}
\url{https://arxiv.org/abs/1312.5776}
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{rvalues}}, \code{\link{PostSummaries}}
}
\examples{
data(MCMCtest)

### For the MCMC output in MCMC_test, the prior assumed for the effect sizes of 
### interest was a mixture of two t-distributions. The function qthetaTMix 
### computes the quantiles for this prior.

qthetaTMix <- function(p) {
    ### function to compute quantiles (for upper tail probabilities) for a 
    ###  mixture of two t-distributions
    mu <- c(.35,-.12)
    sig <- c(.2,.08)
    mix.prop <- c(.25,.75)

    ff <- function(x,pp) {
        prob_less <- 0
        for(k in 1:2) {
            prob_less <- prob_less + pt((x - mu[k])/sig[k],df=4,lower.tail=FALSE)*mix.prop[k] 
        }
        return(prob_less - pp)
    }

    nn <- length(p)
    ans <- numeric(nn)
    for(i in 1:nn) {
        ans[i] <- uniroot(ff,interval=c(-5,5),tol=1e-6,pp=p[i])$root
    }
    return(ans)
}

rvs <- rvaluesMCMC(MCMCtest, qtheta = qthetaTMix)
}
\keyword{htest}
