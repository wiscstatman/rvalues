\name{PostPercentile}
\alias{PostPercentile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior expected percentiles
}
\description{
Computes posterior expected percentiles for both parametric
and nonparametric models.
}
\usage{
PostPercentile(object)
}
\arguments{
  \item{object}{
An object of class "rvals"
}
}
\details{
With parameters of interest \eqn{\theta_1,...,\theta_n} the rank of 
the ith parameter (when we set the ranking so that the largest \eqn{\theta_i}
gets rank 1) is defined as \eqn{ rank(\theta_i) = sum_j(\theta_j \ge \theta_i)} 
and the associated percentile is 
\eqn{ perc(\theta_i) = rank(\theta_i)/(n+1).} The posterior expected percentile 
for the ith unit (see e.g., Lin et. al. (2006)) is simply
the expected value of \eqn{perc(\theta_i)} given the data.

The function \code{PostPercentile} computes an asymptotic version of the 
posterior expected percentile, which is defined as 
\deqn{ P(\theta_i \le \theta|data),}
where \eqn{\theta} has the same distribution as \eqn{\theta_i} and is 
independent of both \eqn{\theta_i} and the data.  
See Henderson and Newton (2014) for additional details.
}
\value{
A vector of estimated posterior expected percentiles.
}
\references{
Henderson, N.C. and Newton, M.A. (2014) Making the Cut: Improved Ranking and Selection 
for Large-Scale Inference. 
\url{http://arxiv.org/abs/1312.5776}

Lin, R., Louis, T.A., Paddock, S.M., and Ridgeway, G. (2006). Loss function based
ranking in two-stage, hierarchical models. \emph{Bayesian Analysis}, \bold{1}, 
915--946.
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{rvalues}}
}
\examples{
n <- 3000
theta <- rnorm(n, sd = 3)
ses <- sqrt(rgamma(n, shape = 1, scale = 1))
XX <- theta + ses*rnorm(n)
dd <- cbind(XX,ses)

rv <- rvalues(dd, family = gaussian)

perc <- PostPercentile(rv)
plot(rv$rvalues, perc)
}
\keyword{ hstat }
\keyword{ models }
