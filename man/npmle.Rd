\name{npmle}
\alias{npmle}
\title{
Maximum Likelihood Estimate of a Mixing Distribution.
}
\description{
Estimates the mixture distribution nonparametrically using an EM algorithm.
The estimate is discrete with the results being returned as a vector of support
points and a vector of associated mixture probabilities.
The available choices for the sampling distribution include: Normal, Poisson, Binomial
and t-distributions.
}
\usage{
npmle(data, family = gaussian, maxiter = 500, tol = 1e-4,
      smooth = TRUE, bass = 0, nmix = NULL)
}
\arguments{
  \item{data}{
      A data frame or a matrix with the number of rows equal
      to the number of sampling units. The first column should
      contain the main estimates, and the second column should
      contain the nuisance terms.
  }
  \item{family}{
     family determining the sampling distribution (see \cite{\link{family}})
  }
  \item{maxiter}{
     the maximum number of EM iterations
  }
  \item{tol}{
     the convergence tolerance
  }
  \item{smooth}{
     logical; whether or not to smooth the estimated cdf
  }
  \item{bass}{
     controls the smoothness level; only relevant if \code{smooth=TRUE}. 
     Values of up to 10 indicate increasing smoothness.
  }
  \item{nmix}{
      optional; the number of mixture components
  }
}
\details{
Assuming the following two-level sampling model
\eqn{ X_i|\theta_i} ~ \eqn{p(x|\theta_i,\eta_i)}
and \eqn{\theta_i} ~ \eqn{F} for \eqn{i = 1,...,n}.
The function \code{npmle} seeks to find an estimate of the mixing distribution
\eqn{F} which maximizes the marginal log-likelihood 
\deqn{
l(F) = \sum_i \int p( X_i |\theta, \eta_i) dF(\theta).
}
The distribution function maximizing \eqn{l(F)} is
known to be discrete; and thus, the estimated mixture distribution is
returned as a set of support points and associated mixture
probabilities.
}
\value{
   An object of class npmix which is a list containing at least 
   the following components 
   \item{ support }{a vector of estimated support points}
   \item{ mix.prop }{a vector of estimated mixture proportions}
   \item{ Fhat }{ a function; obtained through interpolation of
   the estimated discrete cdf}
   \item{ fhat }{ a function; estimate of the mixture density}
   \item{ loglik }{ value of the log-likelihood at each iteration}
   \item{ convergence }{ 0 indicates convergence; 1 indicates
   that convergence was not achieved}
   \item{ numiter }{the number of EM iterations required}
}
\references{Laird, N.M. (1978), Nonparametric maximum likelihood estimation of a 
mixing distribution, \emph{Journal of the American Statistical Association}, \bold{73}, 805--811. 

Lindsay, B.G. (1983), The geometry of mixture likelihoods: a general theory. \emph{The
Annals of Statistics}, \bold{11}, 86--94
}
\author{
   Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{npmixapply}}
}
\examples{
\dontrun{
data(hiv)
npobj <- npmle(hiv, family = tdist(df=6), maxiter = 25)


###  Generate Binomial data with Beta mixing distribution
n <- 3000
theta <- rbeta(n, shape1 = 2, shape2 = 10)
ntrials <- rpois(n, lambda = 10)
x <- rbinom(n, size = ntrials, prob = theta)

###  Estimate mixing distribution 
dd <- cbind(x,ntrials)
npest <- npmle(dd, family = binomial, maxiter = 25)

### compare with true mixture cdf
tt <- seq(1e-4,1 - 1e-4, by = .001)
plot(npest, lwd = 2)
lines(tt, pbeta(tt, shape1 = 2, shape2 = 10), lwd = 2, lty = 2)
}
}
\keyword{ nonparametric }
\keyword{ models }
