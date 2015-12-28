\name{rvalueBoot}
\alias{rvalueBoot}
\title{
Bootstrapped r-values
}
\description{
Estimates a new prior for each bootstrap replications ... (need to add) 
}
\usage{
rvalueBoot(object, statistic = median, R, type = "nonparametric")
}
\arguments{
  \item{object}{
An object of class "rvals"
}
  \item{statistic}{
The statistic used to summarize the bootstrap replicates.
}
  \item{R}{
Number of bootstrap replicates
}
  \item{type}{
Either \code{type="nonparametric"} or \code{type="parametric"}; the nonparametric type
corresponds to the usual bootstrap where units are sampled with replacement.
}
}
\details{
When \code{type="nonparametric"}, the prior is re-estimated (using the resampled data) 
in each bootstrap replication, and r-values are re-computed with respect to this new model.

When \code{type="parametric"}, 
}
\value{
A list with the following two components
  \item{rval.repmat}{ A matrix where each column corresponds to a separate bootstrap replication.    }
  \item{rval.boot}{A vector of r-values obtained by applying the statistic to each row
                   of \code{rval.repmat}.}
}
\references{
Henderson, N.C. and Newton, M.A. (2015) Making the Cut: Improved Ranking and Selection 
for Large-Scale Inference. 
\url{http://arxiv.org/abs/1312.5776}
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{rvalues}}
}
\examples{
\dontrun{
n <- 3000
theta <- rnorm(n, sd = 3)
ses <- sqrt(rgamma(n, shape = 10, rate = 1))
XX <- theta + ses*rnorm(n)
dd <- cbind(XX,ses)

rv <- rvalues(dd, family = gaussian, prior = "conjugate")

rvb <- rvalueBoot(rv, R = 10) 
summary(rvb$rval.repmat[512,])
}
}
\keyword{ htest }
\keyword{ models }