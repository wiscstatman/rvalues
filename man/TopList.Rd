\name{TopList}
\alias{TopList}
\title{
List of Top Units
}
\description{
Returns a list of the top units ranked according
to "r-value" or another specified statistic.
}
\usage{
TopList(object, topnum = 10, sorted.by = c("RValue","PostMean","MLE","PVal"))
}
\arguments{
  \item{object}{
An object of class "rvals"
}
  \item{topnum}{
The length of the top list.
}
  \item{sorted.by}{
The statistic by which to sort; this could be \code{sorted.by = "RValue"},
\code{sorted.by = "PostMean"}, \code{sorted.by = "MLE"}, or
\code{sorted.by = "PVal"}
}
}
\value{
a data frame with \code{topnum} rows and columns containing the
r-value, mle, posterior mean, and p-value rankings.
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{rvalues}}
}
\examples{
n <- 500
theta <- rnorm(n)
ses <- sqrt(rgamma(n,shape=1,scale=1))
XX <- theta + ses*rnorm(n)
dd <- cbind(XX,ses)

rvs <- rvalues(dd, family = gaussian)

TopList(rvs, topnum = 12)
TopList(rvs, topnum = 15, sorted.by = "MLE")

}
\keyword{ hstat }
\keyword{ models }
