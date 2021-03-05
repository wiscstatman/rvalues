\name{Valpha}
\alias{Valpha}
\title{
R-values from a matrix of posterior tail probabilities.
}
\description{
Computes r-values directly from a "Valpha" matrix V where each
column of Valpha contains posterior tail probabilities relative
to a threshold indexed by alpha. 
}
\usage{
Valpha(V, alpha.grid, smooth = "none")
}
\arguments{
  \item{V}{a numeric vector with (i,j) entry: V[i,j] = P(theta_{i} >= theta[alpha_{j}]|data) }
  \item{alpha.grid}{grid of values in (0,1); used for the discrete approximation
       approach for computing r-values.}
  \item{smooth}{either \code{smooth="none"} or \code{smooth} takes 
               a value between 0 and 10; this determines the level of smoothing applied to the
               estimate of \eqn{\lambda(\alpha)}; if \code{smooth} is given a number, the
               number is used as the \code{bass} argument in \code{\link{supsmu}}.}
}
\value{
A list with the following components
\item{rvalues}{a vector of computed r-values}
\item{Vmarginals}{The estimated V-marginals along the alpha grid points}
\item{Vmarginals.smooth}{a function obtained through interpolation
and smoothing (if desired) the Vmarginals; i.e., an estimate of
\eqn{\lambda(\alpha)} (see \code{\link{rvalues}})}
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
\code{\link{rvalues}}  \code{\link{rvaluesMCMC}} 
}
\examples{
\dontrun{
data(fluEnrich)
rvobj <- rvalues(fluEnrich, family = binomial)

Vrvals <- Valpha(rvobj$aux$V, rvobj$aux$alpha.grid)
}
}
\keyword{ hstat}
