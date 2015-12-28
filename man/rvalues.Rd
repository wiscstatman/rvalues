\name{rvalues}
\alias{rvalues}
\title{
R-values
}
\description{
Given data on a collection of units, this
function computes
r-values which are percentiles constructed
to maximize the agreement between the reported
percentiles and the percentiles of the effect of
interest. Additional details about r-values are provided below and can also be found
in the listed references.
}
\usage{
rvalues(data, family = gaussian, hypers = "estimate", prior = "conjugate",
       alpha.grid = NULL, ngrid = NULL, smooth = "none", control = list()) 
}
\arguments{
  \item{data}{A data frame or a matrix with the number of rows equal to the number of sampling units. 
              The first column should contain the main estimates, and the second column should contain 
              the nuisance terms.}
  \item{family}{An argument which determines the sampling distribution; this could be either
                \code{family = gaussian}, \code{family = tdist}, \code{family = binomial},
                 \code{family = poisson}}
  \item{hypers}{values of the hyperparameters; only meaningful when the conjugate prior is used; if set to "estimate", 
  the hyperparameters are found through maximum likelihood; if not set to "estimate" the user
  should supply a vector of length two.}
  \item{prior}{the form of the prior; either \code{prior="conjugate"} or \code{prior="nonparametric"}.}
  \item{alpha.grid}{a numeric vector of points in (0,1); this grid is used in the 
                    discrete approximation of r-values}
  \item{ngrid}{number of grid points for alpha.grid; only relevant when \code{alpha.grid=NULL}}
  \item{smooth}{either \code{smooth="none"} or \code{smooth} takes 
               a value between 0 and 10; this determines the level of smoothing applied to the
               estimate of \eqn{\lambda(\alpha)} (see below for the definition of
               \eqn{\lambda(\alpha)}); if \code{smooth} is given a number, the
               number is used as the \code{bass} argument in \code{\link{supsmu}}.}
  \item{control}{a list of control parameters for estimation of the prior; only used 
                 when the prior is nonparametric}
}
\details{
  The r-value computation assumes the following two-level sampling model
  \eqn{ X_i|\theta_i} ~ \eqn{p(x|\theta_i,\eta_i)}
  and \eqn{\theta_i} ~ \eqn{F}, for \eqn{i = 1,...,n},
  with parameters of interest \eqn{\theta_i}, effect size estimates \eqn{X_i},
  and nuisance terms \eqn{\eta_i}. The form of \eqn{p(x|\theta_i,\eta_i)} is determined
  by the \code{family} argument. When \code{family = gaussian}, it is assumed that
  \eqn{X_i|\theta_i,\eta_i} ~ N(\eqn{\theta_i,\eta_i^{2})}.
  When \code{family = binomial}, the \eqn{(X_i,\eta_i)} represent the number of successes
  and number of trials respectively, and it is assumed that \eqn{X_i|\theta_i,\eta_i} ~
  Binomial\eqn{(\theta_i,\eta_i)}. When \code{family = poisson}, the \eqn{{X_i}} should be
  counts, and it is assumed that \eqn{X_i|\theta_i,\eta_i} ~ Poisson(\eqn{\theta_i * \eta_i)}. 
  
  The distribution of the effect sizes \eqn{F} may be a parametric distribution
  that is conjugate to the corresponding \code{family} argument,
  or it may be estimated nonparametrically. When it is desired that \eqn{F} be
  parametric (i.e., \code{prior = "conjugate"}), the default is to estimate the 
  hyperparameters (i.e., \code{hypers = "estimate"}), but these may be supplied by the 
  user as a vector of length two. To estimate \eqn{F} nonparametrically, one
  should use \code{prior = "nonparametric"} (see \code{\link{npmle}} for
  further details about nonparametric estimation of \eqn{F}).
  
  The \emph{r-value}, \eqn{r_i}, assigned to the ith case of interest is determined by
  \eqn{  r_i = } inf[ \eqn{0 < \alpha < 1: V_\alpha(X_i,\eta_i) \ge \lambda(\alpha) } ]
   where \eqn{V_\alpha(X_i,\eta_i) = P( \theta_i \ge \theta_\alpha|X_i,\eta_i) }
   is the posterior probability that \eqn{\theta_i} exceeds the threshold \eqn{\theta_\alpha},
   and \eqn{\lambda(\alpha)} is the upper-\eqn{\alpha}th quantile associated
   with the marginal distribution of \eqn{V_\alpha(X_i,\eta_i)} (i.e., 
   \eqn{ P(V_\alpha(X_i,\eta_i) \ge \lambda(\alpha)) = \alpha). } Similarly, 
   the threshold \eqn{\theta_\alpha} is the upper-\eqn{\alpha}th quantile of 
   \eqn{F} (i.e., \eqn{P(\theta_i \ge \theta_\alpha) = \alpha} ). 
}
\value{
An object of class "rvals" which is a list containing at least the following components:
\item{main}{a data frame containing the r-values, the r-value rankings along with
the rankings from several other common procedures}
\item{aux}{a list containing other extraneous information}
\item{rvalues}{a vector of r-values}
}
\references{
Henderson, N.C. and Newton, M.A. (2015) \emph{Making the Cut: Improved Ranking and Selection 
for Large-Scale Inference.} 
\url{http://arxiv.org/abs/1312.5776}
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{rvaluesMCMC}}, \code{\link{PostSummaries}}, \code{\link{Valpha}}
}
\examples{
\dontrun{
### Binomial example with Beta prior:
data(fluEnrich)
flu.rvals <- rvalues(fluEnrich, family = binomial)
hist(flu.rvals$rvalues)

### look at the r-values for indices 10 and 2484
fig_indices  <- c(10,2484)
fluEnrich[fig_indices,]

flu.rvals$rvalues[fig_indices]

### Gaussian sampling distribution with nonparametric prior
### Use a maximum of 5 iterations for the nonparam. estimate
data(hiv)
hiv.rvals <- rvalues(hiv, prior = "nonparametric")
}
}
\keyword{hstat}
\keyword{models}
