\name{mroot}
\alias{mroot}
\title{
Multi-dimensional Root (Zero) Finding
}
\description{
For a given multi-dimensional function with both a vector of lower bounds
and upper bounds, \code{mroot} finds a vector such that each 
component of the function is zero.
}
\usage{
mroot(f, lower, upper, \dots, f.lower = f(lower, ...), f.upper = f(upper, ...),
      tol = .Machine$double.eps^0.25, maxiter = 5000) 
}
\arguments{
  \item{f}{
the function for which the root is sought
}
  \item{lower}{
a vector of lower end points 
}
  \item{upper}{
a vector of upper end points
}
  \item{\dots}{
additional arguments to be passed to \code{f}
}
\item{f.lower, f.upper}{
the same as \code{f(lower)} and \code{f(upper)}
}
  \item{tol}{
the convergence tolerance
}
  \item{maxiter}{
the maximum number of iterations
}
}
\details{
The function \code{f} is from \eqn{R^{n}} to \eqn{R^{n}} with 
\eqn{f(x_1,\dots,x_n) = (f_1(x_1),\dots,f_n(x_n))}.

A root \eqn{x = (x_1,\dots,x_n)} of \eqn{f} satisfies
\eqn{f_k(x_k) = 0} for each component \eqn{k}.

\code{lower} \eqn{ = (l_1,\ldots,l_n)} and \code{upper} 
\eqn{= (u_1,\dots,u_n)} are both n-dimensional vectors
such that, for each \eqn{k}, \eqn{f_k} changes sign over the
interval \eqn{ [l_k, u_k]}. 
}
\value{
a vector giving the estimated root of the function
}
\author{
Nicholas Henderson
}
\seealso{
\code{\link{uniroot}}
}
\examples{
ff <- function(x,a) {
    ans <- qnorm(x) - a
    return(ans)
}
n <- 20000
a <- rnorm(n)
low <- rep(0,n)
up <- rep(1,n)

## Find the roots of ff, first using mroot and
## then by using uniroot inside a loop.

system.time(mr <- mroot(ff, lower = low, upper = up, a = a))

ur <- rep(0,n)
system.time({
for(i in 1:n) {
   ur[i] <- uniroot(ff, lower = 0, upper = 1, a = a[i])$root    
}
})

}
\keyword{optimize}
