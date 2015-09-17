\name{tdist}
\alias{tdist}
\title{
t-distribution family object
}
\description{
A t-distribution family object which allows one to specify
a t-density for the sampling distribution. Modeled after \code{\link{family}} objects
often used in the glm function. 
}
\usage{
tdist(df)
}
\arguments{
  \item{df}{
  vector containing the degrees of freedom
}
}
\value{
An object of class \code{"newfam"}, which is a list
containing the following components
\item{family }{The family name}
\item{df }{The degrees of freedom}
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{\code{\link{family}}, \code{\link{glm}}, \code{\link{npmle}} }
\examples{
a <- tdist(df=5)
}
\keyword{stats}
