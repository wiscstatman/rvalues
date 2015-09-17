\name{batavgs}
\alias{batavgs}
\docType{data}
\title{
Batting Averages Data
}
\description{
Data set containing number of at-bats and number of hits for Major League baseball players over
the 2005 season.
}
\usage{data(batavgs)}
\format{
  A data frame with 929 observations on the following 7 variables.
  \describe{
    \item{\code{First.Name}}{factor; player's last name}
    \item{\code{Last.Name}}{factor; player's first name}
    \item{\code{Pitcher}}{numeric vector; an indicator of whether or not the player is a pitcher}
    \item{\code{midseasonAB}}{numeric vector; number of at-bats during the first half of the season}
    \item{\code{midseasonH}}{numeric vector; number of hits during the first half of the season}
    \item{\code{TotalAB}}{numeric vector; total number of at-bats over the season}
    \item{\code{TotalH}}{numeric vector; total number of hits over the season}
  }
}
\details{
The 2005 Major League Baseball season was roughly six months starting from the beginning
of April and ending at the beginning of October. Data from postseason play is not 
included. The midseason data were obtained by only considering the first three months
of the season.
}
\source{
http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.aoas/1206367815
}
\references{
Brown, L. D. (2008), In-Season prediction of batting averages: a field test
of empirical Bayes and Bayes Methodologies, 
\emph{The Annals of Applied Statistics}, \bold{2}, 1, 113--152.
}
\examples{
data(batavgs)
head(batavgs)
}
\keyword{datasets}
