\name{bcwest}
\alias{bcwest}
\docType{data}
\title{
Breast Cancer Gene Expression Data
}
\description{
Effect size estimates and standard errors obtained from 
gene expression measurements on 7129 genes across 49 samples.
}
\usage{data(bcwest)}
\format{
  A data frame with 7129 observations on the following 2 variables.
  \describe{
    \item{\code{estimates}}{a vector of effect size estimates}
    \item{\code{std.err}}{standard errors associated with effect size estimates}
  }
}
\details{
A description of the original data may be found in West et al. (2001).
For each gene, the effect size estimate was computed from the difference
in the mean expression levels of the two groups (i.e., 
mean(bc-positive group) - mean(bc-negative group)).
}
\source{
T. Hothorn, P. Buehlmann, T. Kneib, M. Schmid, and B. Hofner (2013). 
mboost: Model-Based Boosting, R package version 2.2-3, http://CRAN.R-project.org/package=mboost.
}
\references{
  Mike West, Carrie Blanchette, Holly Dressman, Erich Huang,
  Seiichi Ishida, Rainer Spang, Harry Zuzan, John A. Olson Jr.,
  Jeffrey R. Marks and Joseph R. Nevins (2001),
  Predicting the clinical status of human breast cancer by using
  gene expression profiles, \emph{Proceedings of the National Academy of Sciences},
  \bold{98}, 11462-11467.
}
\keyword{datasets}
