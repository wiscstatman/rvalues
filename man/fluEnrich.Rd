\name{fluEnrich}
\alias{fluEnrich}
\docType{data}
\title{
Flu Enrichment Data
}
\description{
Gene-set enrichment for genes that have been identified as
having an effect on influenza-virus replication.
}
\usage{data(fluEnrich)}
\format{
  A data frame with 5719 observations on the following 3 variables.
  \describe{
    \item{\code{nflugenes}}{number of genes both annotated to the given GO term and in the collection of flu genes}
    \item{\code{setsize}}{number of genes annotated to the given GO term}
    \item{\code{GO_terms}}{the GO (gene ontology) term label} 
  }
}
\details{
These data were produced by associating the 984 genes (the collection of flu genes)
identified in the Hao et al. (2013) meta-analysis with gene ontology (GO) gene sets (GO terms). 
In total, 17959 human genes were annotated to at least one GO term and 16572
GO terms were available, though this data set only contains the 5719 terms 
which annotated between 10 and 1000 human genes.
}
\references{
Hao, L. Q. He, Z. Wang, M. Craven, M. A. Newton, and P. Ahlquist (2013). 
Limited agreement of independent RNAi screens for virus-required host genes
owes more to false-negatives than false-positive factors.
\emph{PLoS computational biology}, \bold{9}, 9, e1003235.
}
\keyword{datasets}
