\name{check.norms}
\alias{check.norms}

\title{
Standard errors for norm statistics
}

\description{
The function presents standard errors for the mean, standard deviation, standard scores, stanine boundaries, and percentiles based on a vector of test scores
(Oosterhuis, Van der Ark, and Sijtsma, 2016).  
}
\usage{
check.norms(y)
}

\arguments{
  \item{y}{numerical vector. Typically a numerical vector of length N, representing the test scores of N respondents. Missing values are not allowed}
}

\value{
 list of five components: 
 
 (1) \code{mean}: Sample mean and its standard error (\code{noquote}).
 
 (2) \code{sd}: Sample standard deviation and its standard error (\code{noquote}).
 
 (3) \code{z}: For each unique testscore, the test score, its frequency, the corresponding estimated standard score and its standard error (\code{noquote}).

 (4) \code{sta9}: The estimates of the 8 boundaries of the stanines and their standard error (\code{noquote}).
 
 (5) \code{z}: For each unique testscore, the test score, its frequency, the corresponding estimated percentile rank and its standard error (\code{noquote}).
 
}


\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl} and  H. E. M. Oosterhuis 
}

\references{
 Oosterhuis, H. E. M., van der Ark, L. A., and Sijtsma, K. (2016).
 Standard errors and confidence intervals of norm statistics for educational and psychological tests.
 Manuscript submitted for publication.

 Straat, J. H., van der Ark, L. A., and Sijtsma, K. (2016).
 Using conditional association to identify local dependence.
 Manuscript submitted for publication.
}

\examples{
# Data: A bootstrap sample N = 3111 from the data used by Straat et al. (2016)
data(DS14B)

# Negative affectivity
Na <- DS14B[, 1 : 7] 

# Social inhibition
Si <- DS14B[, 8 : 14]

check.norms(rowSums(Na))
check.norms(rowSums(Si))
}
\keyword{check.norms}
