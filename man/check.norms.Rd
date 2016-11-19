\name{check.norms}
\alias{check.norms}

\title{
Standard errors for norm statistics
}

\description{
The function presents standard errors for the mean, standard deviation, standard scores, stanine boundaries, and percentiles based on a vector of test scores
(Oosterhuis, Van der Ark, and Sijtsma, in press).  
}
\usage{
check.norms(y, nice.output = TRUE)
}

\arguments{
  \item{y}{numerical vector. Typically a numerical vector of length N, representing the test scores of N respondents. Missing values are not allowed}
  \item{nice.output}{Logical: If \code{TRUE}, norm statistics and standard errors are combined in an a single object of class \code{noquote}}
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
  Oosterhuis, H. E. M., Van der Ark, L. A., and Sijtsma, K. (in press).
  Computing standard errors and confidence intervals for norm statistics.
  \emph{Psychometrika}. doi:10.1007/s11336-016-9535-8
}

\seealso{
 \code{\link{DS14}}, \code{\link{recode}}, \code{\link{twoway}}
}

\examples{
  data(DS14)

  # Handle missing data and recode negatively worded items
  X <- DS14[, 3 : 16]
  X <- twoway(X)
  X <- recode(X, c(1, 3))

  # Negative affectivity
  Na <- X[, c(1, 3, 6, 8, 10, 11, 14)] 

  # Social inhibition
  Si <- X[, c(2, 4, 5, 7, 9, 12, 13)] 

  # Norms
  check.norms(rowSums(Na))
  check.norms(rowSums(Si))
}
\keyword{check.norms}
