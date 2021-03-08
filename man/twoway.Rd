\name{twoway}
\alias{twoway}

\title{
Two-way imputation
}

\description{
Returns a single or multiple completed data sets using two-way imputation with normally distributed errors.
}
\usage{
twoway(X, nCompletedDataSets = 1, minX = defaultMinX, maxX = defaultMaxX, seed = FALSE)
}

\details{
 For single imputation (\code{nCompletedDataSets == 1}, default) the function returns an object of the same class as \code{X}, 
 for multiple imputation (\code{nCompletedDataSets > 1}) the function returns a list.
 References for two-way imputation include Bernaards and Sijtsma (2000), Sijtsma and Van der Ark (2003), 
 and Van Ginkel, Van der Ark, and Sijtsma (2007).       
}

\arguments{
  \item{X}{matrix or data frame of integer data 
           containing the score of \code{now(X)} respondents to \code{nicol(X)} items.  
           Typically \code{X} contains missing values.}
  \item{nCompletedDataSets}{Number of completed data sets.}         
  \item{minX}{Minimum item score. By default, the minimum item score is the lowest score found in the data.} 
  \item{maxX}{Maximum item score. By default, the maximum item score is the highest score found in the data.} 
  \item{seed}{Seed for random sampling. If \code{seed = FALSE} (default), no seed is given, otherwise \code{seed} must be 
              a numeric value. Replications having the same seed result in exactly the same outcome value.} 
}

\value{
 The result is \code{X} for which the missing values have been replaced by imputed values. For multiple 
 imputations, the result is a list of matrices/data frames. For single  
 imputations, the result is a matrix/data frame.      
}

\references{

  Bernaards, C. A., & Sijtsma, K. (2000).  
  Influence of simple imputation and EM methods on factor analysis when item nonresponse in questionnaire data is nonignorable
  \emph{Multivariate Behavioral Research, 35}, 321-364.
  \doi{10.1207/S15327906MBR3503_03}

  Sijtsma, K., & Van der Ark, L. A. (2003).
  Investigation and treatment of missing item scores in test and questionnaire data.
  \emph{Multivariate Behavioral Research, 38}, 505-528.
  \doi{10.1207/s15327906mbr3804_4}

  Van Ginkel, J. R., Van dec Ark, L. A., & Sijtsma, K. (2007).
  Multiple imputation of item scores in test and questionnaire data, and influence on psychometric results.
  \emph{Multivariate aBehavioral Research, 42}, 387-414.
  \doi{10.1080/00273170701360803}
}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}

\seealso{\code{\link{DS14}}, \code{\link{recode}}}

\examples{
  data(DS14)

  # Handle missing data and recode negatively worded items
  X <- DS14[, 3 : 16]
  X <- twoway(X)
  X <- recode(X, c(1, 3))
  head(X)
}

\concept{twoway}
