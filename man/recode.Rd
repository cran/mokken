\name{recode}
\alias{recode}

\title{
Recodes negatively worded items
}

\description{
Returns a matrix or data.frame with the indicated items recoded.
}
\usage{
recode(X, items = NULL, values = defaultValues)
}

\details{
 The result is \code{X} for which columns \code{items} have been recoded.      
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
           containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items.  
           Missing values are allowed}
  \item{items}{Vector of integers indicating the items to be recoded}         
  \item{values}{Vector of possible item scores. By default the range of the observed values is taken} 
}

\value{
The result is \code{X} for which columns \code{items} have been recoded.      
}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}

\seealso{
 \code{\link{DS14}}, \code{\link{twoway}}
}

\examples{
  data(DS14)
  
  # Handle missing data and recode negatively worded items
  X <- DS14[, 3 : 16]
  X <- twoway(X)
  X <- recode(X, c(1, 3))
  head(X)
}

\concept{recode}
