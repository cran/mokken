\name{check.bounds}
\alias{check.bounds}

\title{
Check the relative lower bound for scalability coefficients
}

\description{
Returns the relative bounds for Mokken's scalability coefficients for dichotomous items as described by Ellis (2014). 
}
\usage{
check.bounds(X, quant = .90, lower = TRUE, upper = FALSE)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items.  
          Missing values are not allowed}
  \item{quant}{numerical value between 0 and 1 used for the computation of lower bound \code{L2rij}.
          The computation deviates somewhat from the proposal in Ellis (2014) because the \code{stats} function \code{quantile} is used.}
  \item{lower}{Boolean: If \code{TRUE}, the lower bounds are given.}        
  \item{upper}{Boolean: If \code{TRUE}, the upper bounds are given.}        
}


\value{
  List containing two lists \code{UpperBounds} and \code{LowerBounds}, each containing a list of two J x J matrices (J = number of items): 
  \code{L1rij} (overestimator of the lower bound for the correlation) and
  \code{L2rij} (underestimator of the lower bound for the correlation).
}

\references{
  Ellis, J. L. (2014). An inequality for correlations in unidimensional monotone latent variable models for binary variables.
  \emph{Psychometrika, 79}, 303-316.
}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}

\examples{
data(acl)
Communality <- acl[,1:10]
R <- cor(Communality)
res <- check.bounds(Communality, upper = TRUE)
L1rij <- res$LowerBounds$L1rij
L2rij <- res$LowerBounds$L2rij
U1rij <- res$UpperBounds$U1rij
U2rij <- res$UpperBounds$U2rij
# Correlations that meet L1rij (possibly an overestimation of the lower bound). 
R >= L1rij 
# Correlations that meet U1rij (possibly an overestimation of the upper bound). 
R <= U1rij 
# Correlations that meet L2rij (possibly an underestimation of the lower bound). 
R >= L2rij 
# Correlations that meet U2rij (possibly an underestimation of the upper bound). 
R <= U2rij 
}
\keyword{bounds}
