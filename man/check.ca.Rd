\name{check.ca}
\alias{check.ca}

\title{
Check conditional association to indentify local dependence.
}

\description{
The function uses three special cases of conditional association (CA; Holland and Rosenbaum, 1986) to identify positive and negative local dependence in 
the monotone homogeneity model. Straat, Van der Ark, and Sijtsma (2016; also, see Sijtsma, Van der Ark, and Straat, 2015) described the procedure.  
}
\usage{
check.ca(X, Windex = FALSE, MINSIZE = 4, NWEIGHTOPTION = "noweight", 
            COVWEIGHTOPTION = "pnorm", MINGROUP = 4)
}

\arguments{
  \item{X}{Matrix of integers, missing values are not allowed}
  \item{Windex}{Boolean. Should output contain indices W1, W2, and W3?}
  \item{MINSIZE}{Minimum sample size of a rest-score group}
  \item{NWEIGHTOPTION}{Weight of sample size on each conditional covariance. Options: "noweight" (each covariance has weight 1, default in Straat et al., 2016) and "sqrt" (each covariance has weight sqrt(N_k(x)), this option was used in an older, decrepit,  version of Straat et al., 2016)}
  \item{COVWEIGHTOPTION}{Weight of each conditional covariance on the computation of W1, W2, and W3. Options: "pnorm" (weight equals P[cov < 0], default in Straat et al., 2014) and "noweight" (if cov < 0, then weight = 1, and weight = 0, otherwise; this option was used in a previous version of Straat et al., 2014)}
  \item{MINGROUP}{Minimum sample size of the conditioning variable to compute a covariance. Since the term N-3 is used in the computation of the standard error, N = 4 is the default.}
}

\value{
 list of three components: 
 
 (1) \code{InScale} (vector of booleans with length equal to the number of items): indicates whether an item is still in the scale.
 
 (2) \code{Index} (list): Numerical values of indices W1, W2, and W3 (shown only if \code{Windex = TRUE}). \code{Index} has three subcomponents: \code{W1}, \code{W2}, and \code{W3}.
 
 (3) \code{Flagged} (list): Boolean indicating whether a value of W1, W2, and W3 is flagged (1) or not (0) (shown only if \code{Windex = TRUE}) \code{Index} has three subcomponents: \code{F1}, \code{F2}, and \code{F3}, corresponding to the flagging of indices W1, W2, and W3, respectively.
 
 Subcomponents correspond to the iteration. The first subcomponent refers to the situation with all items in the test, the second subcomponent refers to the sitution with the worst item deleted, the third subcomponent refers to the sitution with the two worst items deleted, etc.
}


\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl} and  J. H. Straat 
}

\references{
 Straat, J. H., van der Ark, L. A., and Sijtsma, K. (2016).
 Using conditional association to identify local dependence.
 Manuscript submitted for publication.

 Sijtsma, K., van der Ark, L. A., and Straat, J. H. (2015) 
 Goodness of fit methods for nonparametric IRT models.
 In L. A. van der Ark, D. M. Bolt, W.-C. Wang, J. Douglas, and S.-M. Chow (Eds.),
  \emph{Quantitative psychology research: The 79th Annual Meeting of the Psychometric Society, Madison, Wisconsin, 2014.} (pp. 109 - 120)
  New York, NY: Springer.

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
  
  check.ca(Na, TRUE)
}
\keyword{check.ca}
