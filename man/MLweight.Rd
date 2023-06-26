\name{MLweight}
\alias{MLweight}

\title{
    Weights for Guttman Errors in two-level test data
}

\description{
    Computes weights for Guttman errors in two-level test data (Koopman et al., 2017)
}

\usage{
    MLweight(X, maxx = NULL, minx = NULL, itemstep.order = NULL)
}

\arguments{
  \item{X}{Data matrix with a subject column and two item columns. The subject column is assumed to be the first.}
  \item{maxx}{The highest possible answer category. If not specified it is determined by using the highest item score.}
  \item{minx}{The lowest possible answer category. If not specified it is determined by using the lowest item score.}
  \item{itemstep.order}{The two columns pertaining the two items in question, from a (possibly larger) matrix with 
		number of rows equal to the number of item steps (m) and number of columns equal to the number of items (J). 
		The matrix should consis the integers 1 : (m * J), indicating a predefined order of the items steps with respect to popularity. 
		Value 1 indicates the easiest (most popular) item step, value (m * J) indicates the most difficult item step.}
}

\value{
 Returns a vector with the weights for each item-score pattern of a given item-pair. In case of ties in item popularities the average weights across possible item-orderings are returned.
}

\references{
  
  Koopman, L., Van der Ark, L. A., & Zijlstra, B. J. H. (2017). 
  \emph{Weighted Guttman Errors: Handling Ties and Two-Level Data}. 
  In L. A. Van der Ark, S. Culpepper, J. A. Douglas, W.-C. Wang, & M. Wiberg (Eds.), 
  \emph{Quantitative Psychology: The 81st Annual Meeting of the Psychometric Society, Asheville, North Carolina, 2016} (pp. 183-190). 
  Springer. 
  \doi{10.1007/978-3-319-56294-0_17} 

}

\author{
   L. Koopman \email{V.E.C.Koopman@uva.nl}
   L. A. van der Ark \email{L.A.vanderArk@uva.nl} 
}

\seealso{
    \code{\link{MLcoefH}}
}


\examples{
    smallData <- data.frame(Subs = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3),
                        Xa = c(2, 0, 0, 1, 0, 2, 2, 0, 2, 2, 1, 2, 1, 2, 2), 
                        Xb = c(1, 1, 1, 0, 1, 2, 2, 1, 2, 2, 1, 0, 2, 2, 2), 
                        Xc = c(0, 0, 0, 1, 0, 2, 2, 1, 2, 1, 0, 0, 1, 1, 2))

    # Compute the weights
    I <- 3 
    for(i in 1:(I - 1)) for(j in (i + 1):I) 
    print(MLweight(smallData[, c(1, i + 1, j + 1)], minx = 0, maxx = 2))
}
\concept{weighted Guttman errors} 
\concept{two-level scalability coefficient}
