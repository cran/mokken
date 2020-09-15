\name{ICC}
\alias{ICC}

\title{
    Intraclass correlation
}

\description{
    Computes to intraclass correlation for multilevel data (ICC; 
    Snijders & Bosker, 1999, p. 17) for each item and the total scale of 
    a questionnaire (Koopman, Van der Ark, & Zijlstra, 2020) and the F-test 
    for the null hypothesis that the (total scale) ICC is zero (Snijders & Bosker, 2012, p. 22)
}

\usage{
    ICC(X)
}

\arguments{
  \item{X}{matrix or data frame of numeric data containing a grouping column and the
           item scores of \code{nrow(X)} respondents to \code{ncol(X) - 1} items. 
           Missing values are not allowed}
}

\value{
 \item{itemICC}{The ICC per item}  
 \item{scaleICC}{The ICC for the total scale and the corresponding results for the F-test}
}

\references{
  Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020). 
  \emph{A two-step procedure for scaling multilevel data using Mokken's scalability coefficients}. 
  Manuscript submitted for publication.

  Snijders, T. A. B., & Bosker, R. J. (2012). 
  \emph{Multilevel analysis: An introduction to basic and advanced multilevel modeling} (2nd ed.). Sage.
}

\author{
   Letty Koopman \email{V.E.C.Koopman@uva.nl}
   L. Andries van der Ark \email{L.A.vanderArk@uva.nl} 
}

\seealso{
\code{\link{SWMD}},
\code{\link{MLcoefH}}, 
}


\examples{

# Data example (Koopman et al., 2020)
data(SWMD)

# Compute ICC
ICC(SWMD)

}
\concept{intraclass correlation}
\concept{ICC}
