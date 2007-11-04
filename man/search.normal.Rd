\name{search.normal}
\alias{search.normal}

\title{
Automated Item Selection Algorithm for Mokken Scale Analysis
}

\description{
Returns a vector with as many elements as there are items, indicating the scale an item belongs to
}
\usage{
search.normal(X, lowerbound = 0.3, alpha = 0.05)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
  \item{lowerbound}{numeric scaling criterium; 0 <= \code{lowerbound} < 1}
  \item{alpha}{Type I error level}
}
\details{
The number of Mokken scales cannot exceed \code{ncol(X)/2}.
Procedure may be slow for large data sets.
}

\value{An indicator vector of length J. 
       Each entry refers to an item. 
       Items with same integer belong to the same Mokken scale.
       A zero indicates an unscalable item.
       If \emph{n} is the largest integer, then \emph{n} Mokken scales were found.
}

\references{
  Mokken, R. J. (1971) 
  \emph{A Theory and Procedure of Scale Analysis}. 
  Berlin, Germany: De Gruyter.
  
  Molenaar, I.W. and Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  Groningen, The Netherlands: IEC ProGAMMA.
  
  Sijtsma, K, and Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Thousand Oaks, CA: Sage.
  
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{http://www.jstatsoft.org}
}
\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}

\note{
}

\seealso{
\code{\link{coefH}}, \code{\link{check.monotonicity}}, \code{\link{check.pmatrix}}, \code{\link{check.restscore}} 
}

\examples{
data(acl)
# Partition all 212 items into mokken scales (may take some time).
scale <- search.normal(acl)      

# investigate monotonicity for all items in the first scale.
monotonicity.list <- check.monotonicity(acl[,scale==1])

# summary of the results
summary(monotonicity.list)
}
\keyword{cluster}


