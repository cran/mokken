\name{check.errors}
\alias{check.errors}

\title{
Check the number of Guttman errors for each respondent
}

\description{
Returns an integer vector with the number of Guttman errors (Guttman, 1944) for each respondent (e.g., Meijer, 1994; Mokken, 1971; Molenaar and Sijtsma, 2000; Sijtsma and Moelnaar, 2002). 
}
\usage{
check.errors(X)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items.  
          Missing values are not allowed}
}


\value{
  Vector containing the number of Guttman errors for each respondent
}

\references{
  Guttman, L. (1944) A basis for scaling qualitative data.
  \emph{American Sociological Review, 9}, 139-150.
  
  Meijer, R. R. (1994)
  The number of Guttman errors as a simple and powerful person-fit statistic.
  \emph{Applied Psychological Measurement, 18}, 311-314.

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
  \url{http://www.jstatsoft.org/v20/i11}
}

\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}

\seealso{
\code{\link{check.iio}}, 
\code{\link{check.monotonicity}}, 
\code{\link{check.pmatrix}}, 
\code{\link{check.reliability}} 
\code{\link{coefH}}, 
\code{\link{plot.restscore.class}}, 
\code{\link{summary.restscore.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
Guttman.errors <- check.errors(Communality)
hist(Guttman.errors, breaks = 0:max(Guttman.errors))
}
\keyword{restscore}
