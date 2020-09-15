\name{check.errors}
\alias{check.errors}

\title{
Check the number of Guttman errors (Gplus) and the number of infrequent scores (Oplus) for each respondent
}

\description{
Returns a lost containing outlier scores Gplus (number of Guttman errors; Guttman, 1944) and Oplus for each respondent (Zijlstra, van der Ark & Sijtsma, 2007). 
}
\usage{
check.errors(X, returnGplus = TRUE, returnOplus = FALSE)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items.  
          Missing values are not allowed}
  \item{returnGplus}{Boolean. If \code{TRUE} the output contains outlier score Gplus}        
  \item{returnOplus}{Boolean. If \code{TRUE} the output contains outlier score Oplus}        
}


\value{
  List. Depending on the values of \code{returnGplus} and \code{returnOplus}, the output contains outlier score Gplus (the number of Guttman errors) 
  and Oplusfor each respondent
}

\references{
  Guttman, L. (1944) A basis for scaling qualitative data.
  \emph{American Sociological Review, 9}, 139-150.
    
  Meijer, R. R. (1994)
  The number of Guttman errors as a simple and powerful person-fit statistic.
  \emph{Applied Psychological Measurement, 18}, 311-314.
  \url{https://doi.org/10.1177/014662169401800402}


  Mokken, R. J. (1971) 
  \emph{A Theory and Procedure of Scale Analysis}. 
  De Gruyter.
  
   Molenaar, I.W., & Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  IEC ProGAMMA.
  
  Sijtsma, K., & Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Sage.
 
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{https://www.jstatsoft.org/article/view/v020i11}

  Zijlstra, W. P.,  Van der Ark, L. A., & Sijtsma, K. (e2007).
  Outlier detection in test and questionnaire data.
  \emph{Multivariate Behavioral Research, 42}, 531-555.
  \url{https://doi.org/10.1080/00273170701384340}
}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}

\seealso{
\code{\link{check.ca}}, 
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
Gplus <- check.errors(Communality, TRUE, FALSE)$Gplus
Oplus <- check.errors(Communality, FALSE, TRUE)$Oplus

hist(Gplus, breaks = 0:max(Gplus))
}
\concept{restscore}
