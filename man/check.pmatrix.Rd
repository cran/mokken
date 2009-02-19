\name{check.pmatrix}
\alias{check.pmatrix}

\title{
Check of Nonintersection Using Method Pmatrix
}

\description{
 computes the P(++) and P(- -) matrix and auxiliary information
}

\usage{
check.pmatrix(X, minvi = 0.03)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
  \item{minvi}{minimum size of a violation that is reported}
}

\details{
The output is often numerous. 
Functions \code{plot} and  \code{summary} can be used to summarize the output. 
See Van der Ark (2007) for an example.
}

\value{
 \item{Ppp}{P(++) matrix}
 \item{Pmm}{P(--) matrix}
 \item{I.item}{vector indicating to which items the rows and column the P(++) matrix belong}
 \item{I.step}{the labels of the item steps in order of popularity}
 \item{I.labels}{the item labels}
 \item{Hi}{the item scalability coefficients Hi}
 \item{minvi}{the value of minvi}
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

\seealso{
\code{\link{coefH}}, \code{\link{check.iio}}, \code{\link{check.monotonicity}}, \code{\link{check.restscore}}, \code{\link{check.reliability}}, 
\code{\link{plot.pmatrix.class}}, \code{\link{summary.pmatrix.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
pmatrix.list <- check.pmatrix(Communality)
plot(pmatrix.list)
summary(pmatrix.list)
}

\keyword{univar}
