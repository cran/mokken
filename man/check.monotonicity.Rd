\name{check.monotonicity}
\alias{check.monotonicity}

\title{
Check of Monotonicity
}

\description{
Returns a list (of class monotonicity.class) with results from the investigation of monotonicity.
}

\usage{
check.monotonicity(X, minvi = 0.03, minsize = default.minsize)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
 \item{minvi}{minimum size of a violation that is reported}
 \item{minsize}{minimum size of a rest score group. By default
   \code{minsize} = \eqn{N/10} if \eqn{N \geq 500};
   \code{minsize} = \eqn{N/5 if} \eqn{250 \leq N < 500}; and
   \code{minsize} = \eqn{\max(N/3,50)} if \eqn{N < 250}}
}

\details{
 The output is of class monotonicity.class, and is often numerous. 
 Functions \code{plot} and \code{summary} can be used to summarize the output. 
 See Van der Ark (2007) for an example.      
}

\value{
  \item{results}{A list with as many components as there are items. Each component itself is also
                 a list containing the results of the check of monotonicity.}
  \item{I.labels}{The item labels}
  \item{Hi}{The item scalability coefficients Hi}
  \item{m}{The number of answer categories.}
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
\code{\link{coefH}}, \code{\link{check.iio}}, \code{\link{check.restscore}}, \code{\link{check.pmatrix}}, \code{\link{check.reliability}}, 
 \code{\link{plot.monotonicity.class}}, \code{\link{summary.monotonicity.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
monotonicity.list <- check.monotonicity(Communality)
plot(monotonicity.list)
summary(monotonicity.list)
}
\keyword{univar}
