\name{summary.monotonicity.class}
\alias{summary.monotonicity.class}

\title{
Summarize monotonicity.class objects
}

\description{
S3 Method for \code{summary} of objects of class monotonicity.class.
Summarizes checks of monotonicity
}
\usage{
\method{summary}{monotonicity.class}(object, ...)
}

\arguments{
  \item{object}{list produced by \code{\link{check.monotonicity}}}
  \item{...}{Optional parameters will be ignored}  
}


\value{
 Matrix with \code{ncol(X)} rows and 10 columns, 
 showing for each item a summary of the violations of monotonicity: 
    \code{itemH} = Item-scalability coefficient; 
    \code{#ac} = number of active pairs that were investigated; 
    \code{#vi} = number of violations in which the item is involved;
    \code{#vi/#ac} = propotion of active pairs that is involved in a violation;
    \code{maxvi} = maximum violation;
    \code{sum} = sum of all violations;
    \code{zmax}  = maximum z-value;
    \code{zsig} = number of significant z-values;
    \code{crit} = Crit value (Molenaar & Sijtsma, 2000, pp. 49, 74).
}

\references{
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
}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}

\seealso{
\code{\link{check.monotonicity}}, \code{\link{plot.monotonicity.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
monotonicity.list <- check.monotonicity(Communality)
plot(monotonicity.list)
summary(monotonicity.list)
}

\concept{manifest monotonicity}
\concept{monotonicity}
