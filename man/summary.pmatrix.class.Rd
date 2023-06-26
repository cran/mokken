\name{summary.pmatrix.class}
\alias{summary.pmatrix.class}

\title{
Summarize pmatrix.class objects
}
\description{
S3 Method for \code{summary} of objects of class pmatrix.class.
Summarize checks of nonintersection using method pmatrix.
}
\usage{
\method{summary}{pmatrix.class}(object, ...)
}

\arguments{
  \item{object}{list produced by \code{\link{check.pmatrix}}}
  \item{...}{Optional parameters will be ignored}  
}

\value{
 Matrix with \code{ncol(X)} rows and 10 columns, 
 showing for each item a summary of the violations of nonintersection using method pmatrix: 
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

\details{
For obtaining the P++ and P-- matrix, see examples.
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
    \doi{10.18637/jss.v020.i11}
 
}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}

\seealso{
\code{\link{check.pmatrix}}, \code{\link{plot.pmatrix.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
pmatrix.list <- check.pmatrix(Communality)
plot(pmatrix.list)
summary(pmatrix.list)

# Small example showing how to retrieve the P++ matrix and the P-- matrix
SmallExample <- acl[,1:4]
pmatrix.list <- check.pmatrix(SmallExample)
pmatrix.list$results$Ppp
pmatrix.list$results$Pmm
}

\concept{pmatrix}
