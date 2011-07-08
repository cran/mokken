\name{coefZ}
\alias{coefZ}

\title{
Computation of Z-Values 
}

\description{
Computes Zij-values of item pairs, Zi-values of items, and Z-value of the entire scale,
which are used to test whether Hij, Hi, and H, respectively, are significantly 
greater than zero (Molenaar and Sijtsma, 2000, pp. 59-62; Sijtsma and Molenaar, p. 40; Van der Ark, 2007; 2010).

Used in the function \code{aisp}
}
\usage{
coefZ(X)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
}

\value{
 \item{Zij}{matrix containing the z-values of the item-pairs}
 \item{Zi}{vector containing z-values of the items}
 \item{Z}{real z-value of the entire scale}
}


\references{
 
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

  Van der Ark, L. A. (2010). 
  Getting started with Mokken scale analysis in \pkg{R}. 
  Unpublished manuscript.
  \url{http://cran.r-project.org/web/packages/mokken/vignettes/mokken.pdf}
 }

\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}


\seealso{
\code{\link{coefH}}, \code{\link{aisp}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
coefH(Communality)
coefZ(Communality)
}
\keyword{scalability coefficients}
