\name{summary.iio.class}
\alias{summary.iio.class}

\title{
Summarize iio.class objects
}
\description{
S3 Method for \code{summary} of objects of class iio.class.
Summarize checks of invariant item ordering.
}
\usage{
\method{summary}{iio.class}(object, ...)
}

\arguments{
  \item{object}{list produced by \code{\link{check.iio}}}
   \item{...}{Optional parameters will be ignored}  

}

\value{
 \item{method}{String describing the method used for investigating invariant item ordering: 
   Either "MIIO" (Method Manifest Invariant Item Ordering), 
   "MSCPM" (Method Manifest Scale Cumulative Probability Model), or "IT" (Method Increasingness in Transposition)}
 \item{item.summary}{Matrix with \code{ncol(X)} rows and 9 columns, 
   showing for each item a summary of the violations of an invariant item ordering.}
 \item{backward.selection}{Matrix showing the number of violations for each item (rows) at each step
 of the backward item selection proces (columns).}
 \item{HT}{Numeric: Coefficient HT for the selected items. Given an IIO, coefficient HT expresses the strength of the ordering.}  
}

\references{
  Ligtvoet, R., van der Ark, L. A., The Marvelde, J. M., and Sijtsma (in press) 
  Investigating an invariant item ordering for polytomously scored items. 
  \emph{Educational and Psychological Measurement}.

  Ligtvoet, R., van der Ark, L. A., Bergsma, W. P., and Sijtsma (2009) 
  Polytomous latent scales for the investigation of the ordering of items.
  Manuscript submitted for publication. 


  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{http://www.jstatsoft.org}
}

\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}


\seealso{
\code{\link{check.iio}}
}

\examples{
# Examples from Ligtvoet et al. (2009).
 
data(cavalini)
X1 <- cavalini[,c(3,5,6,7,9,11,13,14)]

# Use Method MIIO and remove items violating MIIO
iio.list1 <- check.iio(X1)
summary(iio.list1)
X2 <- X1[,is.na(charmatch(dimnames(X1)[[2]],names(iio.list1$items.removed)))]

# Use Method MSCPM and remove items violating MSCPM
iio.list2 <- check.iio(X2,method="MSCPM")
summary(iio.list2)
X3 <- X2[,is.na(charmatch(dimnames(X2)[[2]],names(iio.list2$items.removed)))]

# Use Method IT
iio.list3 <- check.iio(X3,method="IT")
summary(iio.list3)
}

\keyword{print}
