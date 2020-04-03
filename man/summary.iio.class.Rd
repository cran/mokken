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
 \item{item.summary}{Matrix with \code{ncol(X)} rows and 10 columns, 
   showing for each item a summary of the violations of an invariant item ordering:
    \code{itemH} = Item-scalability coefficient; 
    \code{#ac} = number of active pairs that were investigated; 
    \code{#vi} = number of violations in which the item is involved;
    \code{#vi/#ac} = propotion of active pairs that is involved in a violation;
    \code{maxvi} = maximum violation;
    \code{sum} = sum of all violations;
    \code{tmax} (for method MIIO), \code{zmax} (for method MSCPM), or \code{xmax} (for method IT) = maximum t-value, z-value, and chi-square value, respectively;
    \code{tsig} (for method MIIO), \code{zsig} (for method MSCPM), or \code{xsig} (for method IT) = number of significant t-values, z-values, and chi-square values, respectively;
    \code{crit} = Crit value (Molenaar and Sijtsma, 2000, pp. 49, 74).}
 \item{backward.selection}{Matrix showing the number of violations for each item (rows) at each step of the backward item selection proces (columns).
    The first column shows the number of violations for each item. Then in an iterative procedure, the item whose removal results in the largest 
    decrease of violations is removed and the number of violations is computated again. If the reduction in the number of violations is undecisive
    then, from the candidate items, the item having the smallest scalability coefficient is removed. The backward selection procedure stops when there
    are no more violations.
 }
 \item{HT}{Numeric: Coefficient HT for the selected items. Given an IIO, coefficient HT expresses the strength of the ordering (Ligtvoet et al., 2010).}  
}

\references{
  Ligtvoet, R., L. A. van der Ark, J. M. te Marvelde, and K. Sijtsma (2010).
  Investigating an invariant item ordering for polytomously scored items. 
  \emph{Educational and Psychological Measurement, 70}, 578-595.

  Ligtvoet, R., L. A. van der Ark, W. P. Bergsma, and K. Sijtsma (2011).
  Polytomous latent scales for the investigation of the ordering of items.
  \emph{Psychometrika, 76}, 200-216.

  Molenaar, I.W. and Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  Groningen, The Netherlands: IEC ProGAMMA.

  Sijtsma, K., R. R. Meijer, and L. A. van der Ark (2011).
  Mokken scale analysis as time goes by: An update for scaling practitioners.
  \emph{Personality and Individual Differences, 50}, 31-37.

  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software, 20} (11), 1-19.
  \url{http://www.jstatsoft.org/v20/i11}

  Van der Ark, L. A. (2012).
  New developments in Mokken scale analysis in R.
 \emph{Journal of Statistical Software, 48}(5), 1-27.
  \url{http://www.jstatsoft.org/v48/i5}

}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}


\seealso{
\code{\link{check.iio},\link{plot.iio.class}}
}

\examples{
# Examples from Ligtvoet et al. (2011).
 
data(cavalini)
X1 <- cavalini[,c(3,5,6,7,9,11,13,14)]

# Use Method MIIO and remove items violating MIIO
iio.list1 <- check.iio(X1)
summary(iio.list1)
plot(iio.list1)
X2 <- X1[,is.na(charmatch(dimnames(X1)[[2]],names(iio.list1$items.removed)))]

# Use Method MSCPM and remove items violating MSCPM
iio.list2 <- check.iio(X2,method="MSCPM")
summary(iio.list2)
X3 <- X2[,is.na(charmatch(dimnames(X2)[[2]],names(iio.list2$items.removed)))]

# Use Method IT
iio.list3 <- check.iio(X3,method="IT")
summary(iio.list3)
}

\concept{iio}
\concept{miio}
\concept{invariant item ordering}
\concept{manifest invariant item ordering}
