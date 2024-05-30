\name{check.iio}
\alias{check.iio}

\title{
Check of Invariant Item Ordering
}

\description{
Returns a list (of class \code{iio.class}) with results from the investigation of invariant item ordering.
Three methods may be used for the investigation of invariant item ordering. 
(1) Method MIIO (manifest invariant item ordering: investigates the manifest item response functions for all pairs of items). For polytomous items, t-tests are used to test violations, for dichotomous items z-tests are used to test are used to test violations.
(2) Method MS-CPM (manifest scale - cumulative probability model: investigates the manifest item step response functions for all pairs of items). Z-tests are used to test violations.
(3) Method IT (increasingness in transposition: investigates all bivariate joint probabilities for all pairs of items). Chi-square tests are used to test violations.

For a complete description of Method MIIO, see Ligtvoet, Van der Ark, Te Marvelde, and Sijtsma (2010); for
a complete description of the Method MS-CPM and Method IT with reference to Method MIIO, see Ligtvoet,
Van der Ark, Bergsma, and Sijtsma (2011). For a discription of investigating the ordering structure of clustered items,
see Koopman & Braeken (2024).

For two-level test data (clustered respondents) argument level.two.var exist, for clustered item data, argument fixed.item.order exist. For both arguments, two lists are returned, containing the results for level 1 (person or item level) and level 2 (group or cluster level), respectively. Only method MIIO is implemented for two-level and clustered-item test data. 
}
\usage{
check.iio(X, method="MIIO", minvi = default.minvi, minsize = default.minsize, 
alpha = .05, item.selection=TRUE, verbose=FALSE, fixed.item.order = NULL,
level.two.var = NULL)
}

\details{
 The output is of class \code{iio.class}, and is often numerous. 
 Functions \code{plot} and \code{summary} can be used to summarize the output. 
 See Van der Ark (2014) for an example. For an example of clustered items,
see Koopman & Braeken (2024).
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
           containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items.  
           Missing values are not allowed}
  \item{method}{Either \code{"MIIO"} (default), \code{"MSCPM"}, or \code{"IT"}. Partial matching is allowed (e.g. \code{method="ms"} is equivalent to \code{method="MSCPM"})}         
  \item{minvi}{minimum size of a violation that is reported. By default 
    \code{minvi} = .03 times the number of item step response functions (m) for Method MIIO;
    \code{minvi} = .03 for Method IT and Method MSCPM}
  \item{minsize}{minimum size of a rest score group.  By default
   \code{minsize} = \eqn{N/10} if \eqn{N \ge 500};
   \code{minsize} = \eqn{N/5 if} \eqn{250 \le N < 500}; and
   \code{minsize} = max\eqn{(N/3,50)} if \eqn{N < 250}}
  \item{alpha}{Nominal Type I error for t test (Method MIIO), z test (Method MSCPM), or McNemar test (Method IT). 
  Default \code{alpha = .05}} 
  \item{item.selection}{Conduct backward item selection procedure (see Ligtvoet et al., 2010).
  Default \code{item.selection=TRUE}}
  \item{verbose}{Show the results of the backward item selection algorithm on screen.
  Default \code{verbose=FALSE}}
  \item{fixed.item.order}{Matrix or vector containing J numeric values to indicate the item ordering from easy to difficult, to perform a confirmatory analysis of manifest invariant item/cluster ordering. For clustered items, the cluster numbers are given, which are repeated for each item of that cluster, such that the length is still J, for example c(1, 1, 2, 2, 3, 3) for three clusters of two items, of which the first cluster is easiest and the last cluster is most difficult (see Koopman & Braeken, 2024).}
  \item{level.two.var}{Add respondent-clustering variable to get results for Level 1 (person level) and Level 2 (cluster level; see Koopman et al., 2023a,b)}.
}

\value{
  \item{results}{A list with as many components as there are item pairs.  
  Each component itself is also
  a list containing the results of the investigation of IIO.}
  \item{violations}{A matrix: Summary of the backward item selection (Corresponds to Table 4 in Ligtvoet et al., 2010,
   and Table 1 in Ligtvoet et al., 2011). 
  The first column gives, for each item, the number of violations of IIO. 
  If the number of violations is nonzero, then the item with the largest number of violations is removed.
  If two or more items have the maximum number of violations, then from those items the item producing the lowest value of Loevinger's H is removed.
  The second column shows the number of violations with one item removed, the third column shows the number of violations with two items removed, etc.}
  \item{items.removed}{List of the items removed in chronological order}
  \item{HT}{Coefficient HT for the remaining items. For the use of coefficient HT see Ligtvoet et al. (2010). If the sample size is extremely large coefficient HT is estimated using a random subsample. For clustered items, coefficient HBT and ratio HBT/HT is also estimated, see Koopman & Braeken (2024).}
  \item{method}{The argument \code{method}}
  \item{item.mean}{The mean item scores}
}

\references{
  Koopman, L. & Braeken, J. (2024). 
  Investigating the Ordering Structure of Clustered Items Using Nonparametric 
  Item Response Theory. Manuscript submitted for publication.
  
  Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2023a). 
  Assumptions and Properties of Two-Level Nonparametric Item Response 
  Theory Models. Manuscript submitted for publication.
  
  Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2023b). 
  Evaluating Model Fit in Two-Level Mokken Scale Analysis. 
  Manuscript submitted for publication.
  
  Ligtvoet, R., L. A. van der Ark, J. M. te Marvelde, & K. Sijtsma (2010).
  Investigating an invariant item ordering for polytomously scored items. 
  \emph{Educational and Psychological Measurement, 70}, 578-595.
  \doi{10.1177/0013164409355697}

  Ligtvoet, R., L. A. van der Ark, W. P. Bergsma, & K. Sijtsma (2011).
  Polytomous latent scales for the investigation of the ordering of items.
  \emph{Psychometrika, 76}, 200-216.
  \doi{10.1007/s11336-010-9199-8}

  Sijtsma, K., R. R. Meijer, & Van der Ark, L. A. (2011).
  Mokken scale analysis as time goes by: An update for scaling practitioners.
  \emph{Personality and Individual Differences, 50}, 31-37.
  \doi{10.1016/j.paid.2010.08.016} 

  Sijtsma, K., & Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Sage.

  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software, 20} (11), 1-19.
    \doi{10.18637/jss.v020.i11}


 Van der Ark, L. A. (2012).
  New developments in Mokken scale analysis in R.
 \emph{Journal of Statistical Software, 48}(5), 1-27.
  \doi{10.18637/jss.v048.i05}

}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}


\seealso{
\code{\link{check.errors}}, 
\code{\link{check.monotonicity}}, 
\code{\link{check.pmatrix}}, 
\code{\link{check.reliability}} 
\code{\link{check.restscore}}, 
\code{\link{coefH}}, 
\code{\link{plot.iio.class}},
\code{\link{summary.iio.class}},
}

\examples{

# Examples from Ligtvoet et al. (2010).

data(acl)            

Order <- acl[,41:50]               
summary(check.iio(Order))
plot(check.iio(Order))

Autonomy <- acl[,91:100]
summary(check.iio(Autonomy))
plot(check.iio(Autonomy))

# Examples from Ligtvoet et al. (2011).
 
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

# Examples for investigating the ordering structure of a clustered item set 
# (Koopman & Braeken, 2024)
 \donttest{     
data("trog")
clusters <- rep(1:20, each = 4)
ico <- check.iio(trog, item.selection = FALSE, fixed.item.order = clusters)
summary(ico)
}

# Compute two-level fit statistics (Koopman et al., 2023a, 2023b)
 \donttest{     
data("autonomySupport")
dat <- autonomySupport[, -1]
groups <- autonomySupport[, 1]
autonomyMIIO <- check.iio(dat, item.selection = FALSE, level.two.var = groups)
summary(autonomyMIIO)
plot(autonomyMIIO)
}

}
\concept{iio}
\concept{miio}
\concept{invariant item ordering}
\concept{manifest invariant item ordering}

