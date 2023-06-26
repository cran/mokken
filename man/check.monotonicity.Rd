\name{check.monotonicity}
\alias{check.monotonicity}

\title{
Check of Monotonicity
}

\description{
Returns a list (of class \code{monotonicity.class}) with results from the investigation of monotonicity (Junker & Sijtsma, 2000; Mokken, 1971; Molenaar & Sijtsma, 2000; Sijtsma & Molenaar, 2002).

For two-level test data (clustered respondents) argument level.two.var exist, such that two lists are returned, containing the results for level 1 (person level) and level 2 (cluster level), respectively. Only method MIIO is implemented for two-level test data. 
}

\usage{
check.monotonicity(X, minvi = 0.03, minsize = default.minsize, level.two.var = NULL)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
 \item{minvi}{minimum size of a violation that is reported}
 \item{minsize}{minimum size of a rest score group. By default
   \code{minsize} = \eqn{N/10} if \eqn{N \ge 500};
   \code{minsize} = \eqn{N/5 if} \eqn{250 \le N < 500}; and
   \code{minsize} = max\eqn{(N/3,50)} if \eqn{N < 250}}
 \item{level.two.var}{Add respondent-clustering variable to get results for Level 1 (person level) and Level 2 (cluster level; see Koopman et al., 2023a,b)}.
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
  Junker, B.W., & Sijtsma, K. (2000). 
  Latent and manifest monotonicity in item response models. 
  \emph{Applied Psychological Measurement, 24}, 65-81.
  \doi{10.1177/01466216000241004}

  Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2023a). 
  Assumptions and Properties of Two-Level Nonparametric Item Response 
  Theory Models. Manuscript submitted for publication.
  
  Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2023b). 
  Evaluating Model Fit in Two-Level Mokken Scale Analysis. 
  Manuscript submitted for publication.
  
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
\code{\link{check.errors}}, 
\code{\link{check.iio}}, 
\code{\link{check.restscore}}, 
\code{\link{check.pmatrix}}, 
\code{\link{check.reliability}}, 
\code{\link{coefH}}, 
\code{\link{plot.monotonicity.class}}, 
\code{\link{summary.monotonicity.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
monotonicity.list <- check.monotonicity(Communality)
plot(monotonicity.list)
summary(monotonicity.list)

# Compute two-level fit statistics (Koopman et al., 2023a, 2023b)
data("autonomySupport")
dat <- autonomySupport[, -1]
groups <- autonomySupport[, 1]
autonomyMM <- check.monotonicity(dat, level.two.var = groups)
summary(autonomyMM)
plot(autonomyMM)

}
\concept{monotonicity}
\concept{manifest monotonicity}
