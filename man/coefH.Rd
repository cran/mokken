\name{coefH}
\alias{coefH}

\title{
Scalability coefficents H
}

\description{
Computes item-pair scalability coefficents Hij, item scalability coefficents Hi, and scale scalability coefficent H 
(Loevinger, 1948; Mokken, 1971, pp. 148-153; Molenaar & Sijtsma, 2000, pp. 11-13; Sijtsma & Molenaar, chap. 4; Van der Ark, 2007; 2010),  
as well as their standard errors (Kuijpers, Van der Ark, & Croon, 2013; also see Van der Ark, Croon, & Sijtsma, 2008)
and possibly confidence intervals (Koopman, Zijlstra, & Van der Ark, 2020a, 2020b). 
Mokken's coefficients and standard errors can also be estimated in two-level data (Koopman, Zijlstra, & Van der Ark, 2020a).
It is also possible to compare scalability coefficients across groups using the item-step ordering of the entire sample 
(cf. \code{CHECK=GROUPS} option in MSP; Molenaar and Sijtsma, 2000). The estimated variance-covariance matrix of the coefficients 
is invisible but can be printed by saving the result, see examples.
 }

\usage{
coefH(X, se = TRUE, ci = FALSE, nice.output = TRUE, level.two.var = NULL, 
     group.var = NULL, fixed.itemstep.order = NULL, type.ci = "WB", 
     results = TRUE)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
  \item{se}{Logical: If \code{TRUE}, the standard errors of the scalability coefficients are given}
  \item{ci}{The confidence level between 0 and 1 of the range-preserving confidence intervals. 
            If \code{FALSE} (default), no confidence intervals are printed (Koopman, Zijlstra, & Van der Ark, 2020b).}
  \item{nice.output}{Logical: If \code{TRUE}, scalability coefficients and standard errors are combined in an a single object of class \code{noquote}}
  \item{level.two.var}{vector of length \code{nrow(X)} or matrix with number of rows equal to \code{nrow(X)} 
          that indicates the level two variable for nested data to get appropriate standard errors (Koopman et al., 2020a.}
  \item{group.var}{vector of length \code{nrow(X)} or matrix with number of rows equal to \code{nrow(X)} to be used as grouping variable}
  \item{fixed.itemstep.order}{matrix with number of rows equal to the number of item steps (m) and number of columns equal to the number of items (J). 
        The matrix should consis the integers 1 : (m * J), indicating a predefined order of the items steps with respect to popularity. 
        Value 1 indicates the easiest (most popular) item step, value (m * J) indicates the most difficult item step.}
  \item{type.ci}{If \code{WB}, Wald-based confidence interval are printed, if \code{RP} range-preserving confidence intervals are printed (Koopman, Zijlstra, & Van der Ark, 2020b, 2020c). Default is \code{WB}. Used only if \code{ci} has been specified.}
  \item{results}{Logical: If TRUE results are printed to the screen. Option \code{FALSE} is useful only for some internal functions}      
}

\details{
 May not work if any of the item variances equals zero. 
 Such items should not be used in a test and removed from the data frame.
 
 If \code{nice.output = TRUE} and \code{se = TRUE}, the result is a list of 3 objects of class \code{noquote};  
 if \code{nice.output = FALSE} and \code{se = TRUE}, the result is a list of 6 matrices (3 for the scalability coefficients and 3 for the standard errors); and 
 if \code{se = FALSE}, the result is a list of 3 matrices (for the scalability coefficients);
 if \code{ci} is specified and \code{se = TRUE} or \code{nice.output = FALSE}, there is one additional matrix for the ci's of the Hij coefficients;
 if \code{level.two.var} is not null the standard errors are adjusted to take the nesting into account;
 if \code{group.var = Y} with Y having K values, an additional element named \code{Groups} is added to the list.
 Element \code{Groups} shows the scalability coefficients per group ordered by means of \code{\link{sort}} (see \code{\link{Sys.getlocale}} for details).
 \code{group.var} returns coefficients for groups containing at least two case.
 Computation of standard errors can be slow for a combination of a large sample size and a large number of items.
}

\value{
 \item{Hij}{scalability coefficients of the item pairs (possibly with standard errors; see details)}  
 \item{Hi}{vector containing scalability coefficients of the items (possibly with standard errors; see details)}  
 \item{H}{scalability coefficient of the entire scale (possibly with standard error; see details)}  
 \item{se.Hij}{standard errors of the scalability coefficients of the item pairs (only if \code{nice.output = FALSE} and \code{se = TRUE}; see details)}  
 \item{se.Hi}{standard errors of the scalability coefficients of the items (see details)}  
 \item{se.H}{standard error of the scalability coefficient of the entire scale (see details)} 
 \item{ci.Hij}{confidence intervals of the scalability coefficients of the item pairs (only if \code{nice.output = FALSE} and/or \code{se = TRUE}; see details)}  
 \item{ci.Hi}{confidence intervals of the scalability coefficients of the items (see details)}  
 \item{ci.H}{confidence intervals of the scalability coefficient of the entire scale (see details)} 
 \item{Groups}{Scalability coefficients for subgroups (see details)} }

\references{
  Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020a). 
  \emph{A two-step, test-guided Mokken scale analysis for nonclustered and clustered data}. 
  Manuscript submitted for publication.

  Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020b). 
  \emph{Range-preserving confidence intervals for scalability coefficients in Mokken scale analysis}. 
  Manuscript in preparation.

  Kuijpers, R. E., Van der Ark, L. A., & Croon, M. A. (2013).
  Standard errors and confidence intervals for scalability coefficients in Mokken scale analysis using marginal models.
  \emph{Sociological Methodology, 43,} 42-69.
  \url{https://doi.org/10.1177/0081175013481958}
  
  Loevinger, J. (1948). 
  The technique of homogeneous tests compared with some aspects of 'scale analysis' and factor analysis.
  \emph{Psychological Bulletin, 45}, 507-530.
     
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
  \emph{Journal of Statistical Software, 20} (11), 1-19.
  \url{https://www.jstatsoft.org/article/view/v020i11}

  Van der Ark, L. A. (2010). 
  Getting started with Mokken scale analysis in \pkg{R}. 
  Unpublished manuscript.
  \url{https://sites.google.com/a/tilburguniversity.edu/avdrark/mokken}

  Van der Ark, L. A., Croon, M. A., & Sijtsma (2008). 
  Mokken scale analysis for dichotomous items using marginal models.
  \emph{Psychometrika, 73}, 183-208.
  \url{https://doi.org/10.1007/s11336-007-9034-z}
 }

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl} 
 L. Koopman}

\seealso{
\code{\link{coefZ}}, \code{\link{search.normal}} }

\examples{
data(acl)
Communality <- acl[, 1:10]

# Compute scalability coefficients and standard errors
Hs <- coefH(Communality)

# Compute scalability coefficients, standard errors, and range-preserving confidence intervals
coefH(Communality, ci = .95)

# Scalability coefficients but no standard errors
coefH(Communality, se = FALSE)

# Scalability coefficients for different groups:
subgroup <- ifelse(acl[,11] < 2,1,2)
coefH(Communality, group.var = subgroup)

# Extract variance-covariance matrices
attributes(Hs)
Hs$covHij
Hs$covHi
Hs$covH

# Nested data:
data(autonomySupport)
scores <- autonomySupport[, -1]
classes <- autonomySupport[, 1]
coefH(scores, level.two.var = classes, ci = .95)

}
\concept{scalability coefficient}
