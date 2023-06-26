\name{coefZ}
\alias{coefZ}

\title{
Computation of Z-Values 
}

\description{
Computes Zij-values of item pairs, Zi-values of items, and Z-value of the entire scale,
which are used to test whether Hij, Hi, and H, respectively, are significantly 
greater than zero using the original method Z
(Molenaar and Sijtsma, 2000, pp. 59-62; Sijtsma and Molenaar, p. 40; Van der Ark, 2007; 2010)
or the Wald-based method (\code{WB}) or range-preserving method (\code{RP}) 
(Kuijpers et al., 2013; Koopman et al., in press a, in press b).
The Wald-based method  and range-preserving method  can also handle nested data and can test other lowerbounds than zero.
Used in the function \code{aisp}
}
\usage{
coefZ(X, lowerbound = 0, type.z = "Z", level.two.var = NULL)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
  \item{lowerbound}{Value of the null hypothesis to which the scalability are compared to compute the Z-score (see details),
          0 <= \code{lowerbound} < 1. The default is 0.}
  \item{type.z}{Indicates which type of z-score is computed:
          "WB": Wald-based z-score based on standard errors as approximated by the delta method 
                (Kuijpers et al., 2013; Koopman et al., in press a);
          "RP": Range-preserving z-score, also based on the delta method (Koopman et al., in press b);
          "Z": uses original Z-test and is only appropriate to test lowerbound = 0 
	  (Mokken, 1971; Molenaar and Sijtsma, 2000; Sijtsma and Molenaar, 2002). The default is "Z".}
  \item{level.two.var}{vector of length \code{nrow(X)} or matrix with number of rows equal to \code{nrow(X)} 
          that indicates the level two variable for nested data (Koopman et al., in press a).}
}

\value{
 \item{Zij}{matrix containing the Z-values of the item-pairs}
 \item{Zi}{vector containing Z-values of the items}
 \item{Z}{Z-value of the entire scale}
}

\details{
For the estimated item-pair coefficient \eqn{Hij} with standard error \eqn{SE(Hij)}, the Z-score  is computed as   
\deqn{Zij = (Hij - lowerbound) / SE(Hij)}  if \code{type.z = "WB"}, and the Z-score  is computed as
\deqn{Zij = -(log(1 - Hij) - log(1 - lowerbound)) / (SE(Hij) / (1 - Hij))} if \code{type.z = "RP"}
(Koopman et al., in press b).
For the estimate item-scalability coefficients \eqn{Hi} and total-scalbility coefficients \eqn{H} a similar procedure is used. 
Standard errors of the Z-scores are not provided.
}

\references{
   Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (in press a).
   A two-step, test-guided Mokken scale analysis for nonclustered and clustered data. 
   \emph{Quality of Life Research}. (advanced online publication)
   \doi{10.1007/s11136-021-02840-2}

   Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (in press b).
   Range-preserving confidence intervals and significance tests for scalability coefficients in Mokken scale analysis.
   In M. Wiberg, D. Molenaar, J. Gonzalez, & Kim, J.-S. (Eds.),
   \emph{Quantitative Psychology; The 1st Online Meeting of the Psychometric Society, 2020}.
   Springer.
   \doi{10.1007/978-3-030-74772-5_16}

 
   Kuijpers, R. E., Van der Ark, L. A., & Croon, M. A. (2013).
  Standard errors and confidence intervals for scalability coefficients in Mokken scale analysis using marginal models.
  \emph{Sociological Methodology, 43,} 42-69.
  \doi{10.1177/0081175013481958}
 
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


  Van der Ark, L. A. (2010). 
  Getting started with Mokken scale analysis in \pkg{R}. 
  Unpublished manuscript.
  \url{https://sites.google.com/a/tilburguniversity.edu/avdrark/mokken}
 }

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
 L. Koopman 
}


\seealso{
\code{\link{coefH}}, \code{\link{aisp}}
}

\examples{
data(acl)
Communality <- acl[,1:10]

# Compute the Z-score of each coefficient
coefH(Communality)
coefZ(Communality)

# Using lowerbound .3
coefZ(Communality, lowerbound = .3, type.z = "WB")

# Z-scores for nested data
data(autonomySupport)
scores <- autonomySupport[, -1]
classes <- autonomySupport[, 1]
coefH(scores, level.two.var = classes)
coefZ(scores, type.z = "WB", level.two.var = classes)


}
\concept{scalability coefficient}
