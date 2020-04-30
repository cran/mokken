\name{aisp}
\alias{aisp}
\alias{search.normal}
\alias{search.extended}
\alias{search.ga}
\alias{search}


\title{
    Automated Item Selection Procedure (AISP) for Mokken Scale Analysis
}

\description{
    Returns a matrix with as many rows as there are items, 
    indicating to which scale an item belongs for each lowerbound.
}

\usage{
    aisp(X, lowerbound=.3, search="normal", alpha=.05, StartSet=FALSE, popsize=20, 
         maxgens=default.maxgens, pxover=0.5, pmutation=0.1, verbose=FALSE, 
         type.se = "delta", test.Hi = FALSE, level.two.var = NULL)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
  \item{search}{Type of item selection procedure: 
            "normal": Mokken's automated item selection procedure 
           (Mokken, 1971; Molenaar and Sijtsma, 2000; Sijtsma and Molenaar, 2002);
           "ga": item selection using a genetic algorithm
           (Straat, van der Ark, and Sijtsma, 2013). The default is "normal".}
  \item{lowerbound}{Value or vector with numeric scaling criteria; 0 <= \code{lowerbound} < 1.  The default is \code{0.3}.}
  \item{alpha}{Type I error level. The default is \code{0.05}.}
  \item{StartSet}{Startset of items for the first scale. Vector of item numbers. If \code{StartSet == FALSE} no startset is provided (default).}
  \item{popsize}{Size of the population of items in genetic. algorithm The default is \code{20}.}
  \item{maxgens}{Number of generations in genetic algorithm. The default is \code{10^(log2(ncol(X)/5)) * 1000}.}
  \item{pxover}{Cross-over probability in genetic algorithm. The default is \code{0.5}.}
  \item{pmutation}{Mutation probability in genetic algorithm. The default is \code{0.1}.}
  \item{verbose}{Logical, indicating whether should output to the screen the results of the model. If \code{FALSE}, no output is produced. The default is \code{TRUE}.}
  \item{type.se}{Indicates which type of standard errors is used in a Z-test whether coefficients meet the scaling criteria:
             "delta": uses standard errors approximated by the delta method (Kuijpers, Van der Ark, Kroon, 2013; Koopman, Zijlstra, Van der Ark, 2020);
             "Z": uses original Z-test (Mokken, 1971; Molenaar and Sijtsma, 2000; Sijtsma and Molenaar, 2002). The default is "delta".}
  \item{test.Hi}{If \code{FALSE}: tests if Hi is significantly larger than zero; 
              If \code{TRUE} tests if Hi is significantly larger than \code{lowerbound}. The default is \code{FALSE}.}
  \item{level.two.var}{vector of length \code{nrow(X)} or matrix with number of rows equal to \code{nrow(X)} 
          that indicates the level two variable for nested data (Koopman et al., 2020).}
}


\details{
    Each scale must consist of at least two items, hence the number of Mokken scales cannot exceed \code{ncol(X)/2}.
    Procedure may be slow for large data sets. Especially if the genetic algorithm is used.
    There is not yet an option \code{search="extended"}.
    \code{aisp} replaces the function \code{search.normal} in earlier versions.
}

\value{An matrix with J rows. 
       Each entry refers to an item. 
       Items with same integer belong to the same Mokken scale.
       A zero indicates an unscalable item.
       If \emph{n} is the largest integer, then \emph{n} Mokken scales were found.
}

\references{
  Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020). 
  \emph{A two-step procedure for scaling multilevel data using Mokken's scalability coefficients}. 
  Manuscript submitted for publication.

  Kuijpers, R. E., Van der Ark, L. A., and Croon, M. A. (2013).
  Standard errors and confidence intervals for scalability coefficients in Mokken scale analysis using marginal models.
  \emph{Sociological Methodology, 43,} 42-69.
  
  Mokken, R. J. (1971) 
  \emph{A Theory and Procedure of Scale Analysis}. 
  Berlin, Germany: De Gruyter.
  
  Molenaar, I.W. and Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  Groningen, The Netherlands: IEC ProGAMMA.
  
  Sijtsma, K, and Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Thousand Oaks, CA: Sage.

   Straat, J. H., Van der Ark, L. A., and Sijtsma, K. (2013).
  Comparing optimization algorithms for item selection in Mokken scale analysis.
  \emph{Journal of Classification, 30,} 72-99.
  
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{http://www.jstatsoft.org/v20/i11}

  Van der Ark, L. A. (2012). 
  New developments in {M}okken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}, 48. 
  doi:10.18637/jss.v048.i05
}
\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl},
 D. van den Bergh
 L. Koopman
 J. H. Straat
}

\seealso{
\code{\link{coefH}}, \code{\link{check.iio}}, \code{\link{check.monotonicity}}, \code{\link{check.pmatrix}}, \code{\link{check.reliability}},\code{\link{check.restscore}} 
}

\examples{
data(acl)

# Select the scale Communality consisting of 10 items.
Communality <- acl[,1:10]

# Partition these 10 items into mokken scales using Mokken's automated item selection procedure.
scale <- aisp(Communality)
coefH(Communality[,scale==1])

# Same but using items 1 and 2 in the startset.
scale <- aisp(Communality, StartSet = c(1, 2), verbose = TRUE)
coefH(Communality[,scale==1])

# Perform aisp for increasing lowerbounds
\donttest{
scales <- aisp(Communality, lowerbound = seq(0, .55, .05))
scales
}

# Use a significant test for criteria Hi > c (rather than the point estimate)
scale <- aisp(Communality, test.Hi = TRUE, verbose = TRUE)
coefH(Communality[,scale==1])

# Partition these 10 items into mokken scales using a genetic algorithm.
scale <- aisp(Communality,search="ga",maxgens=1000) 
coefH(Communality[,scale==1])

# Perform aisp on two-level data
data(autonomySupport)
scores <- autonomySupport[, -1]
classes <- autonomySupport[, 1]
scale <- aisp(scores, level.two.var = classes)
coefH(scores[, scale==1], level.two.var = classes)

}

\concept{aisp}
\concept{automated item selection procedure}
\concept{scalability coefficients}
