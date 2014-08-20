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
Returns a vector with as many elements as there are items, indicating to which scale an item belongs.
}
\usage{
aisp(X, search="normal", lowerbound=.3, alpha=.05, popsize=20, 
maxgens=default.maxgens, pxover=0.5, pmutation=0.1, verbose=FALSE)
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
  \item{lowerbound}{numeric scaling criterium; 0 <= \code{lowerbound} < 1.  The default is \code{0.3}.}
  \item{alpha}{Type I error level. The default is \code{0.05}.}
  \item{popsize}{Size of the population of items in genetic. algorithm The default is \code{20}.}
  \item{maxgens}{Number of generations in genetic algorithm. The default is \code{10^(log2(ncol(X)/5)) * 1000}.}
  \item{pxover}{Cross-over probability in genetic algorithm. The default is \code{0.5}.}
  \item{pmutation}{Mutation probability in genetic algorithm. The default is \code{0.1}.}
  \item{verbose}{Logical, indicating whether should output to the screen the results of the model. If \code{FALSE}, no output is produced. The default is \code{TRUE}.}
}
\details{
Each scale must consist of at least two items, hence the number of Mokken scales cannot exceed \code{ncol(X)/2}.
Procedure may be slow for large data sets. Especially if the genetic algorithm is used.
There is not yet an option \code{search="extended"}.
\code{aisp} replaces the function \code{search.normal} in earlier versions.
}

\value{An indicator vector of length J. 
       Each entry refers to an item. 
       Items with same integer belong to the same Mokken scale.
       A zero indicates an unscalable item.
       If \emph{n} is the largest integer, then \emph{n} Mokken scales were found.
}

\references{
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
}
\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl},
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

# Partition these 10 items into mokken scales using a genetic algorithm.
scale <- aisp(Communality,"ga",maxgens=1000)
coefH(Communality[,scale==1])

}
\keyword{aisp}
\keyword{scalability coefficients}
