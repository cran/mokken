\name{MLcoefH}
\alias{MLcoefH}

\title{
    Two-level scalability coefficents H
}

\description{
    Computes all types of two-level scalability coefficients (Snijders, 2001; Crisan et al., 2016),
    that is, between- and within-rater item-pair scalability coefficents and their ratio (HBij, HWij, and BWij, respectively),
    between- and within-rater item scalability coefficients and their ratio (HBi, HWi, and BWi, respectively),
    and between- and within-rater total scale scalability coefficients and their ratio (HB, HW, and BW, respectively).
    In addition, standard errors are estimated (Koopman, et al., 2020a) and if requested 
    incorporated in confidence intervals (Koopman, Zijlstra, & Van der Ark, 2020b, 2020c). Note that this version is an 
    adaptation of the estimation methods described in Snijders (1999) and in Koopman et al. (2020a), as the group proportions are now by 
    default weighted for group size (Koopman, Zijlstra, & Van der Ark, 2020b). As a result, the estimates for the 
    autonomySupport data differs from the printed table in Koopman et al. (2020a). The estimated variance-covariance matrix 
    of the coefficients can also be printed if requested.
}

\usage{
    MLcoefH(X, se = TRUE, ci = FALSE, nice.output = TRUE, subject = 1, 
            fixed.itemstep.order = NULL, weigh.props = TRUE, 
            type.ci = "WB", cov.mat = FALSE)
}

\arguments{
  \item{X}{matrix or data frame of numeric data containing a subject indicator column and the
           responses of \code{nrow(X)} raters to \code{ncol(X) - 1} items. 
           Missing values are not allowed}
  \item{se}{Logical: If \code{TRUE}, the standard errors are printed alongside the scalability coefficients}
  \item{ci}{The confidence level between 0 and 1 of the range-preserving confidence intervals. 
            If FALSE (default), no confidence intervals are printed (Koopman, Zijlstra, & Van der Ark, 2020b).}
  \item{nice.output}{Logical: If \code{TRUE}, scalability coefficients and standard errors are combined 
           in an a single object of class \code{noquote}. Item-pair ratios BWij are only printed if \code{FALSE}}
  \item{subject}{Represents the subject column. Default is column 1.}
  \item{fixed.itemstep.order}{matrix with number of rows equal to the number of item steps (m) and number of columns equal to the number of items (J). 
           The matrix should consis the integers 1 : (m * J), indicating a predefined order of the items steps with respect to popularity. 
           Value 1 indicates the easiest (most popular) item step, value (m * J) indicates the most difficult item step.}
  \item{weigh.props}{If \code{TRUE}: Use weighted proportions across groups to estimate coefficients and standard errors, if \code{FALSE}: 
           Use averaged proportions across groups to estimate coefficients and standard errors.}
  \item{type.ci}{If \code{WB}, Wald-based confidence interval are printed, if \code{RP} range-preserving confidence intervals are printed 
          (Koopman, Zijlstra, & Van der Ark, 2020b, 2020c). Default is \code{WB}. Only used if \code{ci} is specified.}
  \item{cov.mat}{Logical: If \code{TRUE}, the variance-covariance matrices of the estimated coefficients are printed. Default is \code{FALSE}.}
}

\details{
 If \code{se = TRUE} and \code{nice.output = TRUE}, the result is a list of 3 objects of class \code{noquote};  
 if \code{se = TRUE} and \code{nice.output = FALSE}, the result is a list of 3 matrices, one per set of coefficients; and 
 if \code{se = FALSE}, the result is a list of 3 matrices containing only the scalability coefficients;
 if \code{ci} is specified and \code{nice.output = TRUE}, there is one additional matrix for the ci's of the Hij coefficients.
 Computation of standard errors can be slow for a combination of many subjects and a large number of items;
 if \code{cov.mat = TRUE} three additional matrices are printed with the variance-covariances of Hij, Hi, and H.
}

\value{
 \item{Hij}{Scalability coefficients of the item pairs, upper triangle are the between-rater coefficients, and the
            lower triangle the within-rater coefficients (possibly with standard errors and/or confidence intervals; see details). 
            If \code{nice.output = FALSE} this returns a matrix with the type of coefficients per column}  
 \item{Hi}{Between- and within-rater scalability coefficients of the items and their ratio (possibly with standard errors and/or confidence intervals; see details)}  
 \item{H}{between- and within-rater scalability coefficients of the entire scale and their ratio (possibly with standard error and/or confidence intervals; see details)}  
}

\references{
Crisan, D. R., Van de Pol, J. E., & Van der Ark, L. A. (2016). 
Scalability coefficients for two-level polytomous item scores: An introduction and an application. 
In L. A. Van der Ark, D. M. Bolt, W.-C. Wang, & M. Wiberg (Eds.), 
\emph{Quantitative psychology research}. 
Springer. 
\url{https://doi.org/10.1007/978-3-319-38759-8_11}
  
Koopman, L., Zijlstra, B. J. H. & Van der Ark, L. A., (2020a). 
Standard errors of two-lvel scalability coefficients. 
\emph{British Journal of Statistical and Mathematical Psychology, 73,} 213-236. 
\url{https://doi.org/10.1111/bmsp.12174} 

Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020b). 
\emph{A two-step, test-guided Mokken scale analysis for nonclustered and clustered data}. 
Manuscript submitted for publication.

Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020c). 
\emph{Range-preserving confidence intervals for scalability coefficients in Mokken scale analysis}. 
Manuscript submitted for publication.

Snijders, T. A. B. (2001). Two-level non-parametric scaling for dichotomous data. 
In A. Boomsma, M. A. J. van Duijn, & T. A. B. Snijders (Eds.), 
\emph{Essays on item response theory} (pp. 319-338). 
Springer. 
\url{https://doi.org/10.1007/978-1-4613-0169-1_17}
}

\author{
   Letty Koopman \email{V.E.C.Koopman@uva.nl}
   L. Andries van der Ark \email{L.A.vanderArk@uva.nl} 
}

\seealso{
\code{\link{MLcoefZ}},\code{\link{coefH}}, 
}


\examples{

# Small data example (Koopman et al., 2019)
smallData <- data.frame(Subs = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3),
                        Xa   = c(0, 0, 1, 0, 1, 1, 1, 2, 1, 0, 1, 2, 0, 0, 0), 
                        Xb   = c(0, 0, 1, 0, 2, 2, 2, 1, 2, 1, 2, 2, 1, 1, 0), 
                        Xc   = c(1, 0, 0, 0, 1, 1, 2, 1, 2, 0, 1, 1, 2, 1, 0))

MLcoefH(smallData)

# Compute also the range-preserving confidence intervals
MLcoefH(smallData, ci = .95)

# Print variance-covariance matrices
MLcoefH(smallData, cov.mat = TRUE)

# Load real data example. Note that due to an estimation adaptation (Koopman et al., 2020)
# the results differ from the table in Koopman et al. (2019).
data(autonomySupport)

# Compute scalability coefficients with or without standard errors, range-preserving 
# confidence intervals, nice output
\donttest{          
H.se.nice <- MLcoefH(autonomySupport)
H.se.nice 
}
\donttest{          
H.se.not  <- MLcoefH(autonomySupport, nice.output = FALSE)
H.se.not  
}
\donttest{          
H.se.ci.nice <- MLcoefH(autonomySupport, ci = .95)
H.se.ci.nice
}
\donttest{          
H.se.ci.not <- MLcoefH(autonomySupport, ci = .95, nice.output = FALSE)
H.se.ci.not
}
\donttest{          
H.nice    <- MLcoefH(autonomySupport, se = FALSE)
H.nice    
}
\donttest{          
H.not     <- MLcoefH(autonomySupport, se = FALSE, nice.output = FALSE)
H.not     
}

}
\concept{two-level scalability coefficient}
\concept{scalability coefficient}
