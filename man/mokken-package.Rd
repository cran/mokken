\name{mokken-package}
\alias{mokken-package}
\alias{mokken}
\docType{package}

\title{
Mokken Scale Analysis
}
\description{
Mokken scale analysis (Mokken, 1971; Sijtsma & Molenaar, 2002; Sijtsma & Van der Ark, 2017) is a scaling procedure 
for both dichotomous and polytomous items.
It consists of an item selection algorithm to partition a set of items into Mokken scales
 and several methods to check the assumptions of two nonparametric item response theory
models: the monotone homogeneity model and the double monotonicity model.
The output of this R-package resembles the output of the stand-alone program MSP (Molenaar & Sijtsma, 2000).
}
\details{
\tabular{ll}{
Package:    \tab mokken      \cr
Type:       \tab Package     \cr
Version:    \tab 3.0.6       \cr
Date:       \tab 2021-03-08  \cr
License:    \tab GPL Version 2 or later\cr
}
The package contains principal functions for Mokken scale analysis.

The package contains the following data sets
\tabular{ll}{
\code{\link{acl}}                \tab  Scores on a personality checklist. \cr
\code{\link{autonomySupport}}    \tab  Scores from students on their teacher's autonomy support \cr
\code{\link{balance}}            \tab  Scores on balance taskts \cr
\code{\link{cavalini}}           \tab  Scores on an inventory on industrial malodor \cr
\code{\link{DS14}}               \tab  Scores on a Type D test (bootstrap sample) \cr
\code{\link{mcmi}}               \tab  Scores on some items from the Dutch version of the Millon Clinical Multiaxial Inventory \cr
\code{\link{SWMD}}               \tab  Scores from pupils nested in classrooms on their well-being with teachers \cr
\code{\link{SWMDK}}              \tab  Scores from pupils nested in classrooms on their well-being with teachers and classmates \cr
\code{\link{transreas}}          \tab  Scores on a transitive reasoning test \cr
\code{\link{transreas2}}         \tab  More scores on a transitive reasoning test \cr
}

A guide for Mokken scale analysis in R for people who do not know R (Van der Ark, 2010)  is available as a vignette from \url{https://sites.google.com/a/tilburguniversity.edu/avdrark/mokken}.
The Mokken package was created by Andries van der Ark, who is also the maintainer. Significant parts have been developed by Letty Koopman (multilevel analysis) and Don van den Berg and Hendrik Straat (all C and C++ codes). 
Thanks are due 
 to Michael Allerhand, Geert H. van Kollenburg, Renske E. Kuijpers, Rudy Ligtvoet, Hannah E. M. Oosterhuis, Daniel W. van der Palm, and Max Welz for contributing R code; 
 to Geert H. van Kollenburg, Patrick Mair, and Don van Ravenswaaij for testing the software; 
 to Wijbrandt van Schuur for comments on the vignette;
 to Michael Allerhand, Stephen Cubbellotti, Michael Dewey, Jasmin Durstin, Wilco H. M. Emons, Jue Huang, Michael Kubovy, Ivo Molenaar, Jonathan Rose, Tobias Schlaffer, Klaas Sijtsma, Iris Smits, Jia Jia Syu, Stefan Vermeent, Roger Watson, Max Welz, and Na Yang for reporting comments or bugs;
 to Diederick Stoel (ProfitWise) for financial support, 
 to Samantha Bouwmeester, Pierre Cavalini, Johan Denollet, Gina Rossi, Harrie C. M. Vorst,  for permission to use their data;
 to Robert J. Mokken for lending his last name.

Version 0 was introduced in Van der Ark (2007). It included the functions
\tabular{ll}{
\code{\link{coefH}}              \tab Scalability coefficients \cr
\code{\link{coefZ}}              \tab Test statistics for scalability coefficients \cr
\code{\link{check.monotonicity}} \tab Investigate monotonicity assumptions \cr
\code{\link{check.restscore}}    \tab Investigate nonintersection assumption using Method Restscore \cr
\code{\link{check.pmatrix}}      \tab Investigate nonintersection assumption using Method Pmatrix \cr
\code{\link{search.normal}}      \tab Mokken's automated item selection algorithm \cr
}

The following major modifications have been made.
\tabular{ll}{
\code{\link{aisp}}               \tab  More general automated item selection algorithm. \cr
                                 \tab  Function \code{search} has become obsolete (Version 2.0)\cr
\code{\link{check.reliability}}  \tab  Compute reliability coefficients (Version 2.0)\cr
\code{\link{check.iio}}          \tab  Investigate invariant item orderings (Version 2.4) \cr
\code{\link{coefH}}              \tab  Standard errors for scalability coefficients included (Version 2.6)\cr
}

All updates until version 2.7 are described in Van der Ark (2012). The following modifications have been made in Version 2.7 in comparison to previous versions.
\tabular{ll}{
\code{\link{check.errors}}       \tab  Inclusion new function to compute weighted Guttman errors for each person. \cr
\code{\link{check.iio}}          \tab  \code{plot} has been added. \cr
\code{\link{check.monotonicity}} \tab  Computation of number of active pairs for dichotomous items has been corrected. \cr
\code{\link{check.pmatrix}}      \tab  Summary of the results has been corrected. \cr
\code{\link{check.restscore}}    \tab  Code pertaining to IIO has been deleted. The procedure is now equivalent to MSP. \cr
\code{\link{coefH}}              \tab  Option included to compare scalability coefficients across groups\cr
}

The following modifications have been made in Version 2.7.1 in comparison to previous versions.
\tabular{ll}{
\code{\link{mokken}}             \tab  Some legal issues \cr
}

The following modifications have been made in Version 2.7.2 in comparison to previous versions.
\tabular{ll}{
\code{\link{check.iio}}          \tab  Violations of IIO for dichotomous items are now tested using a z-test rather than a t-test. \cr
}

The following modifications have been made in Version 2.7.3 in comparison to previous versions.
\tabular{ll}{
\code{\link{plot.iio.class}}     \tab  Confidence envelopes around estimated response functions \cr
\code{\link{plot.monotonicity.class}}     \tab  Confidence envelopes around estimated response functions \cr
\code{\link{plot.restscore.class}}     \tab  Confidence envelopes around estimated response functions \cr
}

The following modifications have been made in Version 2.8.1 in comparison to previous versions.
\tabular{ll}{
\code{\link{aisp}}     \tab  Startsets have been added \cr
}

The following modifications have been made in Version 2.8.2 in comparison to previous versions.
\tabular{ll}{
\code{\link{recode}}          \tab  New \cr
\code{\link{check.ca}}        \tab  New \cr
\code{\link{check.norms}}     \tab  New \cr
\code{\link{check.errors}}    \tab  Outlier score O+ has been included \cr
}

The following modifications have been made in Version 2.8.3 in comparison to previous versions.
\tabular{ll}{
\code{\link{twoway}}          \tab  New \cr
\code{\link{DS14}}            \tab  New data set\cr
\code{\link{check.errors}}    \tab  Outlier cutoff scores have been included \cr
}

The following modifications have been made in Version 2.8.4 in comparison to previous versions.
\tabular{ll}{
\code{\link{check.iio}}       \tab  New code for computing HT for large samples\cr
}

The following modifications have been made in Version 2.8.5 in comparison to previous versions.
\tabular{ll}{
\code{\link{MLcoefH}}             \tab New code for computing two-level scalability coefficients and standard errors\cr
\code{\link{autonomySupport}}     \tab New two-level data set.\cr
}

The following modifications have been made in Version 2.8.9 and 2.8.10 in comparison to previous versions.
\tabular{ll}{
\code{\link{coefH}}              \tab Included possibility to include a fixed item-step order\cr
\code{\link{MLcoefH}}            \tab Code updated\cr
\code{\link{check.errors}}       \tab Code updated\cr
}

The following modifications have been made in Version 2.8.11 in comparison to previous versions.
\tabular{ll}{
\code{\link{plot}}               \tab The level of tranparancy of the plotted confidence intervals can be adjusted manually \cr
\code{\link{MLcoefH}}            \tab Code updated\cr
}

The following modifications have been made in Version 2.8.12 in comparison to previous versions.
\tabular{ll}{
\code{\link{check.monotonicity}} \tab Z statistic adjusted (Molenaar & Sijtsma, 2000. p. 72 )\cr
\code{\link{check.norms}}        \tab Z Output corrected for \code{nice.output = FALSE}\cr
}

The following modifications have been made in Version 2.9.0 in comparison to previous versions.
\tabular{ll}{
\code{\link{coefH}}              \tab Z Solution of Koopman et al. (2017) implemented to solve the problem of equal item steps and code updated\cr
\code{\link{MLcoefH}}            \tab Z Solution of Koopman et al. (2017) implemented to solve the problem of equal item steps and code updated\cr
}

The following modifications have been made in Version 3.0.0 in comparison to previous versions.
\tabular{ll}{
\code{\link{aisp}}               \tab Genetic algorithm has been reprogrammed and is now much faster.\cr
                                 \tab Argument \code{lowerbound} can now be a vector, enabling the investigation of several lower bounds simultaneously.\cr
                                 \tab Extra argument to specify which standard errors should be used in Z-test (Koopman et al., 2020b).\cr
                                 \tab Extra argument to indicate which null-hypothesis should be used to test Hi (i.e., Hi = c or Hi = 0)(Koopman et al., 2020b).\cr
                                 \tab Extra argument to handle nested data (Koopman et al., 2020b).\cr
\code{\link{check.iio}}          \tab Computation of Coefficient HT for large samples is now much faster.\cr
\code{\link{coefH}}              \tab New standard errors for nested data (Koopman et al. 2020a).\cr
\code{\link{coefZ}}              \tab Extra argument to compute the Z-score using \code{lowerbound} as the null hypothesis (Koopman et al., 2020b).\cr
                                 \tab Extra argument to compute Z-score using delta method standard error, but the original method remains available (Koopman et al., 2020b).\cr 
                                 \tab Extra argument to compute Z-score in nested data (Koopman et al., 2020b).\cr
\code{\link{ICC}}                \tab New function for ICCs in two-level Mokken scale analysis (Koopman et al. 2020a)\cr
\code{\link{MLcoefH}}            \tab Extra argument for weighted proportions. Reduces bias in two-level standard errors (Koopman et al. 2020a)\cr
                                 \tab Extra argument for a fixed item-step order\cr
\code{\link{SWMD}}               \tab New data file (Koopman et al. 2020a)\cr
}

The following modifications have been made in Version 3.0.3 in comparison to previous versions.
\tabular{ll}{
\code{\link{coefZ}}              \tab Error handling added\cr
                                 \tab The delta test uses range-preserving asymptotic theory (Koopman, Zijlstra, & Van der Ark, 2020c) \cr
\code{\link{coefH}}              \tab Range-preserving confidence intervals added (Koopman et al., 2020c)\cr
                                 \tab Extra argument to print variance-covariance matrices of estimated coefficients\cr
                                 \tab Error handling added for more than 10 response categories\cr
\code{\link{MLcoefH}}            \tab Range-preserving confidence intervals added (Koopman et al., 2020c)\cr
                                 \tab Extra argument to print variance-covariance matrices of estimated coefficients\cr
\code{\link{MLcoefZ}}            \tab New function for z-scores of two-level scalability coefficients (Koopman et al., 2020c)\cr
\code{\link{aisp}}               \tab type.se default adjusted and additional error handling added\cr
\code{\link{check.errors}}       \tab Repaired bug in check.errors\cr
\code{\link{mcmi}}               \tab New data file (Sijtsma & van der Ark, 2020)\cr
All functions                    \tab A warning has been added if items have different numbers of response categories       
}

The following modifications have been made in Version 3.0.4 in comparison to previous versions.
\tabular{ll}{
\code{\link{aisp}}               \tab Argument type.z replaces type.se to accommodate three types of z scores (Mokken's Z, Wald-based, and range-preserving) \cr
\code{\link{search.normal}}      \tab Argument type.z replaces type.se to accommodate three types of z scores (Mokken's Z, Wald-based, and range-preserving) \cr
\code{\link{coefZ}}              \tab Argument type.z replaces type.se to accommodate three types of z scores (Mokken's Z, Wald-based, and range-preserving)\cr
\code{\link{MLcoefZ}}            \tab Argument type.z is added to accommodate two types of z scores (Wald-based and range-preserving) \cr
\code{\link{coefH}}              \tab Argument type.ci is added to accommodate two types of confidence intervals (Wald-based and range-preserving) \cr
				 \tab Argument print.to.screen is replaced by results \cr
\code{\link{MLcoefH}}            \tab Argument type.ci is added to accommodate two types of confidence intervals (Wald-based and range-preserving) \cr
\code{\link{SWMDK}}              \tab New data file (Koopman, Zijlstra, & Van der Ark 2020a)\cr
All functions                    \tab A warning has been added if items have different numbers of response categories       
}

The following modifications have been made in Version 3.0.6 in comparison to previous versions.
\tabular{ll}{
\code{\link{coefH}}               \tab Fixed a bug for Linux environment
}
}

\author{
L. Andries van der Ark & L. Koopman
Maintainer: L. Andries van der Ark <L.A.vanderArk@uva.nl>. 
}

\references{

Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020a). 
\emph{A two-step, test guided Mokken scale analysis for nonclustered and clustered data}. 
Manuscript submitted for publication.

Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2020b).
New advances in Mokken scale analysis.
Paper presented at the online meeting of the Psychometric Society, July, 2020.
  
Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2020c).
Range-preserving confidence intervals for scalability coefficients in Mokken scale analysis.
Manuscript in preparation.
  
Mokken, R. J. (1971) 
\emph{A Theory and Procedure of Scale Analysis}. 
De Gruyter.

Molenaar, I.W., & Sijtsma, K. (2000) 
\emph{User's Manual MSP5 for Windows} [Software manual].
IEC ProGAMMA.
  
Sijtsma, K., & Molenaar, I. W. (2002)
\emph{Introduction to nonparametric item response theory}.
Sage.

Sijtsma, K., & Van der Ark, L. A. (2017).
A tutorial on how to do a Mokken scale analysis on your test and questionnaire data.
\emph{British Journal of Mathematical and Statistical Psychology, 70}, 137-158.
\doi{10.1111/bmsp.12078}
  
Sijtsma, K., & Van der Ark, L. A. (2020),
\emph{Measurement models for psychological attributes.}
Boca Raton, FL: Chapman and Hall/CRC Press.
\url{https://www.routledge.com/Measurement-Models-for-Psychological-Attributes/Sijtsma-Ark/p/book/9780367424527}

Van der Ark, L. A. (2007). 
Mokken scale analysis in \pkg{R}. 
\emph{Journal of Statistical Software, 20}(11), 1-19.
\url{https://www.jstatsoft.org/article/view/v020i11}

Van der Ark, L. A. (2010). 
Getting started with Mokken scale analysis in \pkg{R}. 
Unpuablished manuscript.
\url{https://sites.google.com/a/tilburguniversity.edu/avdrark/mokken}

Van der Ark, L. A. (2012). 
New developments in Mokken scale analysis in \pkg{R}. 
\emph{Journal of Statistical Software, 48}(5), 1-27.
\url{https://www.jstatsoft.org/article/view/v048i05}
}

\examples{

# Personality test
data(acl)               

# Select the items of the scale Communality
Communality <- acl[,1:10]

# Compute scalability coefficients
coefH(Communality)

# Investigate the assumption of monotonicity
monotonicity.list <- check.monotonicity(Communality)
summary(monotonicity.list)
plot(monotonicity.list)

# Investigate the assumption of non-intersecting ISRFs using method restscore
restscore.list <- check.restscore(Communality)
summary(restscore.list)
plot(restscore.list)

# Investigate the assumption of non-intersecting ISRFs using method pmatrix
pmatrix.list <- check.pmatrix(Communality)
summary(pmatrix.list)
plot(pmatrix.list)

# Investigate the assumption of IIO using method MIIO
iio.list <- check.iio(Communality)
summary(iio.list)
plot(iio.list)

# Compute the reliability of the scale
check.reliability(Communality)

# Partition the the scale into mokken scales
aisp(Communality)

}
