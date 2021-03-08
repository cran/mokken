\name{mcmi}
\alias{mcmi}
\docType{data}
\title{Millon Clinical Multiaxial Inventory}

\description{
 Scores of 1208 patients and inmates on 44 dichotomous items from a Dutch version of the Millon Clinical Multiaxial Inventory-III.
}

\usage{
data(mcmi)
}

\format{
A 1208 by 44 data frame containing integers. 
}

\details{
 The data were collected by Gina Rossi (Rossi et al., 2010) as part of a larger project. This subset of 44 items was used to 
 demonstrate diagnostic classification models (cognitive diagnosis models) by de la Torre et al. (2018), Van der Ark et al. 
 (2019), and Sijtsma & Van der Ark (2020). The Q matrix used in the analyses is an attribute. Both the items and the 
 respondents have been deidentified. The deidentified data do allow to replicate the analyses in Sijtsma and Van der Ark (2020, 
 chapter 5), using the code available from \url{https://osf.io/e9jrz}. 
}

\source{
 Data were kindly made available by Gina Rossi from the Vrije Universiteit Brussel, Belgium.
 The original Millon Clinical Multiaxial Inventory-III was developed by Millon (1994).
}

\references{
de la Torre, J., Van der Ark, L. A., & Rossi, G. (2018).
Analysis of clinical data from cognitive diagnosis modelling framework.
\emph{Measurement and Evaluation in Counseling and Development, 51}, 281-296.
\doi{10.1080/07481756.2017.1327286}

Millon, T. (1994). \emph{Manual for the MCMI-III}. Minneapolis, MN: National Computer Systems.

Rossi, G., Elklit, A., & Simonsen, E. (2010). 
Empirical evidence for a four factor framework of personality disorder organization: 
Multigroup confirmatory factor analysis of the million clinical multiaxial inventory-III personality disorders scales across Belgian and Danish data samples.
\emph{Journal of Personality Disorders, 24}, 128-150. 
\doi{10.1521/pedi.2010.24.1.128} 

Sijtsma, K., & Van der Ark, L. A. (2020),
\emph{Measurement models for psychological attributes.} Chapman and Hall/CRC Press.
\url{https://www.routledge.com/Measurement-Models-for-Psychological-Attributes/Sijtsma-Ark/p/book/9780367424527}

Van der Ark, L. A., Rossi, G., & Sijtsma, K. (2019).
Nonparametric item response theory and Mokken scale analysis, with relations to latent class models and cognitive diagnostic models.
In M. von Davier & Y.- S. Lee (Eds.), \emph{Handbook of diagnostic classification models} (pp. 21-45). Springer.
\doi{10.1007/978-3-030-05584-4}
}

\examples{
data(mcmi)

# Q matrix
attributes(mcmi)$Q
}

\keyword{datasets}
