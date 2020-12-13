\name{transreas2}
\alias{transreas2}
\docType{data}
\title{Transitive Reasoning Data}

\description{
 Scores of 606 school children on 16 dichotomous transitive reasoning items.
}

\usage{
data(transreas2)
}

\format{
A 606 by 16 data frame containing integers. 
}

\details{
 The data were collected by Samantha Bouwmeester (Bouwmeester & Sijtsma, 2004). The design of the items is included as attributes. The 
 respondents have been deidentified. The deidentified data do allow to replicate the analyses in Sijtsma and Van der Ark (2020, 
 chapter 4), using the code available from \url{https://osf.io/e9jrz}. Note that the package \code{mokken} includes another transitive reasoning data set. 
}

\source{
 Data were kindly made available by Samantha Bouwmeester from the Erasmus Universiteit Rotterdam, The Netherlands.
}

\references{
Bouwmeester, S., & Sijtsma, K. (2004). 
Measuring the ability of transitive reasoning, using product and strategy information. 
\emph{Psychometrika, 69}, 123-146. 
\url{https://doi.org/10.1007/BF02295843}

Sijtsma, K., & Van der Ark, L. A. (2020),
\emph{Measurement models for psychological attributes.} Chapman and Hall/CRC Press.
\url{https://www.routledge.com/Measurement-Models-for-Psychological-Attributes/Sijtsma-Ark/p/book/9780367424527}

}

\examples{
data(transreas2)

# Create Table 4.2 from Sijtsma & Van der Ark (2020)
data.frame (Item = 1: 16, 
            taskContent = attr(transreas2, which = "taskContent"), 
            presentationForm = attr(transreas2, which = "presentationForm"), 
            taskFormat = attr(transreas2, which = "taskFormat"),
            pValues = round(apply(transreas2, 2, mean), 2))
}
\keyword{datasets}
