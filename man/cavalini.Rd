\name{cavalini}
\alias{cavalini}
\docType{data}
\title{Coping Strategies}

\description{
Data came from 17 polytomous items administered to 828 respondents (Cavalini, 1992) asking them how they coped actively
 with the bad smell from a factory in the neighborhood of their homes.}
\usage{data(cavalini)}
\format{
A 828 by 17 matrix containing integers. \code{dimnames(cavalini)[[2]]} are adjectives
}
\details{
 Items have four ordered answer categories,
\emph{never} (score 0), \emph{seldom} (1), \emph{often} (2), and \emph{always} (3).
The 17 items constitute 4 scales (for detailed information, see Sijtsma and Molenaar, 2002, pp. 82-86).


\tabular{ll}{
Item1  \tab Keep windows closed  \cr
Item2  \tab No laundry outside  \cr   
Item3  \tab Search source of malodor  \cr   
Item4  \tab No blankets outside  \cr   
Item5  \tab Try to find solutions  \cr   
Item6  \tab Go elsewhere for fresh air  \cr   
Item7  \tab Call environmental agency  \cr   
Item8  \tab Think of something else  \cr   
Item9  \tab File complaint with producer  \cr   
Item10 \tab Acquiesce in odor annoyance \cr   
Item11 \tab Do something to get rid of it \cr   
Item12 \tab Say ``it might have been worse'' \cr   
Item13 \tab Experience unrest \cr   
Item14 \tab Talk to friends and family \cr   
Item15 \tab Seek diversion \cr   
Item16 \tab Avoid breathing through the nose \cr   
Item17 \tab Try to adapt to situation \cr   
}
}

\references{

Cavalini, P. M. (1992). 
\emph{It's an ill wind that brings no good. Studies on odour annoyance and the dispersion of odorant concentrations from industries.} 
Unpublished doctoral disseratation. University of Groningen, The Netherlands.

Sijtsma, K. and Molenaar, I. W. (2002). 
\emph{Introduction to nonparametric item response theory}, 
Thousand Oaks, CA: Sage.

}

\seealso{
\code{\link{check.iio}}, 
}

\examples{
data(cavalini)
}

\keyword{datasets}
