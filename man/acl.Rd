\name{acl}
\alias{acl}
\docType{data}
\title{Adjective Checklist Data}

\description{
 Scores of 433 students on 218 items from a Dutch version of the Adjective Checklist.
}
\usage{data(acl)}
\format{
A 433 by 218 matrix containing integers. \code{dimnames(acl)[[2]]} are adjectives
}
\details{
 Each item is an adjective with five ordered answer categories 
 (0 = completely disagree, 1 = disagree, 2 = agree nor disagree, 3 = agree, 4 = completely agree).
 The respondents were instructed to consider whether an adjective described their 
 personality, and mark the answer category that fits best to this description.
The 218 items constitute 22 scales (see table) 
 77 items of the 218 items that constitute the ten scales were negatively worded.
The negatively worded items are indicated by an asterisk in the \code{dimnames}
 and their item scores were reversed.

\tabular{llcll}{
Communality \tab Items 1-10    \tab \tab Change              \tab  Items 111-119 \cr
Achievement \tab Items 11-20   \tab \tab Succorance          \tab  Items 121-129 \cr
Dominance   \tab Items 21-30   \tab \tab Abasement           \tab  Items 131-139 \cr
Endurance   \tab Items 31-40   \tab \tab Deference           \tab  Items 141-149 \cr
Order       \tab Items 41-50   \tab \tab Personal Adjustment \tab  Items 151-159 \cr
Intelligence\tab Items 51-60   \tab \tab Ideal Self          \tab  Items 161-169 \cr
Nurturance  \tab Items 61-70   \tab \tab Critical parent     \tab  Items 171-179 \cr
Affiliation \tab Items 71-80   \tab \tab Nurturant parent    \tab  Items 181-189 \cr
Exhibition  \tab Items 81-90   \tab \tab Adult               \tab  Items 191-199 \cr
Autonomy    \tab Items 91-100  \tab \tab Free Child          \tab  Items 200-209 \cr
Aggression  \tab Items 101-110 \tab \tab Adapted Child       \tab  Items 210-218 }

}
\source{
 Data were kindly made available by H. C. M. Vorst from the University of Amsterdam.
 The original Adjective Checklist was developed by Gough and Heilbrun (1980).
}
\references{
Gough, H. G. and Heilbrun,A. B. (1980) 
 \emph{The Adjective Check List, Manual 1980 Edition}.
 Palo Alto, CA: Consulting Psychologists Press.

Van der Ark, L. A. (2007) 
Mokken scale analysis in \pkg{R}. 
\emph{Journal of Statistical Software}.
\url{http://www.jstatsoft.org}

}

\examples{
data(acl)
}
\keyword{datasets}
