\name{trog}
\alias{trog}
\docType{data}
\title{trog Data}

\description{
    A clustered-item dataset with scores of 210 children (measurement taken at a first time 
point of a larger developmental study; Brinchmann et al., 2019) 
that took the Norwegian adaptation of the Test for Reception of Grammar (TROG; Bishop, 1979). 
The TROG consists of C=20 clusters of Jc=4 dichotomously scored items for all c. 
and all kids were administered each of the J=80 items. This item set is used as a real data example 
to demonstrate how to investigate the ordering structure of a test in Koopman & Braeken (2024).
}

\usage{data(autonomySupport)}
\format{
    A 210 by 80 data frame containing dichotomous data. Each column reflects an item, each row a respondent.
The column labels reflect the cluster (letters) and item within the cluster (numbers). For example,
item a1 is the first item in cluster a, see \code{colnames(trog)}.
}
\details{
    Each item is dichotomously scored, where 1 reflects a correct response and 0 an incorrect response. 
}


\references{
Bishop, D. V. M. (1979). Comprehension in developmental language
disorders. \emph{Developmental Medicine & Child Neurology, 21(2)},
225-238. \doi{10.1111/j.1469-8749.1979.tb01605.x}.

Brinchmann, E. I., Braeken, J., & Lyster, S.-A. H. (2019). Is
there a direct relation between the development of vocabulary and
grammar? \emph{Developmental Science, 22(1)}, 1-13. 
\doi{10.1111/desc.12709}.

Koopman, L. & Braeken, J. (2024). 
  Investigating the Ordering Structure of Clustered Items Using Nonparametric 
  Item Response Theory. Manuscript submitted for publication.


}

\seealso{
\code{\link{check.iio}}, 
}

\examples{
data(trog)
}
\keyword{datasets}
