\name{autonomySupport}
\alias{autonomySupport}
\docType{data}
\title{Autonomy Support Data}

\description{
    A two-level dataset with scores of 14 teachers who are rated by a group of pupils 
    on 7 items from a Dutch Autonomy Support questionnaire 
    (group size ranged between 5 and 39 pupils, total number of pupils is 259).
}
\usage{data(autonomySupport)}
\format{
    A 259 by 8 data frame containing integers. The first column reflects a teacher indicator, 
    the remaining columns the 7 items, see colnames(autonomySupport).
}
\details{
    Each item has five ordered answer categories from \emph{not at all/never} (score 1) 
    to \emph{certainly/always} (score 5). The items reflect several autonomy supportive behaviours from teachers. 
    
    \tabular{lll}{
    \emph{Item} \tab \emph{Short} \tab \emph{Content} \cr
    Item 1  \tab Choose \tab The teacher lets me choose what I am going to do  \cr
    Item 2  \tab Decide \tab The teacher decides which task I will start with (inversely coded)  \cr   
    Item 3  \tab Task \tab I get to choose which task I will start with  \cr   
    Item 4  \tab Listen \tab The teacher listens to me when I disagree with something  \cr   
    Item 5  \tab Help \tab The teacher helps me when I ask for it  \cr   
    Item 6  \tab Accept \tab The teacher accepts me for who I am  \cr   
    Item 7  \tab Understand \tab The teacher helps me when I do not understand a task  \cr   
    }
}

\source{
    The seven items are a subset from a self-constructed 27-item questionnaire on teacher's autonomy support. 
    Data were collected and made available by V. E. C. Koopman from the University of Amsterdam.
}

\references{

Koopman, V. E. C., Van der Ark, L. A., & Zijlstra, B. J. H. (2017). 
\emph{Standard errors of two-level scalability coefficients}. Manuscript in preparation.

}


\seealso{
\code{\link{MLcoefH}}, 
}

\examples{
data(autonomySupport)
}
\keyword{datasets}
