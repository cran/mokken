\name{transreas}
\alias{transreas}
\docType{data}
\title{Transitive Reasoning}

\description{
Data came from 12 dichtomous items administered to 425 children in grades 2 to 6 (Verweij, Sijtsma, and Koops, 1996).
Each item is a transitive reasoning task.}
\usage{data(transreas)}
\format{
A 425 by 13 (grade and scores on 12 items) matrix containing integers. 
  \code{attributes(transreas)} gives details on the items
}
\details{
 Items have two ordered answer categories,
\emph{incorrect} (score 0), \emph{correct} (1).
(for detailed information, see Sijtsma and Molenaar, 2002, p. 33).


\tabular{lrlclr}{
Item \tab Task \tab Property\tab             Format\tab  Objects\tab                   Measures \cr
T09L \tab   9  \tab  length \tab YA = YB < YC = YD \tab  sticks \tab   12.5, 12.5, 13, 13 (cm)  \cr
T12P \tab  12  \tab  pseudo \tab                   \tab         \tab                            \cr
T10W \tab  10  \tab  weight \tab YA = YB < YC = YD \tab   balls \tab      60, 60, 100, 100 (g)  \cr
T11P \tab  11  \tab  pseudo \tab                   \tab         \tab                            \cr
T04W \tab   4  \tab  weight \tab YA = YB = YC = YD \tab   cubes \tab                    65 (g)  \cr
T05W \tab   5  \tab  weight \tab      YA < YB < YC \tab   balls \tab           40, 50, 70 (cm)  \cr
T02L \tab   2  \tab  length \tab YA = YB = YC = YD \tab   tubes \tab                   12 (cm)  \cr
T07L \tab   7  \tab  length \tab      YA > YB = YC \tab  sticks \tab     28.5, 27.5, 27.5 (cm)  \cr
T03W \tab   3  \tab  weight \tab      YA > YB > YC \tab   tubes \tab            45, 25, 18 (g)  \cr
T01L \tab   1  \tab  length \tab      YA > YB > YC \tab  sticks \tab         12, 11.5, 11 (cm)  \cr
T08W \tab   8  \tab  weight \tab      YA > YB = YC \tab   balls \tab            65, 40, 40 (g)  \cr
T06A \tab   6  \tab    area \tab      YA > YB > YC \tab   discs \tab7.5, 7, 6.5 (diameter; cm)  \cr
}
}

\references{

Verweij, A. C., Sijtsma, K., and Koops, W. (1996). 
A Mokken scale for transitive reasoning suited for longitudinal research.
\emph{International Journal of Behavioral Development, 23}, 241-264.

Sijtsma, K. and Molenaar, I. W. (2002). 
\emph{Introduction to nonparametric item response theory}, 
Thousand Oaks, CA: Sage.

}

\examples{
# Construction of Table 3.1 in Sijtsma and Molenaar (2002, p. 33)
data(transreas)
grades <- transreas[,1]
item.scores <- transreas[,-1]
Total.group <- round(apply(item.scores,2,mean),2)
for (i in 2:6) assign(paste("Grade.",i,sep=""),
   round(apply(item.scores[grades==i,],2,mean),2))
Task <- c(9,12,10,11,4,5,2,7,3,1,8,6)
Property <- attributes(transreas)$property
Format <- attributes(transreas)$format
Objects <- attributes(transreas)$objects
Measures <- attributes(transreas)$measures
Table.3.1 <- data.frame(Task,Property,Format,Objects,Measures,
   Total.group,Grade.2,Grade.3,Grade.4,Grade.5,Grade.6)
Table.3.1
}

\keyword{datasets}
