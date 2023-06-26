\name{SWMDK}
\alias{SWMDK}
\docType{data}
\title{SWMDK Data Subset}

\description{
    A subset of the COOL5-18 data (see below) with scores of 639 pupils nested in 30 classes
    on the 7-item Schaal Welbevinden Met Docenten [Scale Well-Being With Teachers] and 6-item 
    Schaal Welbevinden met Klasgenoten [Scale Well-Being With Classmates]
    (Peetsma et al., 2001; Zijsling et al., 2017).
    See \code{attributes(SWMDK)} for the original item labels and content.
    R-code to get this subsample is available as online supplement to 
    Koopman et al. (in press a).
}
\usage{data(SWMDK)}
\format{
    A 639 by 14 data frame containing integers. The first column reflects a classroom indicator, 
    the remaining columns the 13 items, see \code{colnames(SWMDK)}.
}
\details{
    Each item has five ordered answer categories from \emph{not true at all} (score 0) 
    to \emph{completely true} (score 4). 
    
    \tabular{lll}{
    \emph{Item} \tab \emph{Short} \tab \emph{Content} \cr
    Item1  \tab lv_wdo1 \tab The teachers usually know how I feel \cr
    Item2  \tab lv_wdo2 \tab I can talk about problems with the teachers  \cr   
    Item3  \tab lv_wdo3 \tab If I feel unhappy, I can talk to the teachers about it  \cr   
    Item4  \tab lv_wdo4 \tab I feel at ease with the teachers \cr   
    Item5  \tab lv_wdo5 \tab The teachers understand me \cr   
    Item6  \tab lv_wdo6 \tab I have good contact with the teachers \cr   
    Item7  \tab lv_wdo7 \tab I would prefer to have other teachers (inversely coded) \cr   
    Item8  \tab lv_wkl1 \tab I have a lot of contact with my classmates \cr                 
    Item9  \tab lv_wkl2 \tab I would prefer to be in a different class (inversely coded) \cr
    Item10 \tab lv_wkl3 \tab We have a nice class \cr
    Item11 \tab lv_wkl4 \tab I get along well with my classmates \cr                      
    Item12 \tab lv_wkl5 \tab I sometimes feel alone in the class (inversely coded) \cr     
    Item13 \tab lv_wkl6 \tab I enjoy hanging out with my classmates \cr 
    }

    The items were translated from Dutch. For the original items, see pp. 79-83 in Zijsling et al. (2017). The first seven items are also available in dataset \code{\link{SWMD}}.
}

\source{
    Data is a subset of respondents and items from the large-scale cohort study COOL5-18. 
    \url{http://www.cool5-18.nl/} (Zijsling et al., 2017).
    For entire dataset see \doi{10.17026/dans-zfp-egnq}
    Dataset is accessible after login.
}

\references{

   Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (in press a).
   A two-step, test-guided Mokken scale analysis for nonclustered and clustered data. 
   \emph{Quality of Life Research}. (advanced online publication)
   \doi{10.1007/s11136-021-02840-2}

Peetsma, T. T. D., Wagenaar, E., & De Kat, E. (2001). 
   School motivation, future time perspective and well-being of high school students in segregated 
   and integrated schools in the Netherlands and the role of ethnic self-description. 
   In J. Koppen, I. Lunt, & C. Wulf (Eds.), 
   \emph{Education in Europe. Cultures, Values, Institutions in Transition} 
   (pp. 54-74). Waxmann.

Zijsling, D., Keuning, J., Keizer-Mittelhaeuser, M.-A., Naaijer, H., & Timmermans, A. (2017).
   \emph{Cohortonderzoek COOL5-18: Technisch rapport meting VO-3 in 2014}. 
   Onderwijs/Onderzoek. 
}

\seealso{
\code{\link{coefH}},
\code{\link{MLcoefH}}, 
\code{\link{ICC}}, 
\code{\link{SWMD}} 
}

\examples{
# Data example (Koopman et al., in press)
data(SWMDK)

# Item content, see labels
attributes(SWMDK)$labels

# Compute ICC
ICC(SWMDK)
}

\keyword{datasets}
