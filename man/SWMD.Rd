\name{SWMD}
\alias{SWMD}
\docType{data}
\title{SWMD Data Subset}

\description{
    A subset of the COOL5-18 data (see below) with scores of 651 pupils nested in 30 classes
    on the 7-item Schaal Welbevinden Met Docenten [Scale Well-Being With Teachers]
    (Peetsma, Wagenaar, & de Kat, 2001; 
    Zijsling, Keuning, Keizer-Mittelhauser, Naaijer, & Timmermans, 2017).
    See \code{attributes(SWMD)} for the original item labels and content.
    R-code to get this subsample is available as online supplement to 
    Koopman, Zijlstra, & Van der Ark (2020).
}
\usage{data(SWMD)}
\format{
    A 651 by 8 data frame containing integers. The first column reflects a classroom indicator, 
    the remaining columns the 7 items, see \code{colnames(SWMD)}.
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
    }

    The items were translated from Dutch. For the original items, see p. 79 in Zijsling et al. (2017). The scores on these items plus seven additional items are available in dataset \code{\link{SWMDK}}.
    
}

\source{
    Data is a subset of respondents and items from the large-scale cohort study COOL5-18. 
    \url{http://www.cool5-18.nl/} (Zijsling et al., 2017).
    For entire dataset see \doi{10.17026/dans-zfp-egnq}.
    Dataset is accessible after login.
}

\references{

Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020). 
   \emph{A two-step procedure for scaling multilevel data using Mokken's scalability coefficients}. Manuscript submitted 
     for publication.

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
\code{\link{MLcoefH}}, 
\code{\link{ICC}}, 
}

\examples{
# Data example (Koopman et al., 2020)
data(SWMD)

# Item content, see labels
attributes(SWMD)$labels

# Compute ICC
ICC(SWMD)
}

\keyword{datasets}
