\name{DS14}
\alias{DS14}
\docType{data}
\title{DS14}

\description{
 Gender, age, and item scores on the DS14 questionnaire of 541 coronary artery disease patients.
}

\usage{data(DS14)}

\format{A 541 by 16 matrix containing gender, age, and item scores on the DS14 questionnaire.}

\details{
The DS14 (Denollet, 2005) is the most accepted and widely used diagnostic instrument for the assessment of the type-D pattern.
Type D (distressed) is defined as the joint tendency towards negative affectivity (e.g., worry, irritability, gloom) and 
social inhibition (e.g., reticence and a lack of self-assurance).
DS14 contains 14 items, each having five ordered response categories 
(0 = completely disagree, 1 = disagree, 2 = agree nor disagree, 3 = agree, 4 = completely agree).
Items 2, 4, 5, 7, 9, 12, and 13 measure negative affectivity.
Items 1, 3, 6, 8, 10, 11, and 14 measure social inhibition. Items 1 and 3 are negatively worded (indicated by an asterisk in the \code{dimnames}).

The data contain the gender (\code{Male}) of the patients (1 = male, 0 =  female), the age (\code{Age}) of the patients in years,
and the scores to DS14. Ten item scores are missing. Items 1 and 3 must be recoded before the data can beused meaningfully.
The data have been used to investigate predictive value of social inhibition and negative affectivity for 
cardiovascular events and mortality in patients with coronary artery disease (Denollet et al., 2013),
to investigate the relation between Type D and inflammation and endothelial dysfunction (van Dooren et al., 2016),
and to investigate the relation between Type D and increased macrophage activity (Zuccarella-Hackl et al., 2016).
These data have also been analyzed in papers on Mokken scale analysis (Sijtsma & Van der Ark, 2016; Straat et al., 2016).
}

\source{
 Data were kindly made available by J. Denollet from Tilburg University.
}

\references{
 Denollet, J., Pedersen, S. S., Vrints, C. J., & Conraads, V. M. (2013). 
 Predictive value of social inhibition and negative affectivity for cardiovascular events and mortality in patients with coronary artery disease: the Type D personality construct. 
 \emph{Psychosomatic Medicine, 75}, 873-981.
  
 van Dooren, F. E., Verhey, F. R., Pouwer, F., Schalkwijk, C. G., Sep, S. J., Stehouwer, C.D., Henry, R. M., Dagnelie, P. C., 
 Schaper, N. C., Van der Kallen, C. J., Koster, A., Schram, M. T., & Denollet, J. (2016).
 Association of Type D personality with increased vulnerability to depression: Is there a role for inflammation or endothelial dysfunction? - The Maastricht Study. 
 \emph{Journal of Affective Disorders, 189}, 118-125. 
 \url{https://doi.org/10.1016/j.jad.2015.09.028}

 Sijtsma, K., & Van der Ark, L. A. (2017).
 A tutorial on how to do a Mokken scale analysis on your test and questionnaire data.
 \emph{British Journal of Mathematical and Statistical Psychology, 70}, 137-158.
 \url{https://doi.org/10.1111/bmsp.12078}

 Straat, J. H., Van der Ark, L. A., & Sijtsma, K. (2016).
 Using conditional association to identify locally Independent item sets.
 \emph{Methodology, 12}, 117-123.
 \url{https://doi.org/10.1027/1614-2241/a000115}
   
 Zuccarella-Hackl, C., von Kaenel, R., Thomas, L., Kuebler, P., Schmid, J. P., Mattle. H. P., Mono, M. L., Rieben, R., 
 Wiest, R., & Wirtz, P. H. (2016). 
 Higher macrophage superoxide anion production in coronary artery disease (CAD) patients with Type D personality. 
 \emph{Psychoneuroendocrinology, 68}, 186-193.
 \url{https://doi.org/10.1016/j.psyneuen.2016.02.031}
}

\seealso{
 \code{\link{recode}}, \code{\link{twoway}}
}

\examples{
data(DS14)

# Handle missing data and recode negatively worded items
X <- DS14[, 3 : 16]
X <- twoway(X)
X <- recode(X, c(1, 3))
head(X)
}

\keyword{datasets}
