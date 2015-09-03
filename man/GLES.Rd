\name{GLES}
\alias{GLES}
\docType{data}
\title{
German Longitudinal Election Study (GLES)
}
\description{
Data from the German Longitudinal Election Study (GLES), see Rattinger et al. (2014). The GLES is a long-term study of the German electoral process. It collects pre- and post-election data for several federal elections, the data used here originate from the pre-election study for 2013. 
}
\usage{data("GLES")}
\format{
  A data frame with 1155 observations on the following 18 variables.
  \describe{
    \item{\samp{SPD vs FDP}}{Ordinal paired comparison between SPD and FDP}
    \item{\samp{SPD vs Left Party}}{Ordinal paired comparison between SPD and Left Party}
    \item{\samp{SPD vs Greens}}{Ordinal paired comparison between SPD and Greens}
    \item{\samp{SPD vs CDU_CSU}}{Ordinal paired comparison between SPD and CDU/CSU}
    \item{\samp{FDP vs Left Party}}{Ordinal paired comparison between FDP and Left Party}
    \item{\samp{FDP vs Greens}}{Ordinal paired comparison between FDP and Greens}
    \item{\samp{FDP vs CDU_CSU}}{Ordinal paired comparison between FDP and CDU/CSU}
    \item{\samp{Left Party vs Greens}}{Ordinal paired comparison between Left Party and Greens}
    \item{\samp{Left Party vs CDU_CSU}}{Ordinal paired comparison between Left Party and CDU/CSU}
    \item{\samp{Greens vs CDU_CSU}}{Ordinal paired comparison between Greens and CDU/CSU}
    \item{\code{Age}}{Age in years}
    \item{\code{Gender}}{0: male, 1: female}
    \item{\code{EastWest}}{0: West Germany, 1:East Germany}
    \item{\code{PersEcon}}{Personal economic situation, 1: good or very good, 0: else}
    \item{\code{Abitur}}{School leaving certificate, 1: Abitur/A levels, 0: else}
    \item{\code{Unemployment}}{1: currently unemployed, 0: else}
    \item{\code{Church}}{Frequency of attendence in a church/synagogue/mosque/..., 1: at least once a month, 0: else}
    \item{\code{Migration}}{Have you been a German citizen since birth? 1: yes, 0: no}
  }
}
\details{
Variables 1 to 10 represent the response, variables 11 to 18 represent the subject-specific covariates. The response variables are ordinal, with values from 1 to 5. Low values represent string preference of the first-names party, high values represent strong preference of the last-named party. 
}
\source{
\url{http://www.gesis.org/en/elections-home/gles/data-and-documents/}
}
\references{
Rattinger, H., S. Rossteutscher, R. Schmitt-Beck, B. Wessels, and C. Wolf (2014):
Pre-election cross section (GLES 2013). \emph{GESIS Data Archive, Cologne ZA5700
Data file Version 2.0.0.}

Schauberger, Gunther and Tutz, Gerhard (2015): Modelling Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an Application to Party Preference Data, \emph{Department of Statistics, LMU Munich}, Technical Report 183
}
\examples{
data(GLES)
}
\keyword{datasets}
