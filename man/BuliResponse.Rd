% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BTLLasso-package.R
\docType{data}
\name{BuliResponse}
\alias{BuliResponse}
\title{Bundesliga Data Response Data (BuliResponse)}
\format{
A data set containing all information that is necessary to create a response object
for the Bundesliga data \code{link{Buli1516}}
\describe{ 
\item{Result}{Ordinal, 5-categorical results from Bundesliga season 2015/16.}
\item{TeamHome}{Abbreviation of home team.}
\item{TeamAway}{Abbreviation of away team.}
\item{Matchday}{Matchdays from 1 to 34.}
}
 @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
Statistical Software}, to appear

Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
17(3), 223 - 243

Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
Analysis of the importance of on-field covariates in the German Bundesliga, 
\emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
}
\source{
\url{https://www.kicker.de/}
}
\description{
Data from the German Bundesliga from the season 2015/16. The data contain all 
variables from the 306 matches that are necessary to create the respective 
\code{response.BTLLasso} object from the data set \code{\link{Buli1516}}.  The purpose
of the data set is to provide an example how \code{response.BTLLasso} objects can be created.
}
\examples{
\dontrun{
data(BuliResponse)

Y.Buli <- response.BTLLasso(response = BuliResponse$Result, 
                            first.object = BuliResponse$TeamHome,
                            second.object = BuliResponse$TeamAway,
                            subject = BuliResponse$Matchday)
}
}
\keyword{datasets}
