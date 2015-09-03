\name{GLESsmall}
\alias{GLESsmall}
\docType{data}
\title{
Subset of the GLES data set with 200 observations and 4 covariates.
}
\description{
This is a subset of the \code{\link{GLES}} data set from the German Longitudinal Election Study (GLES), see Rattinger et al. (2014). The GLES is a long-term study of the German electoral process. It collects pre- and post-election data for several federal elections, the data used here originate from the pre-election study for 2013. 
}
\usage{data("GLESsmall")}
\format{
  A data frame with 200 observations on the following 14 variables.
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
  }
}
\details{
This is a subset of \code{\link{GLES}}.
Variables 1 to 10 represent the response, variables 11 to 14 represent the subject-specific covariates. The response variables are ordinal, with values from 1 to 5. Low values represent string preference of the first-names party, high values represent strong preference of the last-named party. 
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
\dontrun{
# load data set
data(GLESsmall)

# define response and covariate matrix
X <- scale(GLESsmall[, 11:14])
Y <- as.matrix(GLESsmall[, 1:10])

# vector of subtitles, containing the coding of the single covariates
subs <- c("(in years)","female (1); male (0)",
"East Germany (1); West Germany (0)","(very) good (1); else (0)")

# vector of tuning parameters
lambda <- exp(seq(log(31),log(1),length=50))-1

# compute BTLLasso model
m <- BTLLasso(Y = Y, X = X, lambda = lambda)

op <- par(no.readonly = TRUE) 
par(mar=c(5,4,4,8))

# plot covariate paths
paths(m)

# plot parameter paths
singlepaths(m, subs = subs)

# compute 10-fold cross-validation
set.seed(5)
m.cv <- cv.BTLLasso(Y = Y, X = X, folds = 10, lambda = lambda, cores = 10)

# plot covariate paths, together with cv-optimal model
paths(m.cv)

# plot parameter paths, together with cv-optimal model
singlepaths(m.cv, subs = subs)

# compute bootstrap confidence intervals
m.boot <- boot.BTLLasso(m.cv, B = 100, cores = 25)

# plot bootstrap confidence intervals
par(mar=c(5,5,4,3))
ci.BTLLasso(m.boot, subs = subs)

par(op)
}
}
\keyword{datasets}
