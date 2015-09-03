\name{singlepaths}
\alias{singlepaths}

\title{
Plot parameter paths for BTLLasso
}
\description{
Plots single paths for every parameter of a BTLLasso object or a cv.BTLLasso object. In contrast, to \code{\link{paths}}, one plot per covariate is created, every single parameter is illustrated by one path. For cv.BTLLasso object, the optimal model according to the cross-validation is marked by a vertical dashed line. 
}
\usage{
singlepaths(model, colors = NULL, equal.ranges = FALSE, subs = NULL)
}

\arguments{
  \item{model}{
BTLLasso or cv.BTLLasso object
}
  \item{colors}{
Optional. If specified, vector with length equal to the number of objects. Each object can be represented by another color. 
}
  \item{equal.ranges}{
Should all single plots (for different covariates) have equal ranges on the y-axes. FALSE by default. 
}
  \item{subs}{
Optional vector of subtitles for the single plots. Can be used to note the encoding of the single covariates, especially for dummy variables. 
}
}
\details{
Plots for BTLLasso and cv.BTLLasso objects only differ by the additional vertical line indicating the optimal model according to cross-validation.
}

\author{
Gunther Schauberger\cr
\email{gunther@stat.uni-muenchen.de}\cr
\url{http://www.statistik.lmu.de/~schauberger/}
}

\references{
Schauberger, Gunther and Tutz, Gerhard (2015): Modelling Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an Application to Party Preference Data, \emph{Department of Statistics, LMU Munich}, Technical Report 183
}

\seealso{
\code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}, \code{\link{paths}}
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

# plot parameter paths
singlepaths(m, subs = subs)

# compute 10-fold cross-validation
set.seed(5)
m.cv <- cv.BTLLasso(Y = Y, X = X, folds = 10, lambda = lambda, cores = 10)

# plot parameter paths, together with cv-optimal model
singlepaths(m.cv, subs = subs)

par(op)
}
}

\keyword{BTLLasso}
\keyword{paths}
\keyword{parameter paths}