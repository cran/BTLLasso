\name{paths}
\alias{paths}

\title{
Plot covariate paths for BTLLasso
}
\description{
Plots paths for every covariate of a BTLLasso object or a cv.BTLLasso object. In contrast, to \code{\link{singlepaths}}, only one plot is created, every covariate is illustrated by one path. For cv.BTLLasso object, the optimal model according to the cross-validation is marked by a vertical dashed line.
}
\usage{
paths(model)
}

\arguments{
  \item{model}{
BTLLasso or cv.BTLLasso object
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
\code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}, \code{\link{singlepaths}}
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

# compute 10-fold cross-validation
set.seed(5)
m.cv <- cv.BTLLasso(Y = Y, X = X, folds = 10, lambda = lambda, cores = 10)

# plot covariate paths, together with cv-optimal model
paths(m.cv)

par(op)
}
}

\keyword{BTLLasso}
\keyword{paths}
\keyword{covariate paths}