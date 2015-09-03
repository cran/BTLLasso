\name{ci.BTLLasso}
\alias{ci.BTLLasso}

\title{
Plot confidence intervals for BTLLasso
}
\description{
Plots confidence intervals for every single coefficient. Confidence intervals are separated by covariates, every covariate is plotted separately. Confidence intervals are based on bootstrap, performed by \code{\link{boot.BTLLasso}}.
}
\usage{
ci.BTLLasso(object, subs = NULL, rows = NULL, subset = NULL)
}

\arguments{
  \item{object}{
boot.BTLLasso object
}
  \item{subs}{
Optional vector of subtitles for the single plots. Can be used to note the encoding of the single covariates, especially for dummy variables. 
}
  \item{rows}{
Optional argument for the number of rows, 
}
  \item{subset}{
Optional argument if only a subset of the covariates should be plotted. Vector specifying the numbers of the covariates that should be plotted.
}
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
\code{\link{boot.BTLLasso}}, \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
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
op <- par(no.readonly = TRUE) 
par(mar=c(5,5,4,3))
ci.BTLLasso(m.boot, subs = subs)

par(op)
}
}
\keyword{BTLLasso}
\keyword{confidence interval}
\keyword{bootstrap}
