\name{print.cv.BTLLasso}
\alias{print.cv.BTLLasso}

\title{
Print function for cv.BTLLasso objects
}
\description{
Prints the most important output of cv.BTLLasso objects.
}
\usage{
\method{print}{cv.BTLLasso}(x, ...)
}

\arguments{
  \item{x}{
\code{cv.BTLLasso} object
}
  \item{\dots}{
possible further arguments for print command
}
}

\value{
\item{x}{cv.BTLLasso object}
}
\references{
Schauberger, Gunther and Tutz, Gerhard (2015): Modelling Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an Application to Party Preference Data, \emph{Department of Statistics, LMU Munich}, Technical Report 183
}
\author{
Gunther Schauberger\cr
\email{gunther@stat.uni-muenchen.de}\cr
\url{http://www.statistik.lmu.de/~schauberger/}
}

\seealso{
\code{\link{cv.BTLLasso}}
}
\examples{
\dontrun{
## load data set
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

print(m.cv)
}
}

\keyword{BTLLasso}

