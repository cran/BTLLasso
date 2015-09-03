\name{BTLLasso-package}
\alias{BTLLasso-package}

\docType{package}
\title{
BTLLasso
}
\description{
Performs BTLLasso, a method to model heterogeneity in paired comparison data. Subject-specific covariates are allowd to have an influence on the attractivity/strength of the objects. An L1 penalty on the pairwise differences between the object-specific parameters allows for both clustering of object with regard to covariates and elimination of irrelevant covariates. 
Several additional functions are provided, such as cross-validation, bootstraped confidence intervals, and several plot functions.  
}


\author{
Gunther Schauberger\cr
\email{gunther@stat.uni-muenchen.de}\cr
\url{http://www.statistik.lmu.de/~schauberger/}
}

\references{
Schauberger, Gunther and Tutz, Gerhard (2015): Modelling Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an Application to Party Preference Data, \emph{Department of Statistics, LMU Munich}, Technical Report 183
}

\keyword{ package }
\keyword{BTL}
\keyword{Bradley-Terry}
\keyword{BTLLasso}

\seealso{
\code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}} 
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

