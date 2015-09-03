\name{cv.BTLLasso}
\alias{cv.BTLLasso}

\title{
Cross-validation function for BTLLasso
}
\description{
Performs crossvalidation of BTLLasso, including the BTLLasso algorithm for the whole data set. 
}
\usage{
cv.BTLLasso(Y, X, folds = 10, lambda, control = BTLLasso.ctrl(), 
            cores = folds, trace = TRUE, trace.cv = TRUE)
}

\arguments{
  \item{Y}{
Matrix containing paired comparisons, ordered or binary. For K categories, the colmuns have to be ordered as follows: 
1 vs. 2, 1 vs.3, ..., 1 vs. K, 2 vs. 3, 2 vs. 4, ..., 2 vs. K, ..., K-1 vs. K. One row represents one subject.  
}
  \item{X}{
Matrix containing all subject-specific covariates. One row represents one subject, one column represents one covariate. Has to be standardized. 
}
  \item{folds}{
Number of folds for the crossvalidation. Default is 10. 
}
  \item{lambda}{
Vector of tuning parameters.
}
  \item{control}{
Function for control arguments, mostly for internal use. See also \code{\link{BTLLasso.ctrl}}.
}
  \item{cores}{
Number of cores used for (parallelized) cross-validation. By default, equal to the number of folds.
}
  \item{trace}{
Should the trace of the BTLLasso algorithm be printed?
}
  \item{trace.cv}{
Should the trace fo the cross-validation be printed? If parallelized, the trace is not working on Windows machines. 
}
}
\details{
Cross validation can be performed parallel, default is 10-fold cross validation on 10 cores. Output is a cv.BTLLasso object which can then be used for bootstrap confidence intervalls using \code{\link{boot.BTLLasso}} and \code{\link{ci.BTLLasso}}.
}

  
\value{
\item{coefs}{Matrix containing all (original) coefficients, one row per tuning parameter, one column per coefficient.}
\item{coefs.repar}{Matrix containing all reparameterized (for symmetric side constraint) coefficients, one row per tuning parameter, one column per coefficient.}
\item{logLik}{Vector of log-likelihoods, one per tuning parameter.}
\item{design}{Design matrix, NULL if \code{return.design=FALSE} in \code{\link{BTLLasso.ctrl}}}
\item{Y}{Y matrix}
\item{q}{q=K-1, K is the number of different categories in \code{Y}}
\item{acoefs}{Matrix containing penalties, for internal use.}
\item{response}{Vector containing 0-1 coded response.}
\item{n}{Number of persons/subjects}
\item{I}{Number of paired comparisons}
\item{m}{Number of objects}
\item{p}{Number of covariates}
\item{X}{X matrix}
\item{n.theta}{Number of estimated threshold parameters}
\item{lambda}{Vector of tuning parameters}
\item{deviances}{Vector of cross-validation deviances, one value per tuning parameter.}
\item{folds}{Number of folds in cross validation.}
\item{labels}{Labels of objects, only correct if Y specified correctly by "1 vs. 2","1 vs. 3",...}
\item{epsilon}{Threshold value for convergence of the algorithm, specified in \code{\link{BTLLasso.ctrl}}}
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
\code{\link{BTLLasso}}, \code{\link{boot.BTLLasso}}, \code{\link{singlepaths}}, \code{\link{paths}}
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
\keyword{BTLLasso}
\keyword{cross validation}
