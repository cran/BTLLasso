\name{boot.BTLLasso}
\alias{boot.BTLLasso}

\title{
Bootstrap function for BTLLasso
}
\description{
Performs bootstrap for BTLLasso to get bootstrap confidence intervals. Main input argument is a \code{cv.BTLLasso} object. The bootstrap is performed on level of the cross-validation. Therefore, within every bootstrap iteration the complete cross-validation procedure from the \code{cv.BTLLasso} object is performed.
The resulting \code{boot.BTLLasso} object is put into \code{\link{ci.BTLLasso}} to plot bootstrap confidence intervals.
}
\usage{
boot.BTLLasso(model, B = 500, lambda = NULL, cores = 1, trace = TRUE, trace.cv = TRUE,
              control = BTLLasso.ctrl())
}

\arguments{
  \item{model}{
A \code{cv.BTLLasso} object.
}
  \item{B}{
Number of bootstrap iterations.
}
  \item{lambda}{
Vector of tuning parameters. If not specified (default), tuning parameters from \code{cv.BTLLasso} object are used. See also details.
}
  \item{cores}{
Number of cores for (parallelized) computation.
}
  \item{trace}{
Should the trace of the BTLLasso algorithm be printed?
}
  \item{trace.cv}{
Should the trace fo the cross-validation be printed? If parallelized, the trace is not working on Windows machines.
}
  \item{control}{
Function for control arguments, mostly for internal use. See also \code{\link{BTLLasso.ctrl}}.
}
}
\details{
The method can be highly time-consuming, for high numbers of tuning parameters, high numbers of folds in the crossvalidation and high number of bootstrap iterations B. 
The number of tuning parameters can be reduced by specifying lambda in the \code{boot.BTLLasso} function. You can control if the range of prespecified tuning parameters was to small by looking at the output values \code{lambda.max.alert} and \code{lambda.min.alert}. They are set TRUE if the smallest or largest of the specifed lambda values was chosen in at least one bootstrap iteration. 

}

\value{
\item{cv.model}{\code{cv.BTLLasso} object}
\item{estimatesB}{Matrix containing all B estimates for original parameters. For internal use.}
\item{estimatesBrepar}{Matrix containing all B estimates for reparameterized (symmetric side constraints) parameters.}
\item{lambdaB}{vector of used tuning parameters}
\item{conf.ints}{Bootstrap confidence intervals for original parameters. For internal use.}
\item{conf.ints.repar}{Bootstrap confidence intervals for reparameterized (symmetric side constraint) parameters.}
\item{lambda.max.alert}{Was the largest value of lambda chosen in at least one bootstrap iteration?}
\item{lambda.min.alert}{Was the smallest value of lambda chosen in at least one bootstrap iteration?}
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
\code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}, \code{\link{ci.BTLLasso}}
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
