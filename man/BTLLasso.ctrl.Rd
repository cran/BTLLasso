\name{BTLLasso.ctrl}
\alias{BTLLasso.ctrl}

\title{
Control function for BTLLasso
}
\description{
Contains some additional parameters, mostly for internal use.
}
\usage{
BTLLasso.ctrl(adaptive = TRUE, norm = c("L1", "L2"), epsilon = 1e-04, lambda2 = 1e-04, 
              c = 1e-09, penal.diffs = TRUE,return.design = TRUE)
}

\arguments{
  \item{adaptive}{
Should adaptive lasso be used? Default is TRUE.
}
  \item{norm}{
Specifies the norm used in the penalty term. Currently, only "L1" and "L2" are possible. Default is to "L1", only "L1" allows for clustering and variable selection.
}
  \item{epsilon}{
Threshold value for convergence of the algorithm. 
}
  \item{lambda2}{
Tuning parameter for ridge penalty on all coefficients. Should be small, only used to stabilize results. 
}
  \item{c}{
Internal parameter for the quadratic approximation of the L1 penalty. Should be sufficiently small. For details see \code{\link[gvcm.cat]{cat_control}}.
}
  \item{penal.diffs}{
Should also differences between all parameters be penalized? Default is TRUE. If FALSE, only absolute values are penalized. 
}
  \item{return.design}{
Should the design matrix be returned from the BTLLasso function. Necessary for internal use. 
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

# compute BTLLasso model, increase accuracy of the algorithm
m <- BTLLasso(Y = Y, X = X, lambda = lambda, control = BTLLasso.ctrl(epsilon = 1e-05))

}
}

\keyword{BTLLasso}
\keyword{control}
