\name{BTLLasso}
\alias{BTLLasso}

\title{
Function to perform BTLLasso
}
\description{
Performs BTLLasso, a method to model heterogeneity in paired comparison data. Subject-specific covariates are allowd to have an influence on the attractivity/strength of the objects. An L1 penalty on the pairwise differences between the object-specific parameters allows for both clustering of object with regard to covariates and elimination of irrelevant covariates. 
}
\usage{
BTLLasso(Y, X, lambda, control = BTLLasso.ctrl(), trace = TRUE)
}

\arguments{
  \item{Y}{
Matrix containing paired comparisons, ordered or binary. For K categories, the colmuns have to be ordered as follows: 
1 vs 2, 1 vs 3, ..., 1 vs K, 2 vs 3, 2 vs 4, ..., 2 vs K, ..., K-1 vs K. One row represents one subject.  
}
  \item{X}{
Matrix containing all subject-specific covariates. One row represents one subject, one column represents one covariate. Has to be standardized. 
}
  \item{lambda}{
Vector of tuning parameters.
}
  \item{control}{
Function for control arguments, mostly for internal use. See also \code{\link{BTLLasso.ctrl}}.
}
  \item{trace}{
Should the trace of the BTLLasso algorithm be printed?
}
}
\details{
Y matrix needs a special format. Column names need to be in the style "1 vs 2","1 vs 3",... If so, object labels "1", "2", ... are extracted from the column names. 
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
\item{labels}{Labels of objects, only correct if Y specified correctly by "1 vs 2","1 vs 3",...}
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
\code{\link{cv.BTLLasso}}, \code{\link{boot.BTLLasso}}, \code{\link{singlepaths}}, \code{\link{paths}}
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

par(op)
}
}

\keyword{BTLLasso}

