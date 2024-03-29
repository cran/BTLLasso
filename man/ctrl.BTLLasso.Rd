% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BTLLasso.R
\name{ctrl.BTLLasso}
\alias{ctrl.BTLLasso}
\title{Control function for BTLLasso}
\usage{
ctrl.BTLLasso(
  l.lambda = 30,
  log.lambda = TRUE,
  lambda.min = 0.05,
  adaptive = TRUE,
  scale = TRUE,
  norm = c("L1", "L2"),
  epsilon = 1e-04,
  lambda2 = 1e-04,
  c = 1e-09,
  precision = 3,
  weight.penalties = TRUE,
  include.intercepts = TRUE,
  order.effect = FALSE,
  object.order.effect = FALSE,
  order.center = FALSE,
  name.order = "Order",
  penalize.intercepts = FALSE,
  penalize.X = TRUE,
  penalize.Z2 = FALSE,
  penalize.Z1.absolute = TRUE,
  penalize.Z1.diffs = TRUE,
  penalize.order.effect.absolute = TRUE,
  penalize.order.effect.diffs = FALSE
)
}
\arguments{
\item{l.lambda}{Number of tuning parameters. Applies only if \code{lambda = NULL} in the 
main function.}

\item{log.lambda}{Should the grid of tuning parameters be created on a logarithmic scale 
rather than equidistant. Applies only if \code{lambda = NULL} in the main function.}

\item{lambda.min}{Minimal value for tuning parameter. Applies only if \code{lambda = NULL} in the 
main function.}

\item{adaptive}{Should adaptive lasso be used? Default is TRUE.}

\item{scale}{Should the covariates be scaled so that they are on comparable scales? Default is TRUE.
Variables will be properly scaled AND centered. Please note that results will refer to scaled covariates.
If \code{adaptive = TRUE} scaling is not necessary to keep penalties comparable.}

\item{norm}{Specifies the norm used in the penalty term. Currently, only
'L1' and 'L2' are possible. Default is to 'L1', only 'L1' allows for
clustering and variable selection.}

\item{epsilon}{Threshold value for convergence of the algorithm.}

\item{lambda2}{Tuning parameter for ridge penalty on all coefficients.
Should be small, only used to stabilize results.}

\item{c}{Internal parameter for the quadratic approximation of the L1
penalty. Should be sufficiently small.}

\item{precision}{Precision for final parameter estimates, specifies number of decimals.}

\item{weight.penalties}{Should the penalties across the different model components 
(i.e. intercepts, order effects, X, Z1, Z2) be weighted according to the number of
penalties included? Default is \code{TRUE} to minimize the risk of selection bias
across different model components.}

\item{include.intercepts}{Should intercepts be included in the model?}

\item{order.effect}{Should a global order effect (corresponding to home effect
in sports applications) be included in the model?}

\item{object.order.effect}{Should object-specific order effects (corresponding to home effects
in sports applications) be included in the model?}

\item{order.center}{Should (in case of object-specific order effects) the order effects be centered
in the design matrix? Centering is equivalent to the coding scheme of effect coding instead of 
dummy coding.}

\item{name.order}{How should the order effect(s) be called in plots or prints.}

\item{penalize.intercepts}{Should intercepts be penalized? If \code{TRUE},
all pairwise differences between intercepts are penalized.}

\item{penalize.X}{Should effects from X matrix be penalized? If \code{TRUE},
all pairwise differences corresponding to one covariate are penalized. Can also be used with
a character vector as input. Then, the character vector contains the names of the variables
from X whose parameters should be penalized.}

\item{penalize.Z2}{Should absolute values of effects from Z2 matrix be
penalized? Can also be used with
a character vector as input. Then, the character vector contains the names of the variables
from Z2 whose parameters should be penalized.}

\item{penalize.Z1.absolute}{Should absolute values of effects from Z1 matrix
be penalized? Can also be used with
a character vector as input. Then, the character vector contains the names of the variables
from Z1 whose parameters should be penalized.}

\item{penalize.Z1.diffs}{Should differences of effects from Z1 matrix be
penalized? If \code{TRUE}, all pairwise differences corresponding to one
covariate are penalized. Can also be used with
a character vector as input. Then, the character vector contains the names of the variables
from Z1 whose parameters should be penalized.}

\item{penalize.order.effect.absolute}{Should absolute values of order effect(s) be penalized?
Only relevant if either \code{object.order.effect = TRUE} or \code{order.effect = TRUE}.}

\item{penalize.order.effect.diffs}{Should differences of order effects be
penalized? If \code{TRUE}, all pairwise differences are penalized. Only relevant if 
\code{object.order.effect = TRUE}}
}
\description{
Control parameters for different penalty terms and for tuning the fitting algorithm.
}
\examples{

\dontrun{
op <- par(no.readonly = TRUE)

##############################
##### Example with simulated data set containing X, Z1 and Z2
##############################
data(SimData)

## Specify control argument
## -> allow for object-specific order effects and penalize intercepts
ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
                      penalize.order.effect.diffs = TRUE)

## Simple BTLLasso model for tuning parameters lambda
m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
                  Z2 = SimData$Z2, control = ctrl)
m.sim

par(xpd = TRUE)
plot(m.sim)


## Cross-validate BTLLasso model for tuning parameters lambda
set.seed(1860)
m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
                        Z2 = SimData$Z2, control = ctrl)
m.sim.cv
coef(m.sim.cv)
logLik(m.sim.cv)

head(predict(m.sim.cv, type="response"))
head(predict(m.sim.cv, type="trait"))

plot(m.sim.cv, plots_per_page = 4)


## Example for bootstrap intervals for illustration only
## Don't calculate bootstrap intervals with B = 20!!!!
set.seed(1860)
m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 20, cores = 20)
m.sim.boot
plot(m.sim.boot, plots_per_page = 4)


##############################
##### Example with small version from GLES data set
##############################
data(GLESsmall)

## extract data and center covariates for better interpretability
Y <- GLESsmall$Y
X <- scale(GLESsmall$X, scale = FALSE)
Z1 <- scale(GLESsmall$Z1, scale = FALSE)

## vector of subtitles, containing the coding of the X covariates
subs.X <- c('', 'female (1); male (0)')

## Cross-validate BTLLasso model
m.gles.cv <- cv.BTLLasso(Y = Y, X = X, Z1 = Z1)
m.gles.cv

coef(m.gles.cv)
logLik(m.gles.cv)

head(predict(m.gles.cv, type="response"))
head(predict(m.gles.cv, type="trait"))

par(xpd = TRUE, mar = c(5,4,4,6))
plot(m.gles.cv, subs.X = subs.X, plots_per_page = 4, which = 2:5)
paths(m.gles.cv, y.axis = 'L2')


##############################
##### Example with Bundesliga data set
##############################
data(Buli1516)

Y <- Buli1516$Y5

Z1 <- scale(Buli1516$Z1, scale = FALSE)

ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
                           name.order = "Home", 
                           penalize.order.effect.diffs = TRUE, 
                           penalize.order.effect.absolute = FALSE,
                           order.center = TRUE, lambda2 = 1e-2)

set.seed(1860)
m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
m.buli

par(xpd = TRUE, mar = c(5,4,4,6))
plot(m.buli)


##############################
##### Example with Topmodel data set
##############################
data("Topmodel2007", package = "psychotree")

Y.models <- response.BTLLasso(Topmodel2007$preference)
X.models <- scale(model.matrix(preference~., data = Topmodel2007)[,-1])
rownames(X.models) <- paste0("Subject",1:nrow(X.models))
colnames(X.models) <- c("Gender","Age","KnowShow","WatchShow","WatchFinal")

set.seed(5)
m.models <- cv.BTLLasso(Y = Y.models, X = X.models)
plot(m.models, plots_per_page = 6)

par(op)
}
}
\references{
Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
Statistical Software}, 88(9), 1-29, \doi{10.18637/jss.v088.i09}

Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
17(3), 223 - 243

Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
Analysis of the importance of on-field covariates in the German Bundesliga, 
\emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
}
\seealso{
\code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
}
\author{
Gunther Schauberger\cr \email{gunther.schauberger@tum.de}
}
\keyword{BTLLasso}
\keyword{control}
