% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.BTLLasso.R
\name{cv.BTLLasso}
\alias{cv.BTLLasso}
\title{Cross-validation function for BTLLasso}
\usage{
cv.BTLLasso(
  Y,
  X = NULL,
  Z1 = NULL,
  Z2 = NULL,
  folds = 10,
  lambda = NULL,
  control = ctrl.BTLLasso(),
  cores = folds,
  trace = TRUE,
  trace.cv = TRUE,
  cv.crit = c("RPS", "Deviance")
)
}
\arguments{
\item{Y}{A \code{response.BTLLasso} object created by
\code{\link{response.BTLLasso}}.}

\item{X}{Matrix containing all \bold{subject-specific covariates} that are
to be included with \bold{object-specific effects}. One row represents one
subject, one column represents one covariate. X has to be standardized.}

\item{Z1}{Matrix containing all \bold{object-subject-specific covariates}
that are to be included with \bold{object-specific effects}. One row
represents one subject, one column represents one combination between
covariate and object. Column names have to follow the scheme
'firstvar.object1',...,'firstvar.objectm',...,'lastvar.objectm'. The object
names 'object1',...,'objectm' have to be identical to the object names used
in the \code{response.BTLLasso} object \code{Y}. The variable names and the
object names have to be separated by '.'.  The rownames of the matrix',
Z.name, 'have to be equal to the subjects specified in the response object.
Z1 has to be standardized.}

\item{Z2}{Matrix containing all \bold{object-subject-specific covariates or
object-specific covariates} that are to be included with \bold{global
effects}. One row represents one subject, one column represents one
combination between covariate and object. Column names have to follow the
scheme 'firstvar.object1',...,'firstvar.objectm',...,'lastvar.objectm'. The
object names 'object1',...,'objectm' have to be identical to the object
names used in the \code{response.BTLLasso} object \code{Y}. The variable
names and the object names have to be separated by '.'.  The rownames of the
matrix', Z.name, 'have to be equal to the subjects specified in the response
object. Z2 has to be standardized.}

\item{folds}{Number of folds for the crossvalidation. Default is 10.}

\item{lambda}{Vector of tuning parameters. If \code{NULL}, automatically a grid
of tuning parameters is created.}

\item{control}{Function for control arguments, mostly for internal use. See
also \code{\link{ctrl.BTLLasso}}.}

\item{cores}{Number of cores used for (parallelized) cross-validation. By
default, equal to the number of folds.}

\item{trace}{Should the trace of the BTLLasso algorithm be printed?}

\item{trace.cv}{Should the trace fo the cross-validation be printed? If
parallelized, the trace is not working on Windows machines.}

\item{cv.crit}{Which criterion should be used to evaluate cross-validation. Choice is 
between Ranked probability score and deviance. Only \code{RPS} considers the ordinal
structure of the response.}
}
\value{
\item{coefs}{Matrix containing all (original) coefficients, one row
per tuning parameter, one column per coefficient.} 
\item{coefs.repar}{Matrix
containing all reparameterized (for symmetric side constraint) coefficients,
one row per tuning parameter, one column per coefficient.}
\item{logLik}{Vector of log-likelihoods, one value per tuning parameter.}
\item{design}{List containing design matrix and several additional information like, 
e.g., number and names of covariates.} 
\item{Y}{Response object.} 
\item{penalty}{List containing all penalty matrices and some further information on penalties} 
\item{response}{Vector containing 0-1 coded
response.} 
\item{X}{X matrix containing subject-specific covariates.} 
\item{Z1}{Z1 matrix containing subject-object-specific covariates.} 
\item{Z2}{Z2 matrix containing (subject)-object-specific covariates.} 
\item{lambda}{Vector of tuning parameters.} 
\item{control}{Control argument, specified by \code{\link{ctrl.BTLLasso}}.}
\item{criterion}{Vector containing values of the chosen cross-validation criterion, 
one value per tuning parameter.}
\item{folds}{Number of folds in cross validation.} 
\item{cv.crit}{Cross-validation criterion, either \code{RPS} or \code{Deviance}.}
\item{df}{Vector containing degrees of freedom for all models along the grid 
of tuning parameters.}
}
\description{
Performs cross-validation of BTLLasso, including the BTLLasso algorithm for
the whole data set.
}
\details{
Cross-validation can be performed parallel, default is 10-fold 
cross-validation on 10 cores. Output is a cv.BTLLasso object which can then be
used for bootstrap intervalls using \code{\link{boot.BTLLasso}}.
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
\code{\link{BTLLasso}}, \code{\link{boot.BTLLasso}}, \code{\link{ctrl.BTLLasso}},
\code{\link{plot.BTLLasso}}, \code{\link{paths}}, \code{\link{print.cv.BTLLasso}}, 
\code{\link{predict.BTLLasso}}, \code{\link{coef}}
}
\author{
Gunther Schauberger\cr \email{gunther.schauberger@tum.de}
}
\keyword{BTLLasso}
\keyword{cross}
\keyword{validation}
