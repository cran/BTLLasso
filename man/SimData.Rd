% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BTLLasso-package.R
\docType{data}
\name{SimData}
\alias{SimData}
\title{Simulated data set for illustration}
\format{
A list containing  simulated data for 200 observations. 
The list contains both information on the response (paired comparisons) and different covariates.
\describe{ 
\item{Y}{A response.BTLLasso object with simulated responses including
\itemize{
\item{response: Ordinal paired comparison response vector} 
\item{first.object: Vector containing the first-named object per paired comparison}
\item{second.object: Vector containing the second-named object per paired comparison}
\item{subject: Vector containing a subject identifier per paired comparison}
\item{with.order} Automatically generated vector containing information on order effect. Each paired 
comparison is associated with an order effect.
}}
\item{X}{Matrix containing both subject-specific covariates
\itemize{
\item{X_var1} 
\item{X_var2}
}
}
\item{Z1}{Matrix containing both subject-object-specific covariates
\itemize{
\item{Z1_var1}
\item{Z1_var2}
}
}
\item{Z2}{Matrix containing both object-specific covariates
\itemize{
\item{Z2_var1}
\item{Z2_var2}
}
}
}
}
\description{
This data set is a simulated data set including all possible types of covariates (X, Z1 and Z2)
and is intended to serve for illustration purpose. The data set contains paired comparisons between
four objects with five different response categories from 200 subjects.
}
\examples{

\dontrun{
op <- par(no.readonly = TRUE)

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

par(op)
}
}
\keyword{datasets}
