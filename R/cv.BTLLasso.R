#' Cross-validation function for BTLLasso
#' 
#' Performs crossvalidation of BTLLasso, including the BTLLasso algorithm for
#' the whole data set.
#' 
#' Cross validation can be performed parallel, default is 10-fold cross
#' validation on 10 cores. Output is a cv.BTLLasso object which can then be
#' used for bootstrap confidence intervalls using \code{\link{boot.BTLLasso}}
#' and \code{\link{ci.BTLLasso}}.
#' 
#' @param Y A \code{response.BTLLasso} object created by
#' \code{\link{response.BTLLasso}}.
#' @param X Matrix containing all \bold{subject-specific covariates} that are
#' to be included with \bold{object-specific effects}. One row represents one
#' subject, one column represents one covariate. X has to be standardized.
#' @param Z1 Matrix containing all \bold{object-subject-specific covariates}
#' that are to be included with \bold{object-specific effects}. One row
#' represents one subject, one column represents one combination between
#' covariate and object. Column names have to follow the scheme
#' 'firstvar.object1',...,'firstvar.objectm',...,'lastvar.objectm'. The object
#' names 'object1',...,'objectm' have to be identical to the object names used
#' in the \code{response.BTLLasso} object \code{Y}. The variable names and the
#' object names have to be separated by '.'.  The rownames of the matrix",
#' Z.name, "have to be equal to the subjects specified in the response object.
#' Z1 has to be standardized.
#' @param Z2 Matrix containing all \bold{object-subject-specific covariates or
#' object-specific covariates} that are to be included with \bold{global
#' effects}. One row represents one subject, one column represents one
#' combination between covariate and object. Column names have to follow the
#' scheme 'firstvar.object1',...,'firstvar.objectm',...,'lastvar.objectm'. The
#' object names 'object1',...,'objectm' have to be identical to the object
#' names used in the \code{response.BTLLasso} object \code{Y}. The variable
#' names and the object names have to be separated by '.'.  The rownames of the
#' matrix", Z.name, "have to be equal to the subjects specified in the response
#' object. Z2 has to be standardized.
#' @param folds Number of folds for the crossvalidation. Default is 10.
#' @param lambda Vector of tuning parameters.
#' @param control Function for control arguments, mostly for internal use. See
#' also \code{\link{ctrl.BTLLasso}}.
#' @param cores Number of cores used for (parallelized) cross-validation. By
#' default, equal to the number of folds.
#' @param trace Should the trace of the BTLLasso algorithm be printed?
#' @param trace.cv Should the trace fo the cross-validation be printed? If
#' parallelized, the trace is not working on Windows machines.
#' @param cv.crit Which criterion should be used to evaluate cross-validation. Choice is 
#' between Ranked probability score and deviance. Only \code{RPS} considers the ordinal
#' structure of the response.
#' @return 
#' \item{coefs}{Matrix containing all (original) coefficients, one row
#' per tuning parameter, one column per coefficient.} 
#' \item{coefs.repar}{Matrix
#' containing all reparameterized (for symmetric side constraint) coefficients,
#' one row per tuning parameter, one column per coefficient.}
#' \item{logLik}{Vector of log-likelihoods, one value per tuning parameter.}
#' \item{design}{List containing design matrix and several additional information like, 
#' e.g., number and names of covariates.} 
#' \item{Y}{Response object.} 
#' \item{penalty}{List containing all penalty matrices and some further information on penalties} 
#' \item{response}{Vector containing 0-1 coded
#' response.} 
#' \item{X}{X matrix containing subject-specific covariates.} 
#' \item{Z1}{Z1 matrix containing subject-object-specific covariates.} 
#' \item{Z2}{Z2 matrix containing (subject)-object-specific covariates.} 
#' \item{lambda}{Vector of tuning parameters.} 
#' \item{control}{Control argument, specified by \code{\link{ctrl.BTLLasso}}.}
#' \item{criterion}{Vector containing values of the chosen cross-validation criterion, 
#' one value per tuning parameter.}
#' \item{folds}{Number of folds in cross validation.} 
#' \item{cv.crit}{Cross-validation criterion, either \code{RPS} or \code{Deviance}.}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{BTLLasso}}, \code{\link{boot.BTLLasso}}, \code{\link{ctrl.BTLLasso}},
#' \code{\link{singlepaths}}, \code{\link{paths}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' @keywords BTLLasso cross validation
#' @examples
#' 
#' \dontrun{
#' ##### Example with simulated data set containing X, Z1 and Z2
#' data(SimData)
#' 
#' ## Specify tuning parameters
#' lambda <- exp(seq(log(151),log(1.05),length=30))-1
#' 
#' ## Specify control argument, allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1, 
#'                   Z2 = SimData$Z2, lambda = lambda, control = ctrl)
#' 
#' singlepaths(m.sim, x.axis = "loglambda")
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(5)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1, 
#'                         Z2 = SimData$Z2, lambda = lambda, control = ctrl)
#' 
#' 
#' singlepaths(m.sim.cv, x.axis = "loglambda", plot.order.effect = FALSE, plot.intercepts = FALSE, 
#'             plot.Z2 = FALSE)
#' paths(m.sim.cv, y.axis="L2")
#' 
#' ## Example for bootstrap confidence intervals for illustration only
#' ## Don't calculate bootstrap confidence intervals with B = 10
#' set.seed(5)
#' m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 10, cores = 10)
#' ci.BTLLasso(m.sim.boot)
#' 
#' ##### Example with small version from GLES data set
#' data(GLESsmall)
#' 
#' ## vector of subtitles, containing the coding of the X covariates
#' subs.X <- c("","female (1); male (0)")
#' 
#' ## vector of tuning parameters
#' lambda <- exp(seq(log(61),log(1),length=30))-1
#' 
#' 
#' ## compute BTLLasso model 
#' m.gles <- BTLLasso(Y = GLESsmall$Y, X = GLESsmall$X, Z1 = GLESsmall$Z1, lambda = lambda)
#' 
#' singlepaths(m.gles, x.axis = "loglambda", subs.X = subs.X)
#' paths(m.gles, y.axis="L2")
#' 
#' ## Cross-validate BTLLasso model 
#' m.gles.cv <- cv.BTLLasso(Y = GLESsmall$Y, X = GLESsmall$X, Z1 = GLESsmall$Z1, lambda = lambda)
#' 
#' singlepaths(m.gles.cv, x.axis = "loglambda", subs.X = subs.X)
#' }
#' 
#' @export cv.BTLLasso
cv.BTLLasso <- function(Y, X = NULL, Z1 = NULL, Z2 = NULL, folds = 10, lambda, control = ctrl.BTLLasso(), 
                        cores = folds, trace = TRUE, trace.cv = TRUE, cv.crit = c("RPS","Deviance")){

  cv.crit <- match.arg(cv.crit)


  get.design <- design.BTLLasso(Y = Y, X = X, Z1 = Z1, Z2 = Z2, control = control)
  
  ## exclude missing values
  na.response <- is.na(Y$response)
  na.design <- colSums(matrix(is.na(rowSums(get.design$design)),nrow=Y$q))!=0
  na.total <- na.response | na.design
  Y$response <- Y$response[!na.total]
  Y$first.object <- Y$first.object[!na.total]
  Y$second.object <- Y$second.object[!na.total]
  Y$subject <- Y$subject[!na.total]
  Y$subject.names <- levels(as.factor(Y$subject))
  Y$n <- length(Y$subject.names)
  
  get.design$design <- get.design$design[!rep(na.total,each=Y$q),]
  
  ## create response vector
  if(Y$q==1){
    response <- as.numeric(Y$response)-1
  }else{
    response <- cumul.response(Y)
  }
  
  get.penalties <- penalties.BTLLasso(Y = Y, X = X, Z1 = Z1, Z2 = Z2, control = control)
  
  
  fit <- fit.cv.BTLLasso(response = response, design = get.design$design, penalty = get.penalties, 
                      q = Y$q, m = Y$m, folds = folds, lambda = lambda, control = control, 
                      cores = cores, trace = trace, trace.cv = trace.cv, cv.crit = cv.crit)
  
  
  logLik <- c()
  for(j in 1:nrow(fit$coefs)){
    logLik[j] <- loglik(fit$coefs[j,], Y$response, get.design$design, Y$k)
  }
  
  
  coefs.repar <- round(expand.coefs(fit$coefs, get.design, Y), control$precision)


  ret.list <- list(coefs = fit$coefs, coefs.repar = coefs.repar, logLik = logLik, 
                   design = get.design, Y = Y, penalty = get.penalties,
                   response = response, X = X, Z1 = Z1, Z2 = Z2, lambda = lambda, 
                   control = control, criterion = fit$criterion, folds = folds,
                   cv.crit = cv.crit)
  
  
  class(ret.list) <- c("cv.BTLLasso", "BTLLasso")
  
  return(ret.list)
}
