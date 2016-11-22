#' Control function for BTLLasso
#' 
#' Control parameters for different penalty terms and for tuning the fitting algorithm.
#' 
#' 
#' @param adaptive Should adaptive lasso be used? Default is TRUE.
#' @param scale Should the covariates be scaled so that they are on comparable scales? Default is TRUE.
#' Variables will be properly scaled AND centered. Please note that results will refer to scaled covariates.
#' If \code{adaptive = TRUE} scaling is not necessary to keep penalties comparable.
#' @param norm Specifies the norm used in the penalty term. Currently, only
#' "L1" and "L2" are possible. Default is to "L1", only "L1" allows for
#' clustering and variable selection.
#' @param epsilon Threshold value for convergence of the algorithm.
#' @param lambda2 Tuning parameter for ridge penalty on all coefficients.
#' Should be small, only used to stabilize results.
#' @param c Internal parameter for the quadratic approximation of the L1
#' penalty. Should be sufficiently small. For details see
#' \code{\link[gvcm.cat]{cat_control}}.
#' @param precision Precision for final parameter estimates, specifies number of decimals.
#' @param weight.penalties Should the penalties across the different model components 
#' (i.e. intercepts, order effects, X, Z1, Z2) be weighted according to the number of
#' penalties included? Default is \code{TRUE} to minimize the risk of selection bias
#' across different model components.  
#' @param include.intercepts Should intercepts be included in the model?
#' @param order.effect Should a global order effect (corresponding to home effect
#' in sports applications) be included in the model?
#' @param object.order.effect Should object-specific order effects (corresponding to home effects
#' in sports applications) be included in the model?
#' @param order.center Should (in case of object-specific order effects) the order effects be centered
#' in the design matrix? Centering is equivalent to the coding scheme of effect coding instead of 
#' dummy coding.
#' @param name.order How should the order effect(s) be called in plots or prints.
#' @param penalize.intercepts Should intercepts be penalized? If \code{TRUE},
#' all pairwise differences between intercepts are penalized.
#' @param penalize.X Should effects from X matrix be penalized? If \code{TRUE},
#' all pairwise differences corresponding to one covariate are penalized.
#' @param penalize.Z2 Should absolute values of effects from Z2 matrix be
#' penalized?
#' @param penalize.Z1.absolute Should absolute values of effects from Z1 matrix
#' be penalized?
#' @param penalize.Z1.diffs Should differences of effects from Z1 matrix be
#' penalized? If \code{TRUE}, all pairwise differences corresponding to one
#' covariate are penalized.
#' @param penalize.order.effect.absolute Should absolute values of order effect(s) be penalized?
#' Only relevant if either \code{object.order.effect = TRUE} or \code{order.effect = TRUE}.
#' @param penalize.order.effect.diffs Should differences of order effects be
#' penalized? If \code{TRUE}, all pairwise differences are penalized. Only relevant if 
#' \code{object.order.effect = TRUE}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' @keywords BTLLasso control
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
#' @export ctrl.BTLLasso
ctrl.BTLLasso <- 
  function (adaptive = TRUE, scale = TRUE, norm = c("L1","L2"), epsilon = 1e-4, lambda2 = 1e-4, c = 1e-9,
            precision = 3,  weight.penalties = TRUE, include.intercepts = TRUE, 
            order.effect = FALSE, object.order.effect = FALSE, 
            order.center = FALSE, name.order = "Order", penalize.intercepts = FALSE, penalize.X = TRUE, 
            penalize.Z2 = FALSE,  penalize.Z1.absolute = TRUE, penalize.Z1.diffs = TRUE, 
            penalize.order.effect.absolute = TRUE, penalize.order.effect.diffs = FALSE) 
  { 
    norm <- match.arg(norm)
    
    RET <- list(adaptive = adaptive, scale = scale, norm = norm, epsilon = epsilon, lambda2 = lambda2,
                c = c, penalize.X = penalize.X, penalize.Z1.diffs = penalize.Z1.diffs, 
                penalize.Z2 = penalize.Z2, penalize.Z1.absolute = penalize.Z1.absolute,
                penalize.intercepts = penalize.intercepts, 
                include.intercepts = include.intercepts, order.effect = order.effect, 
                object.order.effect = object.order.effect, order.center = order.center,
                penalize.order.effect.diffs = penalize.order.effect.diffs, 
                penalize.order.effect.absolute = penalize.order.effect.absolute,
                name.order = name.order, precision = precision, 
                weight.penalties = weight.penalties)
    RET
  }




#' Function to perform BTLLasso
#' 
#' Performs BTLLasso, a method to model heterogeneity in paired comparison
#' data. Different types of covariates are allowd to have an influence on the
#' attractivity/strength of the objects. Covariates can be subject-specific, 
#' object-specific or subject-object-specific. L1 penalties are used to reduce the 
#' complexiy of the model by enforcing clusters of equal effects or by elimination of irrelevant
#' covariates.  
#' 
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
#' @param lambda Vector of tuning parameters.
#' @param control Function for control arguments, mostly for internal use. See
#' also \code{\link{ctrl.BTLLasso}}.
#' @param trace Should the trace of the BTLLasso algorithm be printed?
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
#' \item{penalty}{List containing all penalty matrices and some further information on penalties.} 
#' \item{response}{Vector containing 0-1 coded
#' response.} 
#' \item{X}{X matrix containing subject-specific covariates.} 
#' \item{Z1}{Z1 matrix containing subject-object-specific covariates.} 
#' \item{Z2}{Z2 matrix containing (subject)-object-specific covariates.} 
#' \item{lambda}{Vector of tuning parameters.} 
#' \item{control}{Control argument, specified by \code{\link{ctrl.BTLLasso}}.}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{cv.BTLLasso}}, \code{\link{boot.BTLLasso}}, \code{\link{ctrl.BTLLasso}},
#' \code{\link{singlepaths}}, \code{\link{paths}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' @keywords BTLLasso
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
#' @export BTLLasso
BTLLasso <- function(Y, X=NULL, Z1=NULL, Z2=NULL, 
                     lambda, control = ctrl.BTLLasso(), trace = TRUE){


  ## create design matrix
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

  ## create penalty matrix
  get.penalties <- penalties.BTLLasso(Y = Y, X = X, Z1 = Z1, Z2 = Z2, control = control)
  
  ## fit BTLLasso model, with coefficients and degrees of freedom
  fit <- fit.BTLLasso(response, get.design$design, get.penalties, lambda, Y$k, Y$m, control, trace)
  coefs <- fit$coefs
  df <- fit$df
  
  ## calculate log likelihood
  logLik <- c()
  for(j in 1:nrow(coefs)){
    logLik[j] <- loglik(coefs[j,],Y$response,get.design$design,Y$k)
  }

  ## reparameterize coefficients, from reference object to symmetric side constraint
  coefs.repar <- round(expand.coefs(coefs, get.design, Y),control$precision)
  
  ## return stuff
  ret.list <- list(coefs = coefs, coefs.repar = coefs.repar, logLik = logLik, design = get.design, Y = Y, 
                   penalty = get.penalties, response = response, X = X, Z1 = Z1, Z2 = Z2, lambda = lambda, 
                   control = control)
  
  class(ret.list) <- "BTLLasso"
  
  return(ret.list)
  
}
