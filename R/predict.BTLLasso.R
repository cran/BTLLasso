#' Predict function for BTLLasso
#' 
#' Predict function for a \code{BTLLasso} object or a \code{cv.BTLLasso}
#' object. Predictions can be linear predictors, probabilities or the values 
#' of the latent traits for both competitors in the paired comparisons. 
#' 
#' Results are lists of matrices with prediction for every single tuning parameter 
#' for \code{BTLLasso} objects
#' and a single matrix for \code{cv.BTLLasso} objects.
#' 
#' @param object \code{BTLLasso} or \code{cv.BTLLasso} object
#' @param newdata List possibly containing slots Y, X, Z1 and Z2 to use new data for prediction.
#' @param type Type "link" gives the linear predictors for separate categories, 
#' type "response" gives the respective probabilities. Type "trait" gives the estimated latent traits
#' of both competitors/objects in the paired comparisons. 
#' @param ... Further predict arguments.
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, 88(9), 1-29, \doi{10.18637/jss.v088.i09}
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @keywords BTLLasso paths parameter paths
#' @examples
#' 
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' ##############################
#' ##### Example with simulated data set containing X, Z1 and Z2
#' ##############################
#' data(SimData)
#' 
#' ## Specify control argument
#' ## -> allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                   Z2 = SimData$Z2, control = ctrl)
#' m.sim
#' 
#' par(xpd = TRUE)
#' plot(m.sim)
#' 
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(1860)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                         Z2 = SimData$Z2, control = ctrl)
#' m.sim.cv
#' coef(m.sim.cv)
#' logLik(m.sim.cv)
#' 
#' head(predict(m.sim.cv, type="response"))
#' head(predict(m.sim.cv, type="trait"))
#' 
#' plot(m.sim.cv, plots_per_page = 4)
#' 
#' 
#' ## Example for bootstrap intervals for illustration only
#' ## Don't calculate bootstrap intervals with B = 20!!!!
#' set.seed(1860)
#' m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 20, cores = 20)
#' m.sim.boot
#' plot(m.sim.boot, plots_per_page = 4)
#' 
#' 
#' ##############################
#' ##### Example with small version from GLES data set
#' ##############################
#' data(GLESsmall)
#' 
#' ## extract data and center covariates for better interpretability
#' Y <- GLESsmall$Y
#' X <- scale(GLESsmall$X, scale = FALSE)
#' Z1 <- scale(GLESsmall$Z1, scale = FALSE)
#' 
#' ## vector of subtitles, containing the coding of the X covariates
#' subs.X <- c('', 'female (1); male (0)')
#' 
#' ## Cross-validate BTLLasso model
#' m.gles.cv <- cv.BTLLasso(Y = Y, X = X, Z1 = Z1)
#' m.gles.cv
#' 
#' coef(m.gles.cv)
#' logLik(m.gles.cv)
#' 
#' head(predict(m.gles.cv, type="response"))
#' head(predict(m.gles.cv, type="trait"))
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.gles.cv, subs.X = subs.X, plots_per_page = 4, which = 2:5)
#' paths(m.gles.cv, y.axis = 'L2')
#' 
#' 
#' ##############################
#' ##### Example with Bundesliga data set
#' ##############################
#' data(Buli1516)
#' 
#' Y <- Buli1516$Y5
#' 
#' Z1 <- scale(Buli1516$Z1, scale = FALSE)
#' 
#' ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
#'                            name.order = "Home", 
#'                            penalize.order.effect.diffs = TRUE, 
#'                            penalize.order.effect.absolute = FALSE,
#'                            order.center = TRUE, lambda2 = 1e-2)
#' 
#' set.seed(1860)
#' m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
#' m.buli
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.buli)
#' 
#' 
#' ##############################
#' ##### Example with Topmodel data set
#' ##############################
#' data("Topmodel2007", package = "psychotree")
#' 
#' Y.models <- response.BTLLasso(Topmodel2007$preference)
#' X.models <- scale(model.matrix(preference~., data = Topmodel2007)[,-1])
#' rownames(X.models) <- paste0("Subject",1:nrow(X.models))
#' colnames(X.models) <- c("Gender","Age","KnowShow","WatchShow","WatchFinal")
#' 
#' set.seed(5)
#' m.models <- cv.BTLLasso(Y = Y.models, X = X.models)
#' plot(m.models, plots_per_page = 6)
#' 
#' par(op)
#' }
predict.BTLLasso <- function(object, newdata = list(), 
                                type = c("link", "response", "trait"), ...){
  
  type <- match.arg(type)
  
  if(inherits(object,"cv.BTLLasso")){
    coef.cv <- object$coefs.repar[which.min(object$criterion),,drop=FALSE]
    
    ret.list <- predict.help(coef.cv, object = object, newdata = newdata, 
                             type = type)
    ret.list <- ret.list[[1]]
  }else{
    ret.list <- predict.help(object$coefs.repar, object = object, newdata = newdata, 
                             type = type)
  }

   
   ret.list
   
}