#' Print function for BTLLasso objects
#' 
#' Prints the most important output of BTLLasso objects.
#' 
#' @method print BTLLasso
#' @param x \code{BTLLasso} object
#' @param \dots possible further arguments for print command
#' @return \item{x}{BTLLasso object}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2016): Modelling 
#' Football Results in the German Bundesliga Using Match-specific Covariates, 
#' \emph{Department of Statistics, LMU Munich}, Technical Report 197
#' @keywords BTLLasso
#' @examples
#'
#' \dontrun{
#' ##### Example with simulated data set containing X, Z1 and Z2
#' data(SimData)
#' 
#' ## Specify tuning parameters
#' lambda <- exp(seq(log(151), log(1.05), length = 30)) - 1
#' 
#' ## Specify control argument
#' ## -> allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1, 
#'                   Z2 = SimData$Z2, lambda = lambda, control = ctrl)
#' print(m.sim)
#' 
#' singlepaths(m.sim)
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(5)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1, 
#'                         Z2 = SimData$Z2, lambda = lambda, control = ctrl)
#' print(m.sim.cv)
#' 
#' singlepaths(m.sim.cv, plot.order.effect = FALSE, 
#'             plot.intercepts = FALSE, plot.Z2 = FALSE)
#' paths(m.sim.cv, y.axis = 'L2')
#' 
#' ## Example for bootstrap confidence intervals for illustration only
#' ## Don't calculate bootstrap confidence intervals with B = 10!!!!
#' set.seed(5)
#' m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 10, cores = 10)
#' print(m.sim.boot)
#' ci.BTLLasso(m.sim.boot)
#' 
#' 
#' ##### Example with small version from GLES data set
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
#' ## vector of tuning parameters
#' lambda <- exp(seq(log(61), log(1), length = 30)) - 1
#' 
#' 
#' ## compute BTLLasso model 
#' m.gles <- BTLLasso(Y = Y, X = X, Z1 = Z1, lambda = lambda)
#' print(m.gles)
#' 
#' singlepaths(m.gles, subs.X = subs.X)
#' paths(m.gles, y.axis = 'L2')
#' 
#' ## Cross-validate BTLLasso model 
#' m.gles.cv <- cv.BTLLasso(Y = Y, X = X, Z1 = Z1, lambda = lambda)
#' print(m.gles.cv)
#' 
#' singlepaths(m.gles.cv, subs.X = subs.X)
#' }
print.BTLLasso <- function(x, ...) {
  
  m <- x$Y$m
  n <- x$Y$n
  k <- x$Y$q + 1
  n.theta <- x$design$n.theta
  n.intercepts <- x$design$n.intercepts
  if (n.intercepts != 0) {
    n.intercepts <- n.intercepts + 1
  }
  n.order <- x$design$n.order
  p.X <- x$design$p.X
  p.Z1 <- x$design$p.Z1
  p.Z2 <- x$design$p.Z2
  lambda <- x$lambda
  
  vars.X <- x$design$vars.X
  vars.Z1 <- x$design$vars.Z1
  vars.Z2 <- x$design$vars.Z2
  
  labels <- x$Y$object.names
  
  cat("Output of BTLLasso estimation:", "\n")
  
  cat("---", "\n")
  
  cat("Setting:")
  cat("\n", n, "subjects")
  cat("\n", m, "objects")
  cat("\n", k, "response categories")
  cat("\n", p.X, "subject-specific covariate(s)")
  cat("\n", p.Z1, "subject-object-specific covariate(s) with object-specific effects")
  cat("\n", p.Z2, "(subject-)object-specific covariate(s) with global effects")
  if (n.order == m) {
    cat("\n", n.order, "subject-specific order effects")
  }
  if (n.order == 1) {
    cat("\n", "Global order effect")
  }
  if (n.order == 0) {
    cat("\n", "No order effect")
  }
  cat("\n", length(lambda), "different tuning parameters", 
    "\n")
  
  
  invisible(x)
  
  
}

