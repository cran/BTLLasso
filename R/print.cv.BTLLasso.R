#' Print function for cv.BTLLasso objects
#' 
#' Prints the most important output of cv.BTLLasso objects.
#' 
#' @method print cv.BTLLasso
#' @param x \code{cv.BTLLasso} object
#' @param rescale Should the parameter estimates be rescaled for plotting? Only 
#' applies if \code{scale = TRUE} was specified in \code{BTLLasso} or \code{cv.BTLLasso}.
#' @param \dots possible further arguments for print command
#' @return \item{x}{cv.BTLLasso object}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{cv.BTLLasso}}
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
print.cv.BTLLasso <- function(x, rescale = FALSE, ...) {
  
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
  
  cv.crit <- x$cv.crit
  
  
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
  cat("\n Cross-validation criterion:", cv.crit, "\n")
  
  cat("---", "\n")
  
  cat("Parameter estimates after", x$folds, "-", "fold cross-validation", 
    "\n")
  
  cat("\n")
  coefs <- x$coefs.repar[which.min(x$criterion), ]
  
  theta <- intercepts <- order.effects <- gamma.X <- gamma.Z1 <- gamma.Z2 <- c()
  
  if (n.theta > 0) {
    cat("Thresholds:", "\n")
    theta <- coefs[1:n.theta]
    names(theta) <- paste0("theta", 1:n.theta)
    print(theta)
    cat("\n")
  }
  
  if (n.order > 0) {
    cat(paste0(x$control$name.order, ":"), "\n")
    orders <- coefs[(n.theta + 1):(n.theta + n.order)]
    if (n.order == m) {
      names(orders) <- labels
    }
    if (n.order == 1) {
      names(orders) <- NULL
    }
    print(orders)
    cat("\n")
  }
  
  if (n.intercepts > 0) {
    cat("Intercepts:", "\n")
    intercepts <- coefs[(n.theta + n.order + 1):(n.theta + 
      n.order + n.intercepts)]
    names(intercepts) <- labels
    print(intercepts)
    cat("\n")
  }
  
  if (p.X > 0) {
    cat("Object-specific effects for subject-specific covariate(s):", 
      "\n")
    gamma.X <- matrix(coefs[(n.theta + n.order + n.intercepts + 
      1):(n.theta + n.order + n.intercepts + p.X * m)], 
      nrow = p.X, byrow = TRUE)
    if (rescale) {
      gamma.X <- t(t(gamma.X)/rep(x$design$sd.X, each = m))
    }
    colnames(gamma.X) <- labels
    rownames(gamma.X) <- vars.X
    print(gamma.X)
    cat("\n")
  }
  
  if (p.Z1 > 0) {
    cat("Object-specific effects for subject-object-specific covariate(s):", 
      "\n")
    gamma.Z1 <- matrix(coefs[(n.theta + n.order + n.intercepts + 
      p.X * m + 1):(n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m)], nrow = p.Z1, byrow = TRUE)
    if (rescale) {
      gamma.Z1 <- t(t(gamma.Z1)/rep(x$design$sd.Z1, each = m))
    }
    colnames(gamma.Z1) <- labels
    rownames(gamma.Z1) <- vars.Z1
    print(gamma.Z1)
    cat("\n")
  }
  
  if (p.Z2 > 0) {
    
    cat("Global effects for (subject-)object-specific covariate(s):", 
      "\n")
    gamma.Z2 <- coefs[(n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m + 1):(n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m + p.Z2)]
    if (rescale) {
      gamma.Z2 <- t(t(gamma.Z2)/x$design$sd.Z2)
    }
    names(gamma.Z2) <- vars.Z2
    print(gamma.Z2)
    cat("\n")
  }
  
  cat("---", "\n")
  
  cat("\n")
  
  cat("Optimal lambda:", x$lambda[which.min(x$criterion)], 
    "\n")
  
  cat("\n")
  
  cat("Log likelihood:", x$logLik[which.min(x$criterion)], 
    "\n")
  
  
  coef.opt <- list(theta = theta, intercepts = intercepts, 
    order.effects = order.effects, gamma.X = gamma.X, gamma.Z1 = gamma.Z1, 
    gamma.Z2 = gamma.Z2)
  
  invisible(coef.opt)
  
  
}

