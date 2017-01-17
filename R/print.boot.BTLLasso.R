#' Print function for boot.BTLLasso objects
#' 
#' Prints the most important output of boot.BTLLasso objects.
#' 
#' @method print boot.BTLLasso
#' @param x \code{boot.BTLLasso} object
#' @param rescale Should the parameter estimates be rescaled for plotting? Only 
#' applies if \code{scale = TRUE} was specified in \code{BTLLasso} or \code{cv.BTLLasso}.
#' @param \dots possible further arguments for print command
#' @return \item{x}{boot.BTLLasso object}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{boot.BTLLasso}}
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
print.boot.BTLLasso <- function(x, rescale = FALSE, ...) {
  
  model <- x$cv.model
  epsilon <- model$control$epsilon
  accuracy <- -log10(epsilon)
  covariates <- c(model$design$vars.X, model$design$vars.Z1, 
    model$design$vars.Z2)
  conf.ints <- x$conf.ints.repar
  
  m <- model$Y$m
  labels <- model$Y$x.names
  
  n.theta <- model$design$n.theta
  n.order <- model$design$n.order
  n.intercepts <- model$design$n.intercepts
  if (n.intercepts > 0) {
    n.intercepts <- n.intercepts + 1
  }
  p.X <- model$design$p.X
  p.Z1 <- model$design$p.Z1
  p.Z2 <- model$design$p.Z2
  
  estimates <- model$coefs.repar[which.min(model$criterion), 
    ]
  estimates <- round(estimates, accuracy)
  conf.ints <- round(conf.ints, accuracy)
  
  gamma.total <- c()
  
  start <- 1
  end <- n.theta
  
  cat("Output after bootstrap estimation of confidence intervals:\n")
  
  cat("---", "\n")
  
  if (n.theta > 0) {
    cat("Thresholds:", "\n")
    gamma <- matrix(NA, nrow = 3, ncol = n.theta)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    gamma[c(1, 3), ] <- conf.ints[, start:end]
    gamma[2, ] <- estimates[start:end]
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  start <- n.theta + 1
  
  if (n.order > 0) {
    end <- start + n.order - 1
    gamma <- matrix(NA, nrow = 3, ncol = n.order)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat(paste0(model$control$name.order, ":"), "\n")
    gamma[c(1, 3), ] <- conf.ints[, start:end]
    gamma[2, ] <- estimates[start:end]
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  start <- n.theta + n.order + 1
  
  if (n.intercepts > 0) {
    end <- start + n.intercepts - 1
    gamma <- matrix(NA, nrow = 3, ncol = n.intercepts)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat("Intercepts:\n")
    gamma[c(1, 3), ] <- conf.ints[, start:end]
    gamma[2, ] <- estimates[start:end]
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  start <- n.theta + n.order + n.intercepts + 1
  
  if (p.X > 0) {
    end <- start + p.X * m - 1
    if (rescale) {
      est <- estimates[start:end]/rep(model$design$sd.X, 
        each = m)
      est.ci <- t(t(conf.ints[, start:end])/rep(model$design$sd.X, 
        each = m))
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    
    gamma <- matrix(NA, nrow = 3, ncol = p.X * m)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat("Object-specific effects for subject-specific covariate(s):", 
      "\n")
    gamma[c(1, 3), ] <- est.ci
    gamma[2, ] <- est
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  start <- n.theta + n.order + n.intercepts + p.X * m + 1
  
  
  
  if (p.Z1 > 0) {
    end <- start + p.Z1 * m - 1
    if (rescale) {
      est <- estimates[start:end]/rep(model$design$sd.Z1, 
        each = m)
      est.ci <- t(t(conf.ints[, start:end])/rep(model$design$sd.Z1, 
        each = m))
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    gamma <- matrix(NA, nrow = 3, ncol = p.Z1 * m)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat("Object-specific effects for subject-object-specific covariate(s):", 
      "\n")
    gamma[c(1, 3), ] <- est.ci
    gamma[2, ] <- est
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  
  start <- n.theta + n.order + n.intercepts + p.X * m + p.Z1 * 
    m + 1
  
  if (p.Z2 > 0) {
    end <- start + p.Z2 - 1
    if (rescale) {
      est <- estimates[start:end]/model$design$sd.Z2
      est.ci <- t(t(conf.ints[, start:end])/model$design$sd.Z2)
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    gamma <- matrix(NA, nrow = 3, ncol = p.Z2)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat("Global effects for (subject-)object-specific covariate(s):", 
      "\n")
    gamma[c(1, 3), ] <- est.ci
    gamma[2, ] <- est
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  
  invisible(gamma.total)
  
}

