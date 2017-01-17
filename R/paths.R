#' Plot covariate paths for BTLLasso
#' 
#' Plots paths for every covariate of a BTLLasso object or a cv.BTLLasso
#' object. In contrast, to \code{\link{singlepaths}}, only one plot is created,
#' every covariate is illustrated by one path. For cv.BTLLasso object, the
#' optimal model according to the cross-validation is marked by a vertical
#' dashed line.
#' 
#' Plots for BTLLasso and cv.BTLLasso objects only differ by the additional
#' vertical line indicating the optimal model according to cross-validation.
#' 
#' @param model BTLLasso or cv.BTLLasso object
#' @param y.axis Two possible values for the y-axis. Variables can either be plotted
#' with regard to their contribution to the total penalty term (\code{y.axis='penalty'}) or
#' with regard to the $L_2$ norm of the corresponding parameter vector (\code{y.axis='L2'}).
#' @param x.axis Should the paths be plotted against log(lambda+1), against lambda or
#' against the scaled penalty term (between 0 and 1)?
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}},
#' \code{\link{singlepaths}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2016): Modelling 
#' Football Results in the German Bundesliga Using Match-specific Covariates, 
#' \emph{Department of Statistics, LMU Munich}, Technical Report 197
#' @keywords BTLLasso paths covariate paths
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
#' 
#' @export paths
paths <- function(model, y.axis = c("penalty", "L2"), x.axis = c("loglambda", 
  "lambda", "penalty")) {
  
  x.axis <- match.arg(x.axis)
  y.axis <- match.arg(y.axis)
  
  if (y.axis == "penalty") {
    y.text <- "penalty size"
  }
  if (y.axis == "L2") {
    y.text <- "L2 norm"
  }
  
  coefs <- model$coefs
  covar <- c(model$design$vars.X, model$design$vars.Z1, model$design$vars.Z2)
  
  acoefs <- model$penalty$acoefs
  norm <- rowSums(abs(coefs %*% acoefs))
  norm <- norm/max(norm)
  norm.range <- range(norm)
  
  if (x.axis == "lambda") {
    norm <- model$lambda
    norm.range <- rev(range(norm))
  }
  
  if (x.axis == "loglambda") {
    norm <- log(model$lambda + 1)
    norm.range <- rev(range(norm))
  }
  
  
  
  x.axis.name <- x.axis
  if (x.axis == "loglambda") {
    x.axis.name <- expression(log(lambda + 1))
  }
  if (x.axis == "lambda") {
    x.axis.name <- expression(lambda)
  }
  
  m <- model$Y$m
  n.theta <- model$design$n.theta
  n.order <- model$design$n.order
  n.intercepts <- model$design$n.intercepts
  p.X <- model$design$p.X
  p.Z1 <- model$design$p.Z1
  p.Z2 <- model$design$p.Z2
  
  # acoefs.intercepts <- acoefs.X <- acoefs.Z1 <- acoefs.Z2 <-
  # c()
  numpen.order <- model$penalty$numpen.order
  numpen.intercepts <- model$penalty$numpen.intercepts
  numpen.X <- model$penalty$numpen.X
  numpen.Z1 <- model$penalty$numpen.Z1
  numpen.Z2 <- model$penalty$numpen.Z2
  
  labels <- model$Y$object.names
  
  criterion <- model$criterion
  
  order.effects <- intercepts <- gamma.X <- gamma.Z1 <- gamma.Z2 <- c()
  
  index.cols.X <- index.cols.Z1 <- index.cols.Z2 <- c()
  index.rows.X <- index.rows.Z1 <- index.rows.Z2 <- c()
  
  if (n.order > 0) {
    order.effects <- coefs[, (n.theta + 1):(n.theta + n.order)]
  }
  
  if (n.intercepts > 0) {
    intercepts <- coefs[, (n.theta + n.order + 1):(n.theta + 
      n.order + n.intercepts)]
  }
  
  p <- p.X + p.Z1 + p.Z2
  
  paths <- c()
  
  start.row <- n.theta + n.intercepts + n.order
  if (p.X > 0) {
    index <- rep((1:p.X), each = m - 1)
    for (i in 1:p.X) {
      if (y.axis == "penalty") {
        paths <- cbind(paths, rowSums(abs(coefs[, start.row + 
          which(index == i), drop = FALSE] %*% acoefs[start.row + 
          which(index == i), , drop = FALSE])))
      } else {
        paths <- cbind(paths, sqrt(rowSums(coefs[, start.row + 
          which(index == i), drop = FALSE]^2)))
      }
      
    }
    start.row <- start.row + length(index)
  }
  
  if (p.Z1 > 0) {
    index <- rep(1:p.Z1, each = m)
    for (i in 1:p.Z1) {
      if (y.axis == "penalty") {
        paths <- cbind(paths, rowSums(abs(coefs[, start.row + 
          which(index == i), drop = FALSE] %*% acoefs[start.row + 
          which(index == i), , drop = FALSE])))
      } else {
        paths <- cbind(paths, sqrt(rowSums(coefs[, start.row + 
          which(index == i), drop = FALSE]^2)))
      }
      
    }
    start.row <- start.row + length(index)
  }
  
  if (p.Z2 > 0) {
    index <- 1:p.Z2
    for (i in 1:p.Z2) {
      if (y.axis == "penalty") {
        paths <- cbind(paths, rowSums(abs(coefs[, start.row + 
          which(index == i), drop = FALSE] %*% acoefs[start.row + 
          which(index == i), , drop = FALSE])))
      } else {
        paths <- cbind(paths, sqrt(rowSums(coefs[, start.row + 
          which(index == i), drop = FALSE]^2)))
      }
      
    }
  }
  
  
  
  if (!is.null(criterion)) {
    x.axis.min <- norm[which.min(criterion)]
  }
  
  
  
  if (numpen.intercepts > 0) {
    if (y.axis == "penalty") {
      paths <- cbind(rowSums(abs(intercepts %*% acoefs[(n.theta + 
        n.order + 1):(n.theta + n.order + n.intercepts), 
        (numpen.order + 1):(numpen.order + numpen.intercepts)])), 
        paths)
    } else {
      paths <- cbind(sqrt(rowSums(intercepts^2)), paths)
    }
    
    p <- p + 1
    covar <- c("Intercept", covar)
  }
  
  if (numpen.order > 0) {
    if (y.axis == "penalty") {
      paths <- cbind(rowSums(abs(order.effects %*% acoefs[(n.theta + 
        1):(n.theta + n.order), 1:numpen.order])), paths)
    } else {
      paths <- cbind(sqrt(rowSums(order.effects^2)), paths)
    }
    p <- p + 1
    covar <- c(model$control$name.order, covar)
  }
  
  
  plot(norm, paths[, 1], type = "l", ylim = range(paths), ylab = y.text, 
    xlab = x.axis.name, xlim = norm.range, las = 1)
  for (o in 2:p) {
    lines(norm, paths[, o])
  }
  if (!is.null(criterion)) {
    abline(v = x.axis.min, lty = 2, col = 2)
  }
  
  axis(4, at = paths[nrow(paths), ], labels = covar, las = 2)
  
  
  
}
