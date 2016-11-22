#' Print function for boot.BTLLasso objects
#' 
#' Prints the most important output of boot.BTLLasso objects.
#' 
#' @method print boot.BTLLasso
#' @param x \code{boot.BTLLasso} object
#' @param \dots possible further arguments for print command
#' @return \item{x}{boot.BTLLasso object}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{boot.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' @keywords BTLLasso
#' @examples
#' 
#' \dontrun{
#' ## load data set
#' data(GLESsmall)
#' 
#' # define response and covariate matrix
#' X <- scale(GLESsmall[, 11:14])
#' Y <- as.matrix(GLESsmall[, 1:10])
#' 
#' # vector of subtitles, containing the coding of the single covariates
#' subs <- c("(in years)","female (1); male (0)",
#' "East Germany (1); West Germany (0)","(very) good (1); else (0)")
#' 
#' # vector of tuning parameters
#' lambda <- exp(seq(log(31),log(1),length=50))-1
#' 
#' 
#' # compute 10-fold cross-validation
#' set.seed(5)
#' m.cv <- cv.BTLLasso(Y = Y, X = X, folds = 10, lambda = lambda, cores = 10)
#' 
#' print(m.cv)
#' }
print.boot.BTLLasso <- function(x, ...){
  
  model <- x$cv.model
  epsilon <- model$control$epsilon
  accuracy <- -log10(epsilon)
  covariates <- c(model$design$vars.X,model$design$vars.Z1,model$design$vars.Z2)
  conf.ints <- x$conf.ints.repar
  
  m <- model$Y$m
  labels <- model$Y$object.names
  n.theta <- model$design$n.theta
  estimates <- model$coefs.repar[which.min(model$criterion),]
  estimates <- round(estimates,accuracy)
  conf.ints <- round(conf.ints,accuracy)
  
  print.mat <- rbind(conf.ints[1,],estimates,conf.ints[2,])
  rownames(print.mat)[c(1,3)] <- rownames(conf.ints)
  rownames(print.mat)[2] <- ""
  
  cat("Confidence intervals according to bootstrap:","\n")
  
  print(print.mat)
  
  invisible(x)
  
}

