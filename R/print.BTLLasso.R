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
print.BTLLasso <- function(x, ...){

  m <- x$Y$m
  n <- x$Y$n
  k <- x$Y$q +1
  n.theta <- x$design$n.theta
  n.intercepts <- x$design$n.intercepts
  if(n.intercepts!=0){n.intercepts <- n.intercepts + 1}
  n.order <- x$design$n.order
  p.X <- x$design$p.X
  p.Z1 <- x$design$p.Z1
  p.Z2 <- x$design$p.Z2
  lambda <- x$lambda
  
  vars.X <- x$design$vars.X
  vars.Z1 <- x$design$vars.Z1
  vars.Z2 <- x$design$vars.Z2
  
  labels <- x$Y$object.names
  
  cat("Output of BTL-Lasso estimation:","\n")
  
  cat("---","\n")

  cat("Setting:")
  cat("\n", n, "subjects")
  cat("\n", m, "objects")
  cat("\n", k, "response categories")
  cat("\n", p.X, "subject-specific covariate(s)")
  cat("\n", p.Z1, "subject-object-specific covariate(s) with object-specific effects")
  cat("\n", p.Z2, "(subject-)object-specific covariate(s) with global effects")
  if(n.order==m){
    cat("\n", n.order, "subject-specific order effects")
  }
  if(n.order==1){
    cat("\n", "Global order effect")
  }
  if(n.order==0){
    cat("\n", "No order effect")
  }
  cat("\n", length(lambda), "different tuning parameters","\n")
 
  
  invisible(x)
  
  
}

