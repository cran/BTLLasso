#' Print function for cv.BTLLasso objects
#' 
#' Prints the most important output of cv.BTLLasso objects.
#' 
#' @method print cv.BTLLasso
#' @param x \code{cv.BTLLasso} object
#' @param \dots possible further arguments for print command
#' @return \item{x}{cv.BTLLasso object}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{cv.BTLLasso}}
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
print.cv.BTLLasso <- function(x, ...){

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

  cv.crit <- x$cv.crit

  
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
  cat("\n Cross-validation criterion:",cv.crit,"\n")
  
  cat("---","\n")
  
  cat("Parameter estimates after",x$folds,"-","fold cross-validation","\n")
  
  cat("\n")
  coefs <- x$coefs.repar[which.min(x$criterion),]
  
  theta <- intercepts <- order.effects <- gamma.X <- gamma.Z1 <- gamma.Z2 <- c()
  
  if(n.theta>0){
    cat("thresholds:","\n")
    theta <- coefs[1:n.theta]
    names(theta) <- paste0("theta",1:n.theta)
    print(theta)
    cat("\n")
  }
  
  if(n.order>0){
    cat(paste0(x$control$name.order,":"),"\n")
    orders <- coefs[(n.theta+1):(n.theta + n.order)]
    if(n.order==m){
      names(orders) <- labels
    }
    if(n.order==1){
      names(orders) <- NULL
    }
    print(orders)
    cat("\n")
  }
  
  if(n.intercepts>0){
    cat("intercepts:","\n")
    intercepts <- coefs[(n.theta+ n.order+1):(n.theta + n.order+ n.intercepts)]
    names(intercepts) <- labels
    print(intercepts)
    cat("\n")
  }
  
  if(p.X>0){
    cat("object-specific effects for subject-specific covariate(s):","\n")
    gamma.X <- matrix(coefs[(n.theta+ n.order+ n.intercepts + 1): 
                              (n.theta+ n.order + n.intercepts + p.X*m)], nrow = p.X, byrow =TRUE)
    colnames(gamma.X) <- labels
    rownames(gamma.X) <- vars.X
    print(gamma.X)
    cat("\n")
  }
  
  if(p.Z1>0){
    cat("object-specific effects for subject-object-specific covariate(s):","\n")
    gamma.Z1 <- matrix(coefs[(n.theta+ n.order+ n.intercepts + p.X * m + 1): 
                               (n.theta+ n.order+ n.intercepts + p.X * m + p.Z1 * m)], nrow = p.Z1, byrow =TRUE)
    colnames(gamma.Z1) <- labels
    rownames(gamma.Z1) <- vars.Z1
    print(gamma.Z1)
    cat("\n")
  }
  
  if(p.Z2>0){

    cat("global effects for (subject-)object-specific covariate(s):","\n")
    gamma.Z2 <- coefs[(n.theta+ n.order+ n.intercepts + p.X * m + p.Z1 * m + 1): 
                        (n.theta+ n.order+ n.intercepts + p.X * m + p.Z1 * m + p.Z2)]
    names(gamma.Z2) <- vars.Z2
    print(gamma.Z2)
    cat("\n")
  }
  
  cat("---","\n")
  
  cat("\n")
  
  cat("Optimal lambda:",x$lambda[which.min(x$criterion)],"\n")
  
  cat("\n")
  
  cat("log likelihood:",x$logLik[which.min(x$criterion)],"\n")
  
  invisible(x)
  
  
}

