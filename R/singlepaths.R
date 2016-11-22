#' Plot parameter paths for BTLLasso
#' 
#' Plots single paths for every parameter of a BTLLasso object or a cv.BTLLasso
#' object. In contrast, to \code{\link{paths}}, one plot per covariate is
#' created, every single parameter is illustrated by one path. For cv.BTLLasso
#' object, the optimal model according to the cross-validation is marked by a
#' vertical dashed line.
#' 
#' Plots for BTLLasso and cv.BTLLasso objects only differ by the additional
#' vertical line indicating the optimal model according to cross-validation.
#' 
#' @param model BTLLasso or cv.BTLLasso object
#' @param colors Optional. If specified, vector with length equal to the number
#' of objects. Each object can be represented by another color.
#' @param equal.ranges Should all single plots (for different covariates) have
#' equal ranges on the y-axes. FALSE by default.
#' @param plot.X Should paths for variables in \code{X} (if present) be plotted?
#' @param plot.Z1 Should paths for variables in \code{Z1} (if present) be plotted?
#' @param plot.Z2 Should paths for variables in \code{Z2} (if present) be plotted?
#' @param plot.intercepts Should paths for intercepts be plotted separately?
#' @param plot.order.effects Should paths for order effects be plotted separately?
#' @param x.axis Should the paths be plotted against the (scaled) penalty term (between 0 and 1),
#' against lambda or against log(lambda+1)?
#' @param columns How many columns should be used to arrange the plots. If unspecified, plots 
#' are arranged automatically in a quadratic manner.
#' @param subs.X Optional vector of subtitles for variables in \code{X}. Can be used
#' to note the encoding of the single covariates, especially for dummy
#' variables.
#' @param subs.Z1 Optional vector of subtitles for variables in \code{Z1}. Can be used
#' to note the encoding of the single covariates, especially for dummy
#' variables.
#' @param subs.Z2 Optional vector of subtitles for variables in \code{Z2}. Can be used
#' to note the encoding of the single covariates, especially for dummy
#' variables.
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}},
#' \code{\link{paths}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' @keywords BTLLasso paths parameter paths
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
#' @export singlepaths
singlepaths <- function(model, colors = NULL, equal.ranges = FALSE,
                        plot.X = TRUE, plot.Z1 = TRUE, plot.Z2 = TRUE, plot.intercepts = TRUE,
                        plot.order.effects = TRUE, x.axis = c("penalty", "lambda", "loglambda"), 
                        columns = NULL, subs.X = NULL, subs.Z1 = NULL, subs.Z2 = NULL){

  x.axis <- match.arg(x.axis)

  coefs <- model$coefs

  
  acoefs <- model$penalty$acoefs
  norm <- rowSums(abs(coefs%*%acoefs))
  norm <- norm/max(norm)
  norm.range <- range(norm)
  
  if(x.axis == "lambda"){
  norm <- model$lambda
  norm.range <- rev(range(norm))
  }
  
  if(x.axis == "loglambda"){
    norm <- log(model$lambda+1)
    norm.range <- rev(range(norm))
  }

  
  m <- model$Y$m
  n.theta <- model$design$n.theta
  n.order <- model$design$n.order
  n.intercepts <- model$design$n.intercepts
if(n.intercepts>0){
  n.intercepts <- n.intercepts+1
}


p.X <- model$design$p.X
p.Z1 <- model$design$p.Z1
p.Z2 <- model$design$p.Z2
  
  labels <- model$Y$object.names
  
  coefs <- model$coefs.repar

  if(n.theta>0){
    theta <- coefs[,1:n.theta,drop=FALSE]
  }
  
  order.effects <- intercepts <- gamma.X <- gamma.Z1 <- gamma.Z2 <- c()

  if(n.order>0){
    order.effects <- coefs[,(n.theta+1):(n.theta+n.order)]
  }
  
  if(n.intercepts>0){
  intercepts <- coefs[,(n.theta+n.order+1):(n.theta+n.order + n.intercepts),drop=FALSE]
  }
  
  if(p.X>0){
  gamma.X <- coefs[,(n.theta+n.order+ n.intercepts + 1): (n.theta+n.order + n.intercepts + p.X*m)]
  }

  if(p.Z1>0){
  gamma.Z1 <- coefs[,(n.theta+n.order+ n.intercepts + p.X * m + 1): (n.theta+n.order+ n.intercepts + p.X * m + p.Z1 * m)]
  }
  
  if(p.Z2>0){
  gamma.Z2 <- coefs[,(n.theta+n.order+ n.intercepts + p.X * m + p.Z1 * m + 1): 
                      (n.theta+n.order+ n.intercepts + p.X * m + p.Z1 * m + p.Z2), drop = FALSE]
  }

  p <- p.tot <- 0
  gamma <- c()
  covar <- c()
  all.subs <- c()

  if(plot.X){
    gamma <- cbind(gamma, gamma.X)
    p <- p + p.X
    p.tot <- p.tot + p.X
  covar <- c(covar, model$design$vars.X)
  if(is.null(subs.X)){
    subs.X <- rep("",p.X)
  }
  all.subs <- c(all.subs,subs.X)
  }

  if(plot.Z1){
    gamma <- cbind(gamma, gamma.Z1)
    p <- p + p.Z1
    p.tot <- p.tot + p.Z1
    covar <- c(covar, model$design$vars.Z1)
    if(is.null(subs.Z1)){
      subs.Z1 <- rep("",p.Z1)
    }
    all.subs <- c(all.subs,subs.Z1)
  }
  
  if(plot.Z2){
    p.tot <- p.tot + p.Z2
    covar <- c(covar, model$design$vars.Z2)}

  
  if(plot.intercepts & n.intercepts>0){
    covar <- c("Intercept", covar)
    gamma <- cbind(intercepts, gamma)
    p <- 1+ p
    p.tot <- p.tot + 1
    all.subs <- c("",all.subs)
  }
  
  if(plot.order.effects & n.order>0){
    p.tot <- p.tot + 1
    if(n.order>1){
      p <- p+1
    gamma <- cbind(order.effects, gamma)
    covar <- c(model$control$name.order, covar) 
    }
  }
  
if(is.null(columns)){
  cols <- floor(sqrt(p.tot))
}else{
  cols <- columns
}
  rows <- ceiling((p.tot)/cols)
  
  layout(matrix(1:(rows*cols),nrow=rows,byrow=TRUE))

  x.axis.name <- x.axis
  if(x.axis == "loglambda"){
  x.axis.name <- expression(log(lambda+1))
  }
  if(x.axis == "lambda"){
    x.axis.name <- expression(lambda)
  }

  
  y.range <- range(gamma)
  
  if(!is.null(model$criterion)){
    x.axis.min <-norm[which.min(model$criterion)]
  }
  
  if(plot.order.effects & n.order==1){
    plot(norm, order.effects,type="l",main=model$control$name.order, ylab="",
         xlab=x.axis.name, xlim = norm.range)
    if(!is.null(model$criterion)){
      abline(v=x.axis.min,lty=2,col=2)
    }
  }

  if(is.null(colors)){colors <- rep(1,m)}

  index <- 1
  for(i in 1:p){
    if(!equal.ranges){y.range <- range(gamma[,index:(index+m-1)])}
    
    plot(norm, gamma[,index],ylim=y.range,type="l",main="",ylab="",
         xlab=x.axis.name, xlim = norm.range,
         col=colors[1])

    for(u in 1:(m-1)){
      lines(norm,gamma[,index+u],col=colors[u+1])
    }
    axis(4,at=gamma[nrow(gamma),index:(index+m-1)],labels = labels, las=2)
    title(main=covar[i],line=1.2)
    mtext(all.subs[i],side=3,line=0.2,cex=par()$cex)
    if(!is.null(model$criterion)){
    abline(v=x.axis.min,lty=2,col=2)
    }
    index <- index+m
  }

  if(p.Z2>0 & plot.Z2){
    for(i in 1:p.Z2){
      plot(norm, gamma.Z2[,i],type="l",main=model$design$vars.Z2[i], ylab="",
           xlab=x.axis.name,  xlim = norm.range)
      if(!is.null(model$criterion)){
        abline(v=x.axis.min,lty=2,col=2)
      }
      if(!is.null(subs.Z2)){
        mtext(subs.Z2[i],side=3,line=0.2,cex=par()$cex)
      }
    }

  }
  
  layout(1)
}
