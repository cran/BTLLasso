#' Plot confidence intervals for BTLLasso
#' 
#' Plots confidence intervals for every single coefficient. Confidence
#' intervals are separated by covariates, every covariate is plotted
#' separately. Confidence intervals are based on bootstrap, performed by
#' \code{\link{boot.BTLLasso}}.
#' 
#' 
#' @param object boot.BTLLasso object
#' @param plot.X Should confidence intervals for variables in \code{X} (if present) be plotted?
#' @param plot.Z1 Should confidence intervals for variables in \code{Z1} (if present) be plotted?
#' @param plot.Z2 Should confidence intervals for variables in \code{Z2} (if present) be plotted?
#' @param plot.intercepts Should confidence intervals for intercepts be plotted separately?
#' @param plot.order.effects Should confidence intervals for order effects be plotted separately?
#' @param include.zero Should all plots contain zero?
#' @param columns Optional argument for the number of columns in the plot.
#' @param subs.X Optional vector of subtitles for variables in \code{X}. Can be used
#' to note the encoding of the single covariates, especially for dummy
#' variables.
#' @param subs.Z1 Optional vector of subtitles for variables in \code{Z1}. Can be used
#' to note the encoding of the single covariates, especially for dummy
#' variables.
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{boot.BTLLasso}}, \code{\link{BTLLasso}},
#' \code{\link{cv.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' @keywords BTLLasso confidence interval bootstrap
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
#' @export ci.BTLLasso
ci.BTLLasso <- function(object, plot.X = TRUE, plot.Z1 = TRUE, plot.Z2 = TRUE, 
                        plot.intercepts = TRUE, plot.order.effects = TRUE,
                        include.zero = TRUE, columns = NULL, subs.X = NULL, subs.Z1 = NULL ){

  model <- object$cv.model
  epsilon <- model$control$epsilon
  accuracy <- -log10(epsilon)
  covariates <- c(model$design$vars.X,model$design$vars.Z1,model$design$vars.Z2)
  conf.ints <- object$conf.ints.repar

  m <- model$Y$m
  labels <- model$Y$object.names

  n.theta <- model$design$n.theta
  n.order <- model$design$n.order
  n.intercepts <- model$design$n.intercepts
  if(n.intercepts>0){
    n.intercepts <- n.intercepts+1
  }
  p.X <- model$design$p.X
  p.Z1 <- model$design$p.Z1
  p.Z2 <- model$design$p.Z2

  estimates <- model$coefs.repar[which.min(model$criterion),]
  estimates <- round(estimates,accuracy)
  conf.ints <- round(conf.ints,accuracy)
  
  p <- p.global <- 0
  gamma <- gamma.ci <- global <- global.ci <- c()
  covar <- covar.global <- c()
  
  start <- n.theta+1
  
  all.subs <- c()
  
  if(plot.order.effects & n.order>0){
    end <- start+n.order-1
    if(n.order==1){
      global <- c(global,estimates[start:end])
      global.ci <- cbind(global.ci,conf.ints[,start:end])
      covar.global <- c(covar.global, model$control$name.order)
      p.global <- p.global+1
    }
    if(n.order>1){
      gamma <- c(gamma,estimates[start:end])
      gamma.ci <- cbind(gamma.ci,conf.ints[,start:end])
      covar <- c(covar, model$control$name.order)
      p <- p+1
      all.subs <- c(all.subs,"")
    }
  }
  
  start <- n.theta+n.order+1
  
  if(n.intercepts>0 & plot.intercepts){
    end <- start+n.intercepts-1
    covar <- c(covar, "Intercept")
    gamma <- c(gamma, estimates[start:end])
    gamma.ci <- cbind(gamma.ci, conf.ints[,start:end])
    p <- 1+ p
    all.subs <- c(all.subs,"")
  }
  
  start <- n.theta+n.order+n.intercepts+1
  
  if(p.X>0 & plot.X){
    end <- start+p.X*m-1
    covar <- c(covar, model$design$vars.X)
    gamma <- c(gamma, estimates[start:end])
    gamma.ci <- cbind(gamma.ci, conf.ints[,start:end])
    p <- p+ p.X
    if(is.null(subs.X)){
      subs.X <- rep("",p.X)
    }
    all.subs <- c(all.subs,subs.X)
  }
  
  start <- n.theta+n.order+n.intercepts+p.X*m+1
  
  if(p.Z1>0 & plot.Z1){
    end <- start+p.Z1*m-1
    covar <- c(covar, model$design$vars.Z1)
    gamma <- c(gamma, estimates[start:end])
    gamma.ci <- cbind(gamma.ci, conf.ints[,start:end])
    p <- p+ p.Z1
    if(is.null(subs.Z1)){
      subs.Z1 <- rep("",p.Z1)
    }
    all.subs <- c(all.subs,subs.Z1)
  }

  
  start <- n.theta+n.order+n.intercepts+p.X*m+p.Z1*m+1
  
  if(p.Z2>0 & plot.Z2){
    end <- start+p.Z2-1
    covar.global <- c(covar.global, model$design$vars.Z2)
    global <- c(global, estimates[start:end])
    global.ci <- cbind(global.ci, conf.ints[,start:end])
    p.global <- p.global+ p.Z2
  }
  
p.tot <- p
  if(p.global>0){
    p.tot <- p.tot+1
  }
    
    
if(is.null(columns)){
  cols <- floor(sqrt(p.tot))
}else{
  cols <- columns
}
rows <- ceiling((p.tot)/cols)
    
    
      layout(matrix(1:(rows*cols),nrow=rows,byrow=TRUE))
  
    index <- 1
    for(i in 1:p){
      
      xlim <- range(gamma.ci[,index:(index+m-1)])
      if(include.zero){
        xlim <- range(0,xlim)
      }
      plot(gamma[index:(index+m-1)],1:m,xlim=xlim,pch=16,yaxt="n",xlab="",ylab="",
           main="")  
      segments(y0=1:m,x0=gamma.ci[1,index:(index+m-1)],x1=gamma.ci[2,index:(index+m-1)])
      axis(2,at =1:m, labels=labels,las=2)
      title(covar[i],line=1.2)
      mtext(all.subs[i],side=3,line=0.2,cex=par()$cex)
      abline(v=0,lty=2,col="lightgray")
      index <- index + m
   }
      index 
      if(p.global>0){
        xlim <- range(global.ci)
        if(include.zero){
          xlim <- range(0,xlim)
        }
        plot(global,1:p.global,xlim=xlim,
             pch=16,yaxt="n",xlab="",ylab="", main="Global Parameters")  
        segments(y0=1:p.global,x0=global.ci[1,],x1=global.ci[2,])
        axis(2,at =1:p.global, labels=covar.global,las=2)
        abline(v=0,lty=2,col="lightgray")
      }
  
  layout(1)

}

