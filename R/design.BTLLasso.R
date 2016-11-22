design.BTLLasso <- function(Y, X=NULL, Z1=NULL, Z2=NULL, control = ctrl.BTLLasso()){
 
  ####
  ## get all arguments from responseBTLLasso object
  ####
  
  y.ord <- Y$response 
  first.object <- Y$first.object
  second.object <- Y$second.object
  subject <- Y$subject
  withS <- Y$withS
  subject.names <- Y$subject.names  
  object.names <- Y$object.names 
  n <- Y$n 
  m <- Y$m
  k <- Y$k
  q <- Y$q
  
  ## get some control arguments
  include.intercepts <- control$include.intercepts
  order.effect <- control$order.effect 
  object.order.effect <- control$object.order.effect 
  
  
  ####
  ## check X, Z1 and Z2 and initialize p.X, p.Z1 and p.Z2
  ####
  withX <- withZ1 <- withZ2 <- FALSE
  p.X <- p.Z1 <- p.Z2 <- 0
  par.names.X <- par.names.Z1 <- par.names.Z2 <- c()
  vars.X <- vars.Z1 <- vars.Z2 <- c()
  acoefs.X <- acoefs.Z1 <- acoefs.Z2 <- c()
  order.Z1 <- order.Z2 <- order(object.names)
  
  if(!is.null(X)){
    withX <- TRUE
    p.X <- ncol(X)
    vars.X <- colnames(X)
    if(!is.matrix(X))
      stop("X has to be a matrix")
    if(control$scale){
      X <- scale(X)
    }
    par.names.X <- paste(rep(vars.X, each = m-1),object.names[1:(m-1)],sep=".")
  }
  
  
  if(!is.null(Z1)){
    withZ1 <- TRUE
    p.Z1 <- ncol(Z1)/m
    if(ncol(Z1)%%m!=0)
      stop("Number of columns of Z1 has to be a multiple of the number of objects")
    if(!is.matrix(Z1))
      stop("Z1 has to be a matrix")
    if(control$scale){
      Z1 <- matrix(c(scale(matrix(c(Z1),ncol=p.Z1))),ncol=ncol(Z1),
                   dimnames = list(rownames(Z1),colnames(Z1)))
    }
    check.Z1 <- check(Z = Z1, object.names = object.names, subject)
    vars.Z1 <- check.Z1$vars.Z
    order.Z1 <- check.Z1$order.Z
    par.names.Z1 <- paste(rep(vars.Z1, each = m),object.names[1:m],sep=".")
  }
  
  if(!is.null(Z2)){
    withZ2 <- TRUE
    p.Z2 <- ncol(Z2)/m
    check.Z2 <- check(Z = Z2, object.names = object.names, subject)
    vars.Z2 <- check.Z2$vars.Z
    order.Z2 <- check.Z2$order.Z
    if(control$scale){
      Z2 <- matrix(c(scale(matrix(c(Z2),ncol=p.Z2))),ncol=ncol(Z2),
                   dimnames = list(rownames(Z2),colnames(Z2)))
    }
    par.names.Z2 <- vars.Z2
  }
  
  
  
  ## number of intercepts
  n.intercepts <- 0
  par.names.intercepts <- c()
  if(include.intercepts){
    n.intercepts <- m-1
    par.names.intercepts <- object.names[1:(m-1)]
  }
  
  ## number of order effects
  n.order <- 0
  par.names.order <- c()
  if(order.effect){
    n.order <- 1
    par.names.order <- control$name.order
  }
  if(object.order.effect){
    n.order <- m
    par.names.order <- paste(control$name.order, object.names, sep=".")
  }
  
  ####
  ## make design matrix
  ####
  ## design matrix
  design <- create.design(X, Z1, Z2, first.object, second.object, m, subject, control, order.Z1, order.Z2)  

  ## enlarge design matrix so that it fits to the dichotomized cumulative response
  design <- t(matrix(rep(c(design),each=q),nrow=ncol(design),byrow=TRUE))
  colnames(design) <- c(par.names.order, par.names.intercepts, par.names.X, par.names.Z1, par.names.Z2)
  
  
  n.theta <- floor(q/2)
  
  if(k>2){
    theta.design <- matrix(0,ncol=n.theta,nrow=nrow(design))
    colnames(theta.design) <- paste("theta",1:n.theta,sep=".")
    for(i in 1:n.theta){
      vec1 <- rep(0,q)
      vec1[c(i,q-i+1)] <- c(1,-1)
      theta.design[,i] <- rep(vec1,length(y.ord))
    }
    
    design <- cbind(theta.design, design)
  }
  
  
  RET <- list(design = design, p.X = p.X, p.Z1 = p.Z1, p.Z2 = p.Z2, 
              vars.X = vars.X, vars.Z1 = vars.Z1, vars.Z2 = vars.Z2, n.theta = n.theta,
              n.intercepts = n.intercepts, n.order = n.order)
  
  return(RET)
}