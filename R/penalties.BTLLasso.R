penalties.BTLLasso <- function(Y, X=NULL, Z1=NULL, Z2=NULL, control = ctrl.BTLLasso()){
  ####
  ## get arguments from responseBTLLasso object
  ####
  n <- Y$n 
  m <- Y$m
  k <- Y$k
  q <- Y$q
  object.names <- Y$object.names 
  
  ####
  ## extract all control arguments
  ####
  penalize.X <- control$penalize.X
  penalize.Z1.diffs <- control$penalize.Z1.diffs
  penalize.Z1.absolute <- control$penalize.Z1.absolute
  penalize.Z2 <- control$penalize.Z2
  penalize.intercepts <- control$penalize.intercepts
  include.intercepts <- control$include.intercepts
  order.effect  <-  control$order.effect
  object.order.effect <- control$object.order.effect
  penalize.order.effect.diffs <- control$penalize.order.effect.diffs
  penalize.order.effect.absolute <- control$penalize.order.effect.absolute
  
  
  ## number of intercepts
  n.intercepts <- 0
  par.names.intercepts <- c()
  if(include.intercepts){
    n.intercepts <- m-1
    par.names.intercepts <- object.names[1:(m-1)]
  }
  
  ## number of order effects
  n.order <- 0
  if(order.effect){
    n.order <- 1
  }
  if(object.order.effect){
    n.order <- m
  }
  
  ####
  ## create penalty matrices
  ####
  numpen.intercepts <- numpen.X <- numpen.Z1 <- numpen.Z2 <- numpen.order <- 0
  p.X <-p.Z1 <- p.Z2 <- 0
  
  ## penalty matrix for intercepts
  if(include.intercepts & penalize.intercepts){
    acoefs.intercepts <- diag(m-1)
    
    help.pen <- matrix(0,ncol=choose(m-1,2),nrow=m-1)
    combis <- combn(m-1,2)
    for(ff in 1:ncol(combis)){
      help.pen[combis[1,ff],ff] <- 1
      help.pen[combis[2,ff],ff] <- -1
    }
    acoefs.intercepts <- cbind(acoefs.intercepts,help.pen)
    numpen.intercepts <- ncol(acoefs.intercepts)
  }
  
  ## penalty matrix for X
  if(!is.null(X)){
     p.X <- ncol(X) 
     if(penalize.X){
    acoefs.X <- diag(p.X*(m-1))
    help.pen <- matrix(0,ncol=choose(m-1,2),nrow=m-1)
    combis <- combn(m-1,2)
    for(ff in 1:ncol(combis)){
      help.pen[combis[1,ff],ff] <- 1
      help.pen[combis[2,ff],ff] <- -1
    }
    for(pp in 1:p.X){
      m.above <- matrix(rep(matrix(0,ncol=choose(m-1,2),nrow=m-1),pp-1),ncol=choose(m-1,2))
      m.below <- matrix(rep(matrix(0,ncol=choose(m-1,2),nrow=m-1),p.X-pp),ncol=choose(m-1,2))
      acoefs.X <- cbind(acoefs.X,rbind(m.above,help.pen,m.below))
    } 
    numpen.X <- ncol(acoefs.X)
  }
  }

  ## penalty matrix for Z1
  if(!is.null(Z1)){ 
       p.Z1 <- ncol(Z1)/m   
       
       if(penalize.Z1.diffs | penalize.Z1.absolute){
    
    acoefs.Z1 <- c()
    if(penalize.Z1.absolute){
      acoefs.Z1 <- diag(m*p.Z1)
    }
    if(penalize.Z1.diffs){
      help.pen <- matrix(0,ncol=choose(m,2),nrow=m)
      combis <- combn(m,2)
      for(ff in 1:ncol(combis)){
        help.pen[combis[1,ff],ff] <- 1
        help.pen[combis[2,ff],ff] <- -1
      }
      for(pp in 1:p.Z1){
        m.above <- matrix(rep(matrix(0,ncol=choose(m,2),nrow=m),pp-1),ncol=choose(m,2))
        m.below <- matrix(rep(matrix(0,ncol=choose(m,2),nrow=m),p.Z1-pp),ncol=choose(m,2))
        acoefs.Z1 <- cbind(acoefs.Z1,rbind(m.above,help.pen,m.below))
      } 
    }
    numpen.Z1 <- ncol(acoefs.Z1)
  }
  }
  
  ## penalty matrix for Z2
  if(!is.null(Z2)) {
     p.Z2 <- ncol(Z2)/m
     if(penalize.Z2){
    acoefs.Z2 <- diag(p.Z2)
    numpen.Z2 <- ncol(acoefs.Z2)
  }
  }
  
  ## penalty matrix for order effects
  #### if global order effect
  if(order.effect & penalize.order.effect.absolute){
    acoefs.order <- matrix(1,ncol=1,nrow=1)
    numpen.order <- 1
  }
  
  #### if object-specific order effects
  if(object.order.effect & (penalize.order.effect.diffs | penalize.order.effect.absolute)){
    acoefs.order <- c()
    if(penalize.order.effect.absolute){
      acoefs.order <- diag(m)
    }
    if(penalize.order.effect.diffs){
      help.pen <- matrix(0,ncol=choose(m,2),nrow=m)
      combis <- combn(m,2)
      for(ff in 1:ncol(combis)){
        help.pen[combis[1,ff],ff] <- 1
        help.pen[combis[2,ff],ff] <- -1
      }
        acoefs.order <- cbind(acoefs.order,help.pen)
    }
    numpen.order <- ncol(acoefs.order)
  }
  

  
  
  ## total number of penalty
  numpen <- numpen.intercepts + numpen.X + numpen.Z1 + numpen.Z2 + numpen.order
  
  ## initalize total penalty matrix
  acoefs <- matrix(0, ncol=numpen, nrow=n.intercepts + p.X*(m-1) + p.Z1*m + p.Z2 + n.order)
  
  current.row <- 1
  current.col <- 1
  
  ## add penalties for order effects
  if(n.order>0){
    if(numpen.order > 0){
      acoefs[current.row:(current.row+n.order-1),current.col:(current.col+numpen.order-1)] <- acoefs.order
    }
    current.row <- current.row + n.order
    current.col <- current.col + numpen.order
  }
  
  ## add penalties for intercepts
  if(include.intercepts){
    if(penalize.intercepts){
      acoefs[current.row:(current.row+m-2),current.col:(current.col+numpen.intercepts-1)] <- acoefs.intercepts
    }
    current.row <- current.row + m - 1
    current.col <- current.col + numpen.intercepts
  }
  
  ## add penalties for X
  if(!is.null(X)){
    if(penalize.X){
      acoefs[current.row:(current.row+p.X*(m-1)-1), current.col:(current.col+numpen.X-1)] <- acoefs.X
    }
    current.row <- current.row + p.X*(m-1)
    current.col <- current.col + numpen.X
  }
  
  ## add penalties for Z1
  if(!is.null(Z1)){
    if(penalize.Z1.diffs | penalize.Z1.absolute){
      acoefs[current.row:(current.row+p.Z1*m-1), current.col:(current.col+numpen.Z1-1)] <- acoefs.Z1
    }
    current.row <- current.row + p.Z1*m
    current.col <- current.col + numpen.Z1
  }
  
  ## add penalties for Z1
  if(!is.null(Z2)){
    if(penalize.Z2){
      acoefs[current.row:(current.row+p.Z2-1), current.col:(current.col+numpen.Z2-1)] <- acoefs.Z2
    }
    current.row <- current.row + p.Z2
    current.col <- current.col + numpen.Z2
  }
  
  
  
  ## additional rows for thetas, thetas (thresholds) are never penalized
  acoefs <- rbind(matrix(0,nrow=floor(q/2),ncol=ncol(acoefs)),acoefs)
  
  RET <- list(acoefs = acoefs, numpen.intercepts = numpen.intercepts, numpen.X = numpen.X,
              numpen.Z1 = numpen.Z1, numpen.Z2 = numpen.Z2, numpen.order = numpen.order, n.order = n.order,
              p.X = p.X, p.Z1 = p.Z1, p.Z2 = p.Z2, weight.penalties = control$weight.penalties)
  
  return(RET)
}