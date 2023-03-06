
#' Structure learning with Poisson models using OverDispersion Scoring (ODS)
#'  algorithm (Park et al. 2015)
#'
#' This function estimates the adjacency matrix of a Poisson model given a
#' matrix of counts.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @return the estimated adjacency matrix of the graph.
#' @export
#' @importFrom stats coefficients
ODS <- function(X,cpu){
  z <- glm.control(epsilon = 1e-05, maxit = 100, trace = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  ####### Step 1: Estimate candidate parents sets using LPGM
  lpgm.fit <- XMRF(t(X), method="LPGM", N=50, stability = "star",beta = 0.1,
                   nCpus = cpu, parallel = FALSE,  th = 1e-06, 
                   nlams = 20, sth = 0.6,lmin = 0.01)
  adj <- lpgm.fit$network[[lpgm.fit$opt.index]]
  
  ######## Step 2: Estimate the topological ordering using OverDispersion Scoring
  covX <- cov(X)
  meanX <- apply(X,2,mean)
  sigma <- sapply(seq(1,p,1), function(i){
    covX[i,i] - meanX[i] 
  })
  ind <- which (abs(sigma) == min(abs(sigma)))
  eorder <- ind
  jstart <- length(eorder) +1
  for (j in jstart: (p-jstart+1)){
    nind <- setdiff(seq(1,p,1),eorder)
    T <- foreach(t = 1:length(nind), .combine="cbind") %dopar%{
      neigh.ind <- which(adj[, nind[t]]==1)
      cand.pa <- intersect(eorder, neigh.ind)
      if (length(cand.pa)>0){
        patodo <- cand.pa
      }else{
        patodo <- eorder
      }
      datatodo <- X[,c(patodo,nind[t])]
      res <- length(patodo)+1
      fit <- glm(datatodo[,res] ~ datatodo[,-res],family = "poisson",control = z)
      kmean <- fit$fitted.values
      kvar <- (datatodo[,res] - kmean)* t(datatodo[,res] - kmean)
      sjk <- sum(kvar-kmean)/nrow(datatodo)
      return(sjk)
    }
    ind <- which(abs(T) == min(abs(T)))
    eorder <- c(eorder,nind[ind])
  }
  remain <- setdiff(seq(1,p,1),eorder)
  eorder <- c(eorder,remain)
  ######## Step 3: Estimate directed edges using Lasso
  eAdj <- LPGMord(X,eorder,adj,p)
  eAdj[eAdj != 0] <- 1
  return (eAdj)
}

#' Function estimate DAG using Lasso from a matrix of counts
#' and topological ordering
#' @param X the matrix of counts (n times p).
#' @param order the topological ordering of variables (names of nodes)
#' @param adj the matrix (p times p) identifies potential parents (resulting from Step 1)
#' @param the number of random variables
#' @return the estimated adjacency matrix of the graph.
#' @export est_pa function
#' @importFrom glmnet 
LPGMord  <- function(X,order,adj,p){
  T <- foreach(i = 1:p, .combine="cbind",.packages = c("glmnet"),.export = c("est_pa")) %dopar%{
    fit  <- est_pa(i,X=X,order,adj,p)
    return(fit)
  }
  Beta   <- matrix(unlist(T),p,p)
  return(Beta)
}


#' Function estimate parent set using Lasso for each node
#' 
#' @param X the matrix of counts (n times p).
#' @param order the topological ordering of variables (names of nodes)
#' @param adj the matrix (p times p) identifies potential parents (resulting from Step 1)
#' @param p the number of random variables
#' @return the estimated adjacency matrix of the graph.
#' @export 
#' @importFrom glmnet 
est_pa <- function(i,X,order,adj,p){
  Beta <- matrix(0,p,p)
  if (length(colnames(X))>0){
    rownames(Beta) <- colnames(X)
  }else{
    rownames(Beta) <- as.character(seq(1,p,1))
  }
  order     <- as.numeric(order)
  temp      <- which(order==i)
  if(temp==1) {
    Beta[,i] <- 0
  }else{
    ind.Y <- intersect(order[1:(temp-1)],which(adj[,i]==1))
    if (length(ind.Y) >0){
      Y <- X[,ind.Y]
      if(length(ind.Y) ==1){
        fit <- glm(X[,i]~Y,family = "poisson")
        theta <- fit$coefficients  
      }else{
        fit <- cv.glmnet(Y,X[,i],family="poisson",type.measure = "deviance")
        theta <- coef(fit, s="lambda.min")
      }
      Beta[as.character(ind.Y),i] <- round(theta[-1],3)
    }else
      Beta[,i] <- 0
  }
  
  return(Beta[,i])
}



