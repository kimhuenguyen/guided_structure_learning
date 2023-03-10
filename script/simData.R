#' Generate Poisson data
#'
#' This function generates a data matrix (\eqn{n \times p}) under Poisson
#' assumption from the coefficient matrix Pmat, and the topological ordering order.
#' @param n number of samples
#' @param p number of variables (nodes)
#' @param Pmat coefficient matrix
#' @param order the topological ordering of variables (names of nodes)
#' @importFrom
#' @examples
#' p <- 4
#' Pmat <- matrix(c(0, 0.4, -0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0, 0, 0), nrow=p,byrow=T)
#' n <- 10
#' order <- seq(1,p,1)
#' simData(n, p, order,Pmat)
#' @export
simData <- function(n,p,order,Pmat){
  all.simulated <- FALSE
  data <- matrix(0,nrow = n,ncol = p)
  cmax <- apply(abs(Pmat),2,max)
  temp <- which(cmax == 0)
  while(all.simulated == FALSE ){
    for (i in 1: p){
      if (sum(abs(Pmat[,order[i]]))==0)
        data[,order[i]] <- rpois(n,1)
      else
        data[,order[i]] <- rpois(n,exp( data %*% Pmat[,order[i]]))
      }
    all.simulated <- all(!is.na(data))
    }
data
}

