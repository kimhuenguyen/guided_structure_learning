### sample undirected graphs
library(igraph)
################# sample scalefree graph 
graph <- barabasi.game(n=10, m = 2, out.dist = NULL, out.seq = NULL, out.pref = FALSE, 
                       directed=FALSE)
plot(graph)

SAdj <- as_adjacency_matrix(graph)
sum(SAdj)

# save(SAdj, file = "scalefreesample10.RData")

################# sample random graph 

g <-erdos.renyi.game(n=10, p=0.3, type=c("gnp"),
                     directed = FALSE, loops = FALSE)
SAdj <- as_adj(g)
plot(g)
sum(SAdj)/2


# save(SAdj, file = "randomsample10-03.RData")

################# sample hub graph 
library(XMRF)
SAdj <- XMRF.Sim(n = 100, p = 10, model = "LPGM", graph.type = "hub")$B
graph <- graph_from_adjacency_matrix(SAdj, mode =  "undirected", weighted = NULL,
                                     diag = TRUE, add.colnames = NULL, add.rownames = NA)
plot(graph)
sum(SAdj)

# save(SAdj, file = "hubsample100.RData")

#' Sample directed acyclic graphs (DAG) and theirs coefficient matrices from 
#' undirected graphs
################# function to generate coefficient matrix of DAG
#' @param skel the matrix of skeleton
#' @param lb the lower bound of coefficients
#' @param ub the upper bound of coefficients
#' @param ord the topological ordering of nodes
#' @return the coefficient matrix of the graph.
#' @export
#' @importFrom 
#' example
#' sample an order of nodes
#' p <- ncol(SAdj)
#' order <- sample(1:p,p,replace = FALSE)
#' Pmat <- genDAG(SAdj,0.1,0.5,order)
genDAG <-  function(skel,lb,ub,ord){
  # create the adjacency matrix of DAG from set of skeleton and an order of nodes
  p <- ncol (skel)
  for (i in 1:p){
    for (j in 1:p){
      if (which (ord == i)< which(ord==j)) skel[j,i] = 0
    }
  }
  # create a matrix of coefficients
  nedge <- sum(skel)
  g <- function(nedge,lb,ub){
    weight <- c()
    ind <- rbinom(nedge,1,1/4)
    for (i in 1:length(ind)){
      if (ind[i] == 0){
        weight[i] <- -round(runif(1,lb,ub),3)
      }
      else{
        
        weight[i] <- round(runif(1,lb,ub),3)
      }
    }
    return(weight) 
  }
  skel[skel==1] <- g(nedge,lb,ub)
  return(skel)
}


