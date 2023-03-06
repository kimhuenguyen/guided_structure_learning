library('foreach')
library('doParallel')
library('igraph')
library('gbm')
library('bnlearn')
library('pcalg')
library('XMRF')

source("scripts/simData.R")
source("scripts/DiscretizationF.R")
source("scripts/ODS.R")
source("scripts/pdn.R")
source("scripts/result_scores.R")
source("scripts/runalgorithm_guide.R")

#######load Adj matrix
load("data/simscalefreeDAG10.Rdata")

##########
Adj <- Pmat
Adj[Adj!=0] <- 1
p <- ncol(Pmat) 
n <- 1000       
nsim <- 50    
# level of the test
alpha <- 2*pnorm(n^.15,lower.tail=F)
# upper bound for cardinaities of parent sets
npa <- 8        

result.al <- result_algorithms(n, p,Adj,Pmat, nsim, npa, alpha,cpu)
round(apply(result.al,2,mean),2)[c(4,5,6,9,15,16,17,20,26,27,28,31,37,38,39,42,48,49,50,53,59,60,61,64,70,71,72,75,81,82,83,86)]

