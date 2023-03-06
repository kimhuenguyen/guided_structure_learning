#### discretize continuous measurements
library(mclust)
############################
##### model based clustering
#############################
##################################################
discMBC <- function(dataC){
  resClass <- apply(dataC, 2,
                    function(x) Mclust(x, G=2:3)$classification)
  resClass
}

##################################
#### equal interval classification
###################################
##################################################
discI <- function(dataC, ncat=2){
  resClass <- apply(dataC, 2, function(x) cut(x, breaks=ncat,
                                              labels=FALSE))
  resClass
}
##################################################


#### example
#Sigma <- matrix(c(1,0.56,0.21,0.56,1,0.42,0.21,0.42,1), nrow=3)
#dataC <- rmvnorm(12, sigma=Sigma)
#dataD <- disc(dataC)


#dataD <- discI(dataC, ncat=2)
#dataD2 <- discI(dataC, ncat=3)
