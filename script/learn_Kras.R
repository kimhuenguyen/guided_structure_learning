
#### Analysis for selected 8 genes that belong to a connected component of the 
#### non-small cell lung cancer pathway
library(gRbase)
library(org.Hs.eg.db)
library(foreach)
library(igraph)
library('doParallel')
cl <- makeCluster(5)
registerDoParallel(cl)


## Load count matrix
load('data/Kras_dataset.RData')
counts <- t(counts)
dim(counts)
# The original data has 5618 genes (in the pathways), and 5505 observations.

## Genes of interest
map.genes

## Topological ordering of 8 genes from KEGG
newORD <- c("1839","1956","5781","2885","6654","3845","5290","5894")

## Select dataset consists of 8 interested genes
ind.gene <- which(colnames(counts) %in% newORD)
data_gpath <- counts[, ind.gene]
dim(data_gpath)
# The data is considered for further analysis has 8 genes and 5505 observations

## Filter data falls more than 3 standard deviantions from the mean
mean_data <- mean(data_gpath)
sd_data <- sd(data_gpath)
ind_skip <- c()
pa_th <- mean_data + 3*sd_data
for (i in 1:dim(data_gpath)[1]) {
  if(sum(data_gpath[i,]>pa_th)>0)
    ind_skip <- c(ind_skip,i)
}

X <- data_gpath[-ind_skip,]
dim(X)

##  Estimating gene networks by Or-PPGM 
library(learn2count)
maxcard <- 6
alpha <- 0.05

### Or-PPGM
res <-  pois.ord(X, maxcard, alpha, newORD)

## plot the resulting graph
colnames(res) <- rownames(res)
g.res <- as_graphnel(graph_from_adjacency_matrix(res))
name.gen <- c()
name.numb <- colnames(res)
for(i in 1:length(name.numb )){
  temp <- which(map.genes[,2]==name.numb [i])
  name.gen <- c(name.gen,map.genes[,1][ temp])   
}
colnames(res) <- rownames(res) <- name.gen
graph.res <- as(res,"graphNEL")
plot(graph.res)
