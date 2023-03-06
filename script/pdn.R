#usePackage <- function(p) {
  #if (!is.element(p, installed.packages()[,1]))
  #  install.packages(p, dep = TRUE)
 # require(p, character.only = TRUE)
#}


#usePackage("igraph")
#usePackage('foreach')
#àusePackage('doParallel')
#usePackage('doMC')
#registerDoMC(detectCores()-1)
#usePackage('gbm')
#usePackage('rpart')
#usePackage('rgexf')
#usePackage('reshape2')
#usePackage('rCharts')


learnPDN = function(pdndataframe, families = "poisson", method = "gbm",
                    evidenceMask=rep(0,ncol(pdndataframe))){
  
  #if(sum(apply(pdndataframe,2,var)==0) > 0){
  #  stop("error in input data, a column has no variance")
  #}
  
  cols = colnames(pdndataframe)
  variables = cols[1:ncol(pdndataframe)*!evidenceMask]
  evidence  = cols[1:ncol(pdndataframe)*evidenceMask]
  
  graph = graph.empty(n=0, directed=TRUE)
  graph <- graph + vertices(unlist(cols))
  V(graph)$evidence = TRUE
  
  models = foreach(i = 1:length(variables)) %do% {
      
    if(length(families) > 1){
      family = families[i]
    }else{
      family = families
    }
        
    pdndata = pdndataframe[!is.na(pdndataframe[[variables[i]]]),]
    neighbornames = c(variables[-i],evidence)
    formula = as.formula(paste(paste(variables[i], paste(neighbornames, collapse="+"), sep="~"), "-1", sep=""))
    
    vertex_size = sum(pdndata[[variables[i]]])
    
    if(method == "glm"){
      if(family == "poisson"){
        fam = poisson(link = "log")
      }else{
        fam = family
      }
      
      #print(formula)
      
      model = glm(formula = formula, data = pdndata, family = fam)
      #weights = 1-drop1(model, test="LRT")[[5]][-1]
      weights = model$coefficients / sum(model$coefficients)
      
    }else if(method == "gbm"){
      #print(formula)
      #print(family)
      model = gbm(formula = formula, data= pdndata, distribution = family ,
                  n.trees = 20, bag.fraction = 1, n.minobsinnode = 1, keep.data=FALSE, 
                  interaction.depth = 1)

      weights = summary(model, plotit=FALSE, normalize=FALSE, order = FALSE)$rel.inf 
      weights = sqrt(weights)
    }else if(method == "rpart"){
      model = rpart(formula = formula, data = pdndata, method = family, control = rpart.control(cp = 0.001, maxdepth=30))
      #model = rpart(formula = formula, data = pdndata, method = family)
      model = prune(model, cp=   model$cptable[which.min(model$cptable[,"xerror"]),"CP"])
      
      weights = matrix(data = rep(0, length(neighbornames)), nrow = 1)
      colnames(weights) <- neighbornames
      weights[1,names(model$variable.importance)] = model$variable.importance
      
    }else if(method == "gbmo"){
      model = gbm(formula = formula, data= pdndata, distribution = family , n.trees = 100, 
                  bag.fraction = 1, n.minobsinnode = 2)
      #s = summary(model, plotit=FALSE, normalize=TRUE, method=permutation.test.gbm)
      #weights =  s$rel.inf[match(neighbornames, s$var)]
      weights = summary(model, plotit=FALSE, normalize=FALSE, order = FALSE, method=permutation.test.gbm)$rel.inf 
    } else if(method == "mi"){
      model = mutinformation(pdn$pdndata)[i,]
      weights = model[-i]
    }else if(method == "rc"){
      model = hoeffd(as.matrix(pdn$pdndata))
      weights = (model$D[i,]*(model$P[i, ] < 0.05)+0.5)[-i]
    }
    
    return(list(model = model, weights = weights, variable = variables[i], vertex_size = vertex_size, formula = formula))
    
  }
  
  for(i in 1:length(variables)){
    
    weights = models[[i]]$weights
    neighbornames = variables[-i]
    
    for(j in 1:length(neighbornames)){
      if(abs(weights[j]) > 0){
        
        graph <- graph + edges(c(neighbornames[[j]], variables[[i]]), weight=abs(weights[j]))
        # the other order seems to be wrong, hence we reverse the edge direction
        #graph <- graph + edges(c( variables[[i]], neighbornames[[j]]),weight=abs(weights[j]))
        
        if(method == "mi" || method == "rc"){
          graph <- graph + edges(c( variables[[j]], neighbornames[[i]]),weight=abs(weights[j]))
        }
        
      }
    }
    
    vertex_size = models[[i]]$vertex_size
    #graph = set.vertex.attribute(graph, "size", index=c(i), 1)
    #graph = set.vertex.attribute(graph, "size", index=c(i), vertex_size)
    V(graph)[name == variables[i]]$evidence = FALSE
    V(graph)[name == variables[i]]$size = vertex_size
    

  }
  
  weights = E(graph)$weight
  weights = abs(weights)
  weights = linnorm(weights, 0, 1)
  weights = round(weights,2)
  E(graph)$weight = weights

  V(graph)$label = V(graph)$name
  
  adjacencyMatrix = as.data.frame(as.matrix(get.adjacency(graph, attr="weight")))
  
  return(list("models" = models, "graph" = graph, "method" = method, "families" = families, "evidenceMask" = evidenceMask, "adjacencyMatrix" = adjacencyMatrix))
}


pdn.gibbs.samples = function(pdn, evidence, evidenceMask=rep(0,ncol(pdndataframe)), iter=100, burnin = 10){
  last = evidence
  samples = evidence
  
  variables = colnames(evidence)[evidenceMask == 0]
  
  models <- new.env(hash=T, parent=emptyenv())
  for(m in pdn$models){
    models[[m$variable]] <- m
  }
  
  for(i in 1:iter){
    for(v in variables){
      
      m = models[[v]]
      
      lambda = predict(m$model, last, n.trees = 100)
      lambda = abs(lambda)
      
      prediction = rpois(1, lambda)
      
      sample = last
      sample[[v]] = prediction
      
      samples = rbind(samples, sample)
      
      last = sample
    }
  }
  
  samples = samples[burnin:iter,]
  rownames(samples)<-NULL
  return(samples)
}


plot.am = function(graph){
  corrmatrix = get.adjacency(graph, attr="weight")
  corrdata=as.data.frame(as.matrix(corrmatrix))
  corrdata$Variable1=names(corrdata)
  corrdatamelt=melt(corrdata,id="Variable1")
  names(corrdatamelt)=c("Variable1","Variable2","CorrelationCoefficient")
  corrmatplot = rPlot(Variable2 ~ Variable1, color = 'CorrelationCoefficient', data = corrdatamelt, type = 'tile', height = 600)
  corrmatplot$addParams(height = 400, width=800)
  corrmatplot$guides("{color: {scale: {type: gradient2, lower: 'red', middle: 'white', upper: 'blue',midpoint: 0}}}")
  corrmatplot$guides(y = list(numticks = length(unique(corrdatamelt$Variable1))))
  corrmatplot$guides(x = list(numticks = 3))
  corrmatplot
}

plot.pdn.contingency = function(graph, accuracy, score){
  
  l <- t(as.matrix(rbind(score,accuracy)))
  # Define edge widths:
  E(graph)$width <- E(graph)$weight*5+1
  #V(graph)$label <-""
  E(graph)$arrow.width <- 1
  
  E(graph)$curved <- 0.2
  
  
  plot.igraph(graph,layout=l,vertex.color=rainbow(length(V(graph))),
              vertex.size=2,rescale=FALSE,asp=0,xlim=c(min(c)-0.1,max(c)+0.1),
              vertex.label = "", 
              ylim=c(min(acc)-0.1,max(acc)+0.1),axes=TRUE,xlab="social score",
              ylab="accuracy score")
  textxy(c,acc,V(graph)$name,cex=1.2)
  grid(5, 5, lwd = 2)
  
  
}

plot.pdn = function(graph, plotWeights=TRUE){
  
  
  weights = E(graph)$weight
  weights = abs(weights)
  weights = linnorm(weights, 0.5, 5)
  weights = round(weights,2)
  
  #exp(2*weights+0.00001)
  
  #weights = normalize(weights, method="range")
  set.seed(107)
  
  vertex_size = V(graph)$size
  vertex_size = linnorm(vertex_size, 8,20)
  #vertex_size = 1-(max(vertex_size)-vertex_size)/(max(vertex_size)-min(vertex_size))
  
  par(bg='white') 
  if(plotWeights == TRUE){
    
    plot.igraph(graph,vertex.label=V(graph)$name,layout=layout.fruchterman.reingold(graph, niter=10000, area=50*vcount(graph)^2),  vertex.size = vertex_size, vertex.label.color="black",edge.color="gray", edge.width=weights, edge.label=weights, edge.curved=TRUE, edge.arrow.size=0.25, edge.arrow.width=2, edge.curved=0.15, edge.label.cex=0.8)
    #plot.igraph(graph,vertex.label=V(graph)$name,layout=layout.circle,  vertex.size = vertex_size, vertex.label.color="black",edge.color="gray",edge.width=weights, edge.arrow.size=1,edge.curved=TRUE, edge.label=round(E(graph)$weight,2),edge.arrow.size=0.25,edge.arrow.width=0.25, edge.curved=0.0, edge.label.cex=0.5)
  }else{
    plot.igraph(graph,vertex.label=V(graph)$name,layout=layout.circle, vertex.size = vertex_size, vertex.label.color="black",edge.color="gray",edge.width=exp(2*(weights)), edge.arrow.size=1,edge.curved=TRUE)
  }
}

linnorm = function(x, lower, upper){
  
  i = min(x)
  a = max(x)
  
  if(a == i)
    return(x)
  
  num = upper*(x - i) + lower*(a - x)
  den = (a-i)
  
  return(num/den)
  
}


saveAsGEXF = function(g, filepath="converted_graph.gexf")
{
  if(is.null(V(g)$label))
    V(g)$label <- as.character(V(g)$name)
  
  if(is.null(E(g)$weight))
    E(g)$weight <- rep.int(1, ecount(g))
  
  nodes <- data.frame(cbind(V(g), V(g)$label))
  edges <- t(Vectorize(get.edge, vectorize.args='id')(g, 1:ecount(g)))
  
  vAttrNames <- setdiff(list.vertex.attributes(g), "label")
  nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&",get.vertex.attribute(g, attr))), stringsAsFactors = FALSE)
  
  eAttrNames <- setdiff(list.edge.attributes(g), "weight")
  edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&",get.edge.attribute(g, attr))), stringsAsFactors = FALSE)
  
  output <- write.gexf(nodes, edges, edgesWeight=E(g)$weight, edgesAtt = edgesAtt, nodesAtt = nodesAtt)
  print(output, filepath, replace=T)
}

