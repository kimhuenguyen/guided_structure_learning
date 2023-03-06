
result_algorithms <- function(n, p,Adj,Pmat, nsim, npa, alpha,cpu){
  W <- matrix(NA,nsim,77)
  set.seed(123)
  for (i in 1:nsim){
    data <- simData(n,p,order,Pmat)
    X <- data
    colnames(X) <- seq(1,p,1)
    
    # PK2 with BIC
    PKBC.time <- system.time(gPB <- try(Poisk2(X, order, criterion="BIC", npa),silent=TRUE))
    if (unique(class (gPB)!="try-error")){
      PKBC <- result(Adj,gPB$a)
    }else PKBC <- rep(NA,6)
    
    # LPGM with order
    LPGM.time <- system.time(LP <- try(LPGM.ord(X,order,alpha),silent=TRUE))
    if (unique(class (LP)!="try-error")){
      LPGM <- result(Adj,LP)
    }else LPGM <- rep(NA,6)
    
    # OR-PPGM
    PCLPGM.time <- system.time(PCP <- try(pois.ord (X,npa,alpha,order),silent=TRUE))
    if (unique(class (PCP)!="try-error")){
      PCLPGM <- result(Adj,PCP )
    }else PCLPGM <- rep(NA,6)
    
    # PDN
    Y <- data.frame(X)
    PDN.time <- system.time(fit <- try(learnPDN(Y, families = "poisson"),silent=TRUE))
    if (unique(class (fit)!="try-error")){
      beta <- as.matrix(fit$adjacencyMatrix)
      beta[beta != 0] <- 1
      PDN <- result(Adj, beta)
    }else PDN <- rep(NA,6)
    
    # ODS
    ODS.time <- system.time(  O <- try(ODS(X,cpu),silent=TRUE))
    if (unique(class (O)!="try-error")){
      ods <- result(Adj, O)
    } else ods <- rep(NA,6)   
    
    # MMHC with mixtured models
    dataD <- try(discMBC(log(1+X)),silent=TRUE)
    if (unique(class (dataD)!="try-error")){
      dataD <- data.frame(dataD)
      MMHC.time <- system.time(  try( fit <- mmhc(dataD),silent=TRUE))
      if (unique(class (fit)!="try-error")){
        adjM <- amat(fit)
        MMHC <- result(Adj, adjM)
      }else   MMHC <- rep(NA,6)  
    }else   MMHC <- rep(NA,6)

    # PC with log transform data
    dataL <- log(1+X)
    V <- colnames(dataL)
    PC.time <- system.time(pc.fit <- try(pc(suffStat = list(C = cor(dataL), n = n),
                                            indepTest = gaussCItest, ## indep.test: partial correlations
                                            alpha, labels = V),silent=TRUE))
    if (unique(class (pc.fit)!="try-error")){
      gpc <- pcalg::pdag2dag(pc.fit@graph, keepVstruct=TRUE)
      adj.PC <- as(gpc$graph, "matrix")
      PC <- result(Adj,adj.PC)
    }else PC <- rep(NA,6)
    
    
    W[i,] <-c(PKBC,PKBC.time,PCLPGM,PCLPGM.time,LPGM,LPGM.time,PDN,PDN.time,ods,ODS.time,PC,PC.time,MMHC,MMHC.time)
  }
  return(W)
}


