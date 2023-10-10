pvbsrBootstrap = function(y,x,nsamp=100,cores=8){
  library(dplyr)
  n <- length(y)
  replicateMatrix = sample(1:(n*nsamp),replace=TRUE) %>%
    matrix(n,nsamp)
  replicateMatrix = replicateMatrix%%n +1
  #print(replicateMatrix[1:5,])
  fxn1 <- function(shuf,y,x){
    #library(utilityFunctions) return(utilityFunctions::fastlmbeta(y[shuf],x[shuf,])) library(vbsr) y=y[shuf] 
    #x=x[shuf,]
    res <- vbsr::vbsr(y[shuf],x[shuf,])
    
    ###identify significant features return unpenalized betas
    baz = rep(0,ncol(x)+1)
    names(baz) = c('intercept',colnames(x))
    whichSig = which(res$pval < 0.05/ncol(x))
    if(length(whichSig)>0){
      #write a fastlmbeta function that also returns correlation
      baz2 = c()
      try(baz2 <- spike::fastlmbeta2(y[shuf],x[shuf,whichSig],colnames(x)[whichSig]),silent=T)
      if(length(baz2)>0){
        baz[names(baz2)] <- baz2
      }
    }
    
    #return(c(res$alpha,res$beta))
    return(baz)
  }
  
  cl <- parallel::makeCluster(cores)
  registerDoParallel(cl)
  #betaMatrix <- foreach(i=1:nsamp,.combine='rbind') %dopar% fxn1(replicateMatrix[,i],y,data.matrix(x))
  betaMatrix <- t(parallel::parApply(cl,replicateMatrix,2,fxn1,y,data.matrix(x)))
  #betaMatrix <- t(apply(replicateMatrix,2,fxn1,y,data.matrix(x)))
  parallel::stopCluster(cl)
  #colnames(betaMatrix) <- c(colnames(x))
  return(apply(betaMatrix!=0,2,mean))
  #return(betaMatrix)
  
}
