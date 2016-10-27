#' Create WRS Object
#'
#' This function generates WRS object to feed in later test function
#'
#' @param nA the number of sample in the first group
#' @param nB the number of sample in the secend group
#' @param inclPermMat a logical vector indicating TRUE when include a permutation matrix, FALSE otherwise.
#'
#' @examples CreateWRSobj(7,7)
#'
#' @export


CreateWRSobj=function(nA,nB,inclPermMat=T){
  # nA=10; nB=10

  N=nA+nB
  WRS.lwr=sum(1:nA)
  WRS.upr=sum((N-nA+1):N)
  WRS=rep(NA,WRS.upr)

  # brute force it, because we assume that nA and nB small enough that we can..
  PermMat=combn(N,nA)
  PermWRS=colSums(PermMat)
  for(i in WRS.lwr:WRS.upr){
    idx=which(PermWRS==i)[1]
    WRS[i]=WRSpvals=wilcox.test(PermMat[,idx],setdiff(1:N,PermMat[,idx]))$p.value
  }

  WRSobj=list(WRS=WRS,N=N,nA=nA,nB=nB)
  if(inclPermMat) WRSobj$PermMat=PermMat

  return(WRSobj)
}



