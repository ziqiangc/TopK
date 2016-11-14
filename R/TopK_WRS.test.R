#' TopK exact test based on WRS p-values
#'
#' This function perform TopK exact test on a two group data set (Case-Control study).
#' The univariate p-values based on Wilcoxon rank sum test.
#'
#' @param X a data matrix containing features as rows and samples as columns.
#' @param nA number of samples in group 1
#' @param nB number of samples in group 2
#' @param Kvals a numeric vector indicating how many "K" we choose.
#' @param B the number of inner permutation, default is "0" for exact test.
#' @param alternative a character string specifying the alternative hypothesis,
#'          must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
#'          You can specify just the initial letter.
#' @param ReturnType character(1) specify how to return the results.
#'          Could be \code{"list"}, \code{"vector"} or rich format \code{"TopK"}. See detail.
#' @param vrb a logical vector indicating TRUE when show some words while running, FALSE otherwise.
#'
#' @examples TopK_WRS.test(X,7,7)
#'
#' @import MASS
#'
#' @export

TopK_WRS.test=function(X,nA,nB,
                             Kvals=c(1,2,3,4,5,10,25),
                             B=0,# set B=0 for exact test
                             alternative = c("two.sided", "less", "greater"),
                             ReturnType="TopK",vrb=T){
    # hrep=1; SimType="null1"; vrb=T; Kvals=c(1,2,3,4,5,10,25);B=0
    # alternative <- match.arg(alternative)
    WRSobj = CreateWRSobj(nA, nB, alternative=alternative)

    # Grab important objects from WRSobj
    N=WRSobj$N
    nA=WRSobj$nA
    nB=WRSobj$nB
    WRS=WRSobj$WRS

    if(length(WRSobj$PermMat)>0){
        PermMat=WRSobj$PermMat
        if(B!=0){ # if B!=0 and PermMat provided then keep 1st perm (observed data)
            #  and add B-1 additional permutations - sample without replacement
            PermMat=PermMat[,c(1,sample(2:choose(N,nA),B-1))] }
    } # end conditional that  Permmat is present

    if(length(WRSobj$PermMat)==0){
        # randomly generate permutations
        PermMat=replicate(B-1,sample(1:N,nA))
        PermMat=cbind(1:nA,PermMat)
    } # end conditional that  Permmat is not present

    nPerm=dim(PermMat)[2]

    # cat("nPerm=",nPerm,fill=T)

    nK=length(Kvals)

    #---------------------------------------------------------------
    # Get the Exact WRS p-values under all permutations
    #---------------------------------------------------------------

    # populate rank matrix
    RX=t(apply(X,1,rank))

    if(vrb) cat("   Getting WRS values.",fill=T)

    hMap=function(h) WRS[rowSums(RX[,PermMat[,h]])]
    ExactWRS=mapply(hMap,1:nPerm)


    #---------------------------------------------------------------
    # Get the TopK  p-values under each permutation
    #---------------------------------------------------------------


    GetTopK=function(K){
        jgetK=function(j) ExactWRS[order(ExactWRS[,j])[1:K],j]
        if(K>1) return(mapply(jgetK,1:nPerm))
        if(K==1) return(matrix(mapply(jgetK,1:nPerm),nrow=1))
    }


    TopK=list()
    for(k in Kvals){
        TopK=c(TopK,list(GetTopK(k)))
    }


    #---------------------------------------------------------------
    # Get the eCDF estimate for each TopK value..
    #---------------------------------------------------------------

    # have to doubly correct the the top2
    GetTopKcdf=function(topK){
        K=dim(topK)[1]
        jCDFK=function(j) mean(colSums((topK-topK[,j])<=0)==K)
        return(mapply(jCDFK,1:nPerm))
    }

    TopKcdf=list()
    for(k in 1:nK){
        TopKcdf=c(TopKcdf,list(GetTopKcdf(TopK[[k]])))
    }


    #---------------------------------------------------------------
    # Adjust eCDF estimate for each TopK value, so that we have
    # a legitimate p-value
    #---------------------------------------------------------------

    getTadj=function(ik) mean(TopKcdf[[ik]]<=(TopKcdf[[ik]][1]))
    Tadj=mapply(getTadj,1:nK)



    #---------------------------------------------------------------
    # Get minP adjusted p-value across all considered
    # choices of K
    #---------------------------------------------------------------

    # Let's get the minP across them all

    CDFmat=t(matrix(unlist(TopKcdf),ncol=nK))

    ipval=function(pvec) rank(pvec,ties.method="max")
    Kpvals=apply(CDFmat,1,ipval)
    Kpvals=Kpvals/(dim(Kpvals)[1])
    KminP=apply(Kpvals,1,min)


    #---------------------------------------------------------------
    # Get the BPO p-values
    # BPO: best possible outcome
    # BPO p-value keeps track of the number of entries at the BPO
    # !! Have to think about ramifications of two-sided tests !!
    #---------------------------------------------------------------
    p.BPO=min(WRS,na.rm=T)

    jBPO=function(pvec) sum(pvec==p.BPO)
    prm.BPO=apply(ExactWRS,2,jBPO)
    BPO.pvalues= 1-rank(prm.BPO,ties="max")/length(prm.BPO)

    # consider a k dependent BPO adjustment:
    #  if # BPO > K then use min(BPO.pvalue, Kpval)

    KpvalsBPO=Kpvals

    for(k in 1:nK){
        swapDX=which(prm.BPO>Kvals[k])
        if(length(swapDX)>0){
            swapDX=swapDX[BPO.pvalues[swapDX]<Kpvals[swapDX,k]]
            KpvalsBPO[swapDX,k]=BPO.pvalues[swapDX]
        }
    }
    # plot(-log10(Kpvals),-log10(KpvalsBPO));abline(0,1,col=2)

    KminP.BPO=apply(KpvalsBPO,1,min)


    #---------------------------------------------------------------
    #  format output
    #---------------------------------------------------------------

    #  p.values
    p.values=c(Tadj, mean(KminP<=min(Tadj)))
    names(p.values)=c(paste("top",Kvals,sep=""),"minP.allk")

    # actual K
    # note: due to ties, can end up with more than K values tied at the top.
    # have some numerical precision / rounding issues!
    # as a work around, only keep the 8 most significant digits
    PrecExactWRS=round(ExactWRS[,1],8)
    unqP=sort(unique(PrecExactWRS))
    K.counts=summary(factor(PrecExactWRS,levels=unqP),maxsum=1e6)
    K.counts=rbind(K.counts,cumsum(K.counts))

    # now, match up the numbers of each of the K.cumsum values that are the alt model
    K.alt=numeric()
    for(ik in 1:length(unqP)) K.alt=c(K.alt,sum(Mdl[which(PrecExactWRS==unqP[ik])]))
    K.counts=rbind(K.counts,K.alt)
    K.counts=rbind(K.counts,cumsum(K.alt))


    if(vrb) print(round(p.values,4))

    if(ReturnType=="list") return(list(p.values=p.values,K.counts=K.counts,K.values=Kvals,
                                       N=N,nA=nA,nB=nB,nPerm=nPerm)
    )



    if(ReturnType=="vector"){
        return(c(p.values,Kvals,K.counts[1,],K.counts[2,],K.counts[3,]))
    }

    if(ReturnType=="TopK") {
        res <- list(p.values=p.values,
                    K.counts=K.counts,
                    TopK = TopK,
                    TopKcdf = TopKcdf,
                    Tadj = Tadj,
                    K.values=Kvals,N=N,nA=nA,nB=nB,nPerm=nPerm
                    )
        class(res) <- c(class(res), "TopK")
        return(res)
    }

}



CreateWRSobj=function(nA, nB, alternative= c("two.sided", "less", "greater"), inclPermMat=T){
    # nA=10; nB=10
    # alternative = match.arg(alternative)
    N=nA+nB
    WRS.lwr=sum(1:nA)
    WRS.upr=sum((N-nA+1):N)
    WRS=rep(NA,WRS.upr)

    # brute force it, because we assume that nA and nB small enough that we can..
    PermMat=combn(N,nA)
    PermWRS=colSums(PermMat)
    for(i in WRS.lwr:WRS.upr){
        idx=which(PermWRS==i)[1]
        WRS[i]=WRSpvals=wilcox.test(PermMat[,idx],setdiff(1:N,PermMat[,idx]), alternative=alternative)$p.value
    }

    WRSobj=list(WRS=WRS,N=N,nA=nA,nB=nB)
    if(inclPermMat) WRSobj$PermMat=PermMat

    return(WRSobj)
}




