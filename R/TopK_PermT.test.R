#' TopK exact test based on permutation T p-values
#'
#' This function perform TopK exact test on a two group data set (Case-Control study).
#' The univariate p-values based on permutation t-test.
#'
#' @param X a data matrix containing features as rows and samples as columns.
#' @param nA number of samples in group 1
#' @param nB number of samples in group 2
#' @param Kvals a numeric vector indicating how many "K" we choose.
#' @param B the number of inner permutation, default is "0" for exact test.
#' @param ReturnType character(1) specify how to return the results.
#'          Could be "list" or "vector". See detail.
#' @param vrb a logical vector indicating TRUE when show some words while running, FALSE otherwise.
#'
#' @examples TopK_PermT.test(X,WRSobj)
#'
#' @importFrom broman perm.test
#'
#' @export

TopK_PermT.test=function(X,nA,nB,
                       Kvals=c(1,2,3,4,5,10,25),
                       B=0,# set B=0 for exact test
                       ReturnType="list",vrb=T){
    # hrep=1; SimType="null1"; vrb=T; Kvals=c(1,2,3,4,5,10,25);B=0

    N=nA+nB
    PermMat=combn(N,nA)
    nPerm=dim(PermMat)[2]

    # cat("nPerm=",nPerm,fill=T)

    nK=length(Kvals)

    #---------------------------------------------------------------
    # Get the exact permutation T p-values under all permutations
    #---------------------------------------------------------------

    # OPTION 1
    # adapt the intermediate function binary.v in package 'broman'
    binary.v <- function (n)
    {
        x <- 1:(2^n)
        mx <- max(x)
        digits <- floor(log2(mx))
        ans <- 0:(digits - 1)
        lx <- length(x)
        x <- matrix(rep(x, rep(digits, lx)), ncol = lx)
        (x%/%2^ans)%%2
    }


    perm.test <- function (x, y, var.equal = TRUE, pval = TRUE)
    {
        kx <- length(x)
        ky <- length(y)
        n <- kx + ky
        X <- c(x, y)
        # z <- rep(1:0, c(kx, ky))

        if(pval){
            pobs <- t.test(x, y, var.equal = var.equal)$p.value

            o <- binary.v(n)
            o <- o[, apply(o, 2, sum) == kx]
            nc <- choose(n, kx)
            allp <- 1:nc
            for (i in 1:nc) {
                xn <- X[o[, i] == 1]
                yn <- X[o[, i] == 0]
                allp[i] <- t.test(xn, yn, var.equal = var.equal)$p.value
            }

            return(allp)

        } else {

            tobs <- t.test(x, y, var.equal = var.equal)$statistic

            o <- binary.v(n)
            o <- o[, apply(o, 2, sum) == kx]
            nc <- choose(n, kx)
            allt <- 1:nc
            for (i in 1:nc) {
                xn <- X[o[, i] == 1]
                yn <- X[o[, i] == 0]
                allt[i] <- t.test(xn, yn, var.equal = var.equal)$statistic
            }

            return(allt)
        }
    }

    if(vrb) cat("   Getting permutation t values.",fill=T)
    ExactWRS=t(mapply(permScan, 1:nrow(X)))


    # OPTION 2
    # t.test(xn, yn, var.equal=TRUE)$p.value


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
    p.values=c(Tadj, mean(KminP<=min(Tadj)),
               BPO.pvalues[1],mean(KminP.BPO<=KminP.BPO[1]))
    names(p.values)=c(paste("top",Kvals,sep=""),"minP.allk","BPO.pvalue","minP.allk.BPO")

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
                                       N=N,nA=nA,nB=nB,nPerm=nPerm,
                                       BPO.pvalue=BPO.pvalues[1])
    )



    if(ReturnType=="vector"){
        return(c(p.values,Kvals,K.counts[1,],K.counts[2,],K.counts[3,]))
    }

}


