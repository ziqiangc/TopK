#' TopK exact test based on permutation T p-values
#'
#' This function perform TopK exact test on a two group data set (Case-Control study).
#' The univariate p-values based on permutation t-test.
#'
#' @param X a data matrix containing features as rows and samples as columns.
#' @param nA number of samples in group 1
#' @param nB number of samples in group 2
#' @param Kvals a numeric vector indicating how many \code{"K"} we choose.
#' @param B the number of inner permutation, default is "0" for exact test.
#' @param alternative a character string specifying the alternative hypothesis,
#'          must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
#'          You can specify just the initial letter.
#' @param ReturnType character(1) specify how to return the results.
#'          Could be \code{"list"}, \code{"vector"} or rich format \code{"TopK"}. See detail.
#' @param vrb a logical vector indicating TRUE when show some words while running, FALSE otherwise.
#'
#' @examples TopK_PermT.test(X,7,7)
#'
#' @importFrom broman perm.test
#'
#' @export

TopK_PermT.test=function(X,nA,nB,
                       Kvals=c(1,2,3,4,5,10,25),
                       B=0,# set B=0 for exact test
                       alternative = c("two.sided", "less", "greater"),
                       ReturnType="TopK",vrb=T){
    # hrep=1; SimType="null1"; vrb=T; Kvals=c(1,2,3,4,5,10,25);B=0
    # alternative <- match.arg(alternative)
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

    # Put the functions outside the parents function

    if(vrb) cat("   Getting permutation t values.",fill=T)

    permScan <- function(i) perm.test(X[i, 1:nA], X[i, (nA+1):(nA+nB)], alternative = alternative)
    ExactWRS=t(mapply(permScan, 1:nrow(X)))


    # OPTION 2

    # allp <- array(NA, c(nrow(X),nPerm))
    # for (h in 1:nPerm) {
    #     for (i in 1:nrow(X)) {
    #         idx <- PermMat[,h]
    #         allp[i,h] <- t.test(X[i,idx], X[i,setdiff(1:N, idx)], var.equal = T)$p.value
    #     }
    # }

    #--- The code below runs 40% slower
    # tScan <- function(i,idx) t.test(X[i,idx], X[i,setdiff(1:N, idx)], var.equal = T)$p.value
    # tScan_h=function(h) mapply(tScan, 1:nrow(X),MoreArgs=list(idx=PermMat[,h]))
    #
    # res=mapply(tScan_h,1:nPerm)



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
    #  format output
    #---------------------------------------------------------------

    #  p.values
    p.values=c(Tadj, mean(KminP<=min(Tadj)))
    names(p.values)=c(paste("top",Kvals,sep=""),"minP.allk")

    # Tobs
    getTobs = function(i) TopKcdf[[i]][1]
    Tobs = mapply(getTobs, 1:length(Kvals))
    names(Tobs)=c(paste("top",Kvals,sep=""))
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

    if(ReturnType=="list") return(list(p.values=p.values,Tobs=Tobs,K.counts=K.counts,K.values=Kvals,
                                       N=N,nA=nA,nB=nB,nPerm=nPerm)
    )



    if(ReturnType=="vector"){
        return(c(p.values,Tobs,Kvals,K.counts[1,],K.counts[2,],K.counts[3,]))
    }

    if(ReturnType=="TopK") {
        res <- list(p.values=p.values,
                    Tobs=Tobs,
                    K.counts=K.counts,
                    TopK = TopK,
                    TopKcdf = TopKcdf,
                    Tadj = Tadj,
                    K.values=Kvals,N=N,nA=nA,nB=nB,nPerm=nPerm,
                    method = "Permutation t-test"
        )
        class(res) <- c(class(res), "TopK")
        return(res)
    }

}



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



perm.test <- function (x, y, alternative = c("two.sided", "less", "greater"),
                       var.equal = TRUE, pval = TRUE)
{
    # alternative <- match.arg(alternative)
    kx <- length(x)
    ky <- length(y)
    n <- kx + ky
    X <- c(x, y)
    # z <- rep(1:0, c(kx, ky))

    if(pval){
        pobs <- t.test(x, y, var.equal = var.equal, alternative=alternative)$p.value

        o <- binary.v(n)
        o <- o[, apply(o, 2, sum) == kx]
        nc <- choose(n, kx)
        allp <- 1:nc
        for (i in 1:nc) {
            xn <- X[o[, i] == 1]
            yn <- X[o[, i] == 0]
            allp[i] <- t.test(xn, yn, var.equal = var.equal, alternative = alternative)$p.value
        }

        return(allp)

    } else {

        tobs <- t.test(x, y, var.equal = var.equal, alternative = alternative)$statistic

        o <- binary.v(n)
        o <- o[, apply(o, 2, sum) == kx]
        nc <- choose(n, kx)
        allt <- 1:nc
        for (i in 1:nc) {
            xn <- X[o[, i] == 1]
            yn <- X[o[, i] == 0]
            allt[i] <- t.test(xn, yn, var.equal = var.equal, alternative = alternative)$statistic
        }

        return(allt)
    }
}


