#' Plot triangle figure for a TopK object
#'
#' A triangle plot illustrate the results of TopK test. This is a plot only for the scenario of \code{K=2}.
#' The x-axis denotes the \code{-log10 p_(1)} -value and the y-axis denotes \code{-log10 p_(2)} -value.
#' The p-values were tabulated across all permutations and displayed as counts (for WRS test);
#' as scattered and smoothed points (for permutation t-test).
#' The red circle denotes the value combination that was observed in the actual data.
#'
#'@param x \code{TopK} object, typically result of TopK.test
#'@param WobsDX the index of observed data in all permutations. Default value is \code{1}
#'@param ... other parameters to be passed through to plotting functions.


#' @export
plot.TopK <- function(x, WobsDX=1, ...) {
    stopifnot(any(class(x)=="TopK"))
    stopifnot(any(names(x)=="TopK"))
    stopifnot(is.list(x$TopK))
    stopifnot(any(x$K.values==2))

    TopK = x$TopK
    Top2=TopK[[2]]
    lgTop2=-log10(Top2)
    WobsDX=WobsDX


    if (x$method == "WRS") {
        plot(-log10(Top2[1,]),-log10(Top2[2,]),pch=" ",axes=F,
             main="Top-K FWER Example for K=2",
             xlab=expression(paste("-log"[10] , "p"[K:1])),
             ylab=expression(paste("-log"[10] , "p"[K:2])))
        axis(1); axis(2)

        ttltcnt=0
        for(x1 in unique(-log10(Top2[1,]))){
            for(x2 in unique(-log10(Top2[2,]))){
                tcnt=sum(((-log10(Top2[1,]))==x1)
                         &((-log10(Top2[2,]))==x2) )
                ttltcnt=ttltcnt+tcnt
                if(tcnt>0){
                    if(x1>=lgTop2[1,WobsDX]) text(x1,x2,tcnt,cex=0.8,col="forestgreen")
                    if(x1<lgTop2[1,WobsDX]) text(x1,x2,tcnt,cex=0.65,col="gray30")
                }
            }}
        text(-log10(Top2[1,WobsDX]),-log10(Top2[2,WobsDX]),"O",cex=3,col=2)

        #
        # graphical representation of what the final minP correction does
        #


        CumPobs=mean((Top2[1,]<=Top2[1,WobsDX])&(Top2[2,]<=Top2[2,WobsDX]))

        ttlTobstcnt=0.0
        for(x1 in unique(-log10(Top2[1,]))){
            for(x2 in unique(-log10(Top2[2,]))){
                if(x2<=x1){
                    Px1x2=mean((lgTop2[1,]>=x1)&(lgTop2[2,]>=x2))
                    tcnt=sum(((lgTop2[1,])==x1)
                             &((lgTop2[2,])==x2) )
                    #cat(x1,x2,tcnt,Px1x2,fill=T)

                    if((tcnt>0)&(Px1x2<=CumPobs)){

                        if((x1>=lgTop2[1,WobsDX])&(x2>=lgTop2[2,WobsDX])){
                            text(x1,x2,"[ ]",cex=2,col=4)
                        }
                        else{text(x1,x2,"[ ]",cex=2,col="purple")}

                        ttlTobstcnt=ttlTobstcnt+tcnt
                    }
                }}}

    }

    if(x$method == "PermT") {
        X1.grey = c()
        X2.grey = c()
        X1.green = c()
        X2.green = c()
        ttltcnt=0
        ulgtop2.1 = unique(-log10(Top2[1,]))
        ulgtop2.2 = unique(-log10(Top2[2,]))

        for(x1 in ulgtop2.1){
            for(x2 in ulgtop2.2){
                tcnt=sum((lgTop2[1,]==x1)
                         &(lgTop2[2,]==x2) )
                ttltcnt=ttltcnt+tcnt
                if(tcnt>0){
                    if(x1>=lgTop2[1,WobsDX]) {
                        X1.green = c(X1.green, x1)
                        X2.green = c(X2.green, x2)
                    }
                    if(x1<lgTop2[1,WobsDX]) {
                        X1.grey = c(X1.grey, x1)
                        X2.grey = c(X2.grey, x2)
                    }
                }
            }}

        X1.blue = c()
        X2.blue = c()
        X1.purple = c()
        X2.purple = c()


        CumPobs=mean((Top2[1,]<=Top2[1,WobsDX])&(Top2[2,]<=Top2[2,WobsDX]))

        ttlTobstcnt=0.0
        for(x1 in ulgtop2.1){
            for(x2 in ulgtop2.2){
                if(x2<=x1){
                    Px1x2=mean((lgTop2[1,]>=x1)&(lgTop2[2,]>=x2))
                    tcnt=sum(((lgTop2[1,])==x1)
                             &((lgTop2[2,])==x2) )
                    # cat(x1,x2,tcnt,Px1x2,fill=T)

                    if((tcnt>0)&(Px1x2<=CumPobs)){

                        if((x1>=lgTop2[1,WobsDX])&(x2>=lgTop2[2,WobsDX])){
                            X1.blue = c(X1.blue, x1); X2.blue = c(X2.blue, x2)
                        }
                        else{
                            X1.purple = c(X1.purple, x1); X2.purple = c(X2.purple, x2)
                        }

                        ttlTobstcnt=ttlTobstcnt+tcnt
                    }
                }}}

        plot(-log10(Top2[1,]),-log10(Top2[2,]),pch=" ",axes=F,
             main="Top-K FWER Example for K=2",
             xlab=expression(paste("-log"[10] , "p"[K:1])),
             ylab=expression(paste("-log"[10] , "p"[K:2])))
        axis(1); axis(2)
        smoothScatter(X1.grey, X2.grey, colramp=colorRampPalette(c(rgb(1, 1, 1, 0), "grey33"), alpha=TRUE), add=T)
        points(X1.green, X2.green, pch="+", col="green")
        points(X1.blue, X2.blue, pch="o", col="blue")
        points(X1.purple, X2.purple, pch=7, col="purple")
        text(-log10(Top2[1,WobsDX]),-log10(Top2[2,WobsDX]),"o",cex=3,col=2)


    }


}
