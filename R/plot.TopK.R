#' @export
plot.TopK <- function(x, WobsDX=1, ...) {
    stopifnot(any(class(x))=="TopK")
    stopifnot(any(names(x))=="TopK")
    stopifnot(is.list(x$TopK))
    stopifnot(any(x$K.values)==2)

    TopK = x$TopK
    Top2=TopK[[2]]
    lgTop2=-log10(Top2)
    WobsDX=WobsDX

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
                cat(x1,x2,tcnt,Px1x2,fill=T)

                if((tcnt>0)&(Px1x2<=CumPobs)){

                    if((x1>=lgTop2[1,WobsDX])&(x2>=lgTop2[2,WobsDX])){
                        text(x1,x2,"[ ]",cex=2,col=4)
                    }
                    else{text(x1,x2,"[ ]",cex=2,col="purple")}

                    ttlTobstcnt=ttlTobstcnt+tcnt
                }
            }}}

}
