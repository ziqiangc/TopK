#' TopK exact test (main function)
#'
#' This function perform TopK exact test on a two group data set (Case-Control study).
#' The univariate p-values were calculated based on mulitiple choices of test.
#'
#' @param x a data matrix containing features as rows and samples as columns.
#' @param g a vector or factor object giving the group for the corresponding samples of x.
#' @param Kvals a numeric vector indicating how many \code{"K"} we choose.
#' @param method a character string specifying which method used in TopK test,
#'          must be one of \code{"WRS"} (default), "t.test". See detail.
#' @param B the number of inner permutation, default is "0" for exact test.
#' @param alternative a character string specifying the alternative hypothesis,
#'          must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
#'          You can specify just the initial letter.
#' @param ReturnType character(1) specify how to return the results.
#'          Could be \code{"list"}, \code{"vector"} or rich format \code{"TopK"}. See detail.
#' @param vrb a logical vector indicating TRUE when show some words while running, FALSE otherwise.
#'
#' @examples TopK.test(x,g,method="WRS")
#'
#' @import MASS
#'
#' @export
#'
TopK.test <- function(x, g,
                      Kvals=c(1,2,3,4,5,10,25),
                      method = c("WRS", "t.test"),
                      B=0,# set B=0 for exact test
                      alternative = c("two.sided", "less", "greater"),
                      ReturnType="TopK",vrb=T) {

    method <- match.arg(method)
    alternative <- match.arg(alternative)
    if (!is.numeric(x))
        stop("'x' must be numeric")

    if (missing(g))
        stop("'g' cannot be missing")

    if (dim(x)[2] != length(g))
        stop("'x' and 'g' must have the same length in terms of the number of samples")

    if (!is.numeric(Kvals))
        stop("'Kvals' must be numeric")

    if (max(Kvals)>=dim(x)[1])
        stop("The number of features must be greater than 'K' ")

    if (length(unique(g))<2)
        stop("The number of group must be greater than 2")

    if ((length(unique(g))>2) && (method == "WRS" || method == "t.test")  ) {
        warning("'WRS' and 't.test' are not available for the data have more than two groups, use other methods instead")
    }

    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))

    X = x [ ,order(g) ]



    if (method == "WRS") {
        RVAL <- TopK_WRS.test(X=X, nA = nA, nB = nB,
                             Kvals = Kvals,
                             B = B,
                             alternative = alternative,
                             ReturnType=ReturnType,vrb=vrb)
    }

    if (method == "t.test") {
        RVAL <- TopK_PermT.test(X=x, nA = nA, nB = nB,
                                Kvals = Kvals,
                                B = B,
                                alternative = alternative,
                                ReturnType=ReturnType,vrb=vrb)
    }

    RVAL$data.name = DNAME

    return(RVAL)


}



