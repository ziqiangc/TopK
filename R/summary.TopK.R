#' Summarizing the TopK test results
#'
#' \code{summary} method for class \code{"lm"}.
#'
#' @param x \code{TopK} object, typically result of TopK.test
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
summary.TopK <- function(object, ...) {
    z <- object
    pval <- z$p.values
    Kval <- z$K.values
    mthd <- z$method
    tobs <- z$Tobs
    DNAME <- z$data.name

    cat("\n")
    cat(strwrap("TopK exact test", prefix = "\t"), sep = "\n")
    cat("\n")
    cat("Data:  ", DNAME, "\n", sep = "")
    cat("Method:  ", mthd, "\n", sep = "")
    cat(cat("K values:  "), cat(Kval,sep = ", "), "\n")
    cat("\n")
    cat("Statistics (unadjusted p-values):  \n")
    print(tobs)
    cat("\n")
    cat("TopK p-values:  \n")
    print(pval)
    cat("\n")

    invisible(object)

}
