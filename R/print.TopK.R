#' Print Methods for TopK test
#'
#' Pringting objects of class \code{"TopK"}.
#'
#'@param x \code{TopK} object, typically result of TopK.test.
#'@param ... further arguments to be passed to or from methods.
#'
#'@export

print.TopK <- function(x, ...) {
    cat("\n")
    cat(strwrap("TopK exact test", prefix = "\t"), sep = "\n")
    cat("\n")
    cat("Data:  ", x$data.name, "\n", sep = "")
    # cat("method:  ", x$method, "\n", sep = "")
    cat("\n")
    cat("TopK p-values:  \n")
    print(x$p.values)
    cat("\n")
    invisible(x)
}
