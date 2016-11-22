% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TopK.test.R
\name{TopK.test}
\alias{TopK.test}
\title{TopK exact test (main function)}
\usage{
TopK.test(x, g, Kvals = c(1, 2, 3, 4, 5, 10, 25), method = c("WRS",
  "t.test"), B = 0, alternative = c("two.sided", "less", "greater"),
  ReturnType = "TopK", vrb = T)
}
\arguments{
\item{g}{a vector or factor object giving the group for the corresponding samples of x.}

\item{Kvals}{a numeric vector indicating how many "K" we choose.}

\item{method}{a character string specifying which method used in TopK test,
must be one of \code{"WRS"} (default), "t.test". See detail.}

\item{B}{the number of inner permutation, default is "0" for exact test.}

\item{alternative}{a character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
You can specify just the initial letter.}

\item{ReturnType}{character(1) specify how to return the results.
Could be \code{"list"}, \code{"vector"} or rich format \code{"TopK"}. See detail.}

\item{vrb}{a logical vector indicating TRUE when show some words while running, FALSE otherwise.}

\item{X}{a data matrix containing features as rows and samples as columns.}
}
\description{
This function perform TopK exact test on a two group data set (Case-Control study).
The univariate p-values were calculated based on mulitiple choices of test.
}
\examples{
TopK.test(x,g,method="WRS")

}
