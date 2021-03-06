% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TopK_WRS.test.R
\name{TopK_WRS.test}
\alias{TopK_WRS.test}
\title{TopK exact test based on WRS p-values}
\usage{
TopK_WRS.test(X, nA, nB, Kvals = c(1, 2, 3, 4, 5, 10, 25), B = 0,
  alternative = c("two.sided", "less", "greater"), ties.method = c("random",
  "min", "max", "average"), ReturnType = "TopK", vrb = T)
}
\arguments{
\item{X}{a data matrix containing features as rows and samples as columns.}

\item{nA}{number of samples in group 1}

\item{nB}{number of samples in group 2}

\item{Kvals}{a numeric vector indicating how many \code{"K"} we choose.}

\item{B}{the number of inner permutation, default is "0" for exact test.}

\item{alternative}{a character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
You can specify just the initial letter.}

\item{ties.method}{a character string specifying how ties are treated, see ‘Details’; can be abbreviated.}

\item{ReturnType}{character(1) specify how to return the results.
Could be \code{"list"}, \code{"vector"} or rich format \code{"TopK"}. See detail.}

\item{vrb}{a logical vector indicating TRUE when show some words while running, FALSE otherwise.}
}
\description{
This function perform TopK exact test on a two group data set (Case-Control study).
The univariate p-values based on Wilcoxon rank sum test.
}
\examples{
TopK_WRS.test(X,7,7)

}

