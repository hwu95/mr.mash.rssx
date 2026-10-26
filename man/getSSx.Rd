\name{getSSx}
\alias{getSSx}
\title{Construct sufficient statistics from GWAS summary statistics}
\usage{
getSSx.test(Bhat, Shat, LD, YY, n)
}
\arguments{
\item{Bhat}{p by R coefficient matrix.}
\item{Shat}{p by R matrix of positive standard errors.}
\item{LD}{A p by p LD correlation matrix or list of R such matrices.}
\item{YY}{Outcome sums of squares.}
\item{n}{Sample sizes, each greater than two.}
}
\description{Converts GWAS summary statistics and LD into cross-products in the same predictor order.}
\value{A list containing XX, XY, YY and n.}
