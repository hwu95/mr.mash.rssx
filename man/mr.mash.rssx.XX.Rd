\name{mr.mash.rssx.XX}
\alias{mr.mash.rssx.XX}
\title{Fit multivariate regression using summary statistics}
\usage{
mr.mash.rssx.XX(XX, XY, YY, n, B0=NULL, S0,
                         w0=rep(1/length(S0), length(S0)),
                         maxIter=50, tol=1e-4,
                         update_scheme=c("coordinate", "blockwise"),
                         acceleration=c("none", "squarem", "daarem"),
                         version=c("R", "Rcpp"),
                         daarem_order=10,
                         daarem_damping=0.5,
                         block_size=256,
                         block_tol=tol,
                         max_block_sweeps=100,
                         negative_gain_patience=2,
                         verbose=TRUE)
}
\arguments{
\item{XX}{List of R predictor cross-product matrices, each p by p.}
\item{XY}{List of R predictor-outcome cross-product vectors, each of length p.}
\item{YY}{Numeric vector of R outcome sums of squares.}
\item{n}{Numeric vector of R sample sizes.}
\item{B0}{Optional p by R initial coefficient matrix; defaults to zero.}
\item{S0}{List of K R by R prior covariance matrices, possibly including a zero matrix.}
\item{w0}{Initial nonnegative mixture weights; normalized to sum to one.}
\item{maxIter}{Maximum number of outer iterations.}
\item{tol}{Absolute ELBO convergence tolerance.}
\item{update_scheme}{Update scheme: coordinate or blockwise.}
\item{acceleration}{Mixture-weight acceleration: none, squarem, or daarem.}
\item{version}{Implementation: R or Rcpp.}
\item{daarem_order}{Positive integer history length for DAAREM.}
\item{daarem_damping}{DAAREM damping in (0, 1].}
\item{block_size}{Number of predictors per block.}
\item{block_tol}{Within-block gain tolerance; defaults to tol.}
\item{max_block_sweeps}{Maximum number of sweeps per block.}
\item{negative_gain_patience}{Number of consecutive negative block gains allowed before deactivation.}
\item{verbose}{Print iteration progress.}
}
\description{Fits mixture-of-normal-priors regression with coordinate or blockwise updates. Predictors require positive cross-product diagonals. Accelerated R and Rcpp trajectories may differ slightly due to floating-point arithmetic.}
\value{A list with posterior means B, residual variances v, mixture weights w0, posterior weights w1, covariance matrices Sj, ELBO and its trace, iteration count, and algorithm diagnostics.}
