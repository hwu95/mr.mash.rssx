list2array.rssx=function(x) {
    nr=nrow(x[[1]])
    nc=ncol(x[[1]])
    out=array(0, dim=c(nr, nc, length(x)))

    for (i in seq_along(x)) {
        out[,,i]=x[[i]]
    }

    out
}

Sj2array.rssx=function(Sj, R) {
    out=array(0, dim=c(R, R, length(Sj)))

    for (j in seq_along(Sj)) {
        if (!is.null(Sj[[j]])) {
            out[,,j]=Sj[[j]]
        }
    }

    out
}

ELBO.rssx=function(v, ERSS, KL, n) {
    -0.5*sum(n*log(2*pi*v))-0.5*sum(ERSS/v)-sum(KL)
}

#' Construct sufficient statistics from GWAS summary statistics
#'
#' Converts GWAS summary statistics (coefficients and standard errors) and LD into
#' cross-products for mr.mash.rssx.XX. Sample sizes must exceed two,
#' standard errors must be positive, and all inputs must use the same
#' predictor order and compatible scaling.
#'
#' @param Bhat p by R matrix of coefficients.
#' @param Shat p by R matrix of standard errors.
#' @param LD A p by p LD correlation matrix, or a list of R such matrices.
#' @param YY Numeric vector of R outcome sums of squares.
#' @param n Numeric vector of R sample sizes.
#' @return A list with XX, XY, YY and n, suitable as input to mr.mash.rssx.XX.
#' @export
getSSx=function(Bhat, Shat, LD, YY, n) {
    p=nrow(Bhat)
    R=ncol(Bhat)

    stopifnot(all(dim(Shat)==c(p, R)), length(YY)==R, length(n)==R)

    if (is.matrix(LD)) {
        stopifnot(all(dim(LD)==c(p, p)))
        LD_list=rep(list(LD), R)
    } else if (is.list(LD)) {
        stopifnot(length(LD)==R)

        for (r in 1:R) {
            stopifnot(all(dim(LD[[r]])==c(p, p)))
        }

        LD_list=LD
    } else {
        stop("LD must be either a p x p matrix or a list of R p x p matrices")
    }

    XX=vector("list", R)
    XY=vector("list", R)

    for (r in 1:R) {
        vY=YY[r]/(n[r]-1)
        zhat=Bhat[,r]/Shat[,r]
        a=(n[r]-1)/(n[r]-2+zhat^2)
        dx=vY*a/(Shat[,r]^2)

        XX[[r]]=LD_list[[r]]*outer(sqrt(dx), sqrt(dx))
        XY[[r]]=a*zhat*vY/Shat[,r]
    }

    list(XX=XX,
         XY=XY,
         YY=YY,
         n=n)
}
