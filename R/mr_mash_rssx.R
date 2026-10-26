#' Fit multivariate regression using summary statistics
#'
#' Fits a mixture-of-normal-priors regression with trait-specific residual
#' variances using coordinate or blockwise updates. Predictors must have
#' positive cross-product diagonals and prior covariances must be compatible
#' with the outcome dimension. Accelerated R and Rcpp trajectories may differ
#' slightly due to floating-point arithmetic.
#'
#' @param XX List of R predictor cross-product matrices, each p by p.
#' @param XY List of R predictor-outcome cross-product vectors, each of length p.
#' @param YY Numeric vector of R outcome sums of squares.
#' @param n Numeric vector of R sample sizes.
#' @param B0 Optional p by R initial coefficient matrix; defaults to zero.
#' @param S0 List of K R by R prior covariance matrices, possibly including a zero matrix.
#' @param w0 Initial nonnegative mixture weights; normalized to sum to one.
#' @param maxIter Maximum number of outer iterations.
#' @param tol Absolute ELBO convergence tolerance.
#' @param update_scheme Update scheme: "coordinate" or "blockwise".
#' @param acceleration Mixture-weight acceleration: "none", "squarem", or "daarem".
#' @param version Implementation: "R" or "Rcpp".
#' @param daarem_order Positive integer history length for DAAREM.
#' @param daarem_damping DAAREM damping in (0, 1].
#' @param block_size Number of predictors per block.
#' @param block_tol Within-block gain tolerance; defaults to tol.
#' @param max_block_sweeps Maximum number of sweeps per block.
#' @param negative_gain_patience Number of consecutive negative block gains allowed before deactivation.
#' @param verbose Print iteration progress.
#' @return A list containing posterior means B (p by R), residual variances v,
#'   mixture weights w0, posterior component weights w1, posterior covariance
#'   matrices Sj, ELBO, ELBO_trace, iterations, update_scheme, acceleration,
#'   and version. Additional fields record acceleration and block diagnostics.
#' @useDynLib mr.mash.rssx, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
mr.mash.rssx.XX=function(XX, XY, YY, n, B0=NULL, S0,
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
                         verbose=TRUE) {

    update_scheme=match.arg(update_scheme)
    acceleration=match.arg(acceleration)
    version=match.arg(version)

    if (length(daarem_order)!=1L | !is.finite(daarem_order) |
        daarem_order<1 | daarem_order!=as.integer(daarem_order)) {
        stop("daarem_order must be a positive integer")
    }

    if (length(daarem_damping)!=1L | !is.finite(daarem_damping) |
        daarem_damping<=0 | daarem_damping>1) {
        stop("daarem_damping must be in (0,1]")
    }

    daarem_order=as.integer(daarem_order)

    # Dimensions
    p=nrow(XX[[1]])
    R=length(XX)
    K=length(S0)

    if (version=="Rcpp") {
        S0_rcpp=list2array.rssx(S0)
        XY_rcpp=do.call(cbind, XY)
        XX_diag_rcpp=vapply(XX, diag, numeric(p))
    }

    # Initial posterior mean
    if (is.null(B0)) {
        B=matrix(0, nrow=p, ncol=R)
    } else {
        if (!is.matrix(B0) | !all(dim(B0)==c(p, R))) {
            stop("B0 must be a p x R matrix")
        }
        if (any(!is.finite(B0))) {
            stop("All elements of B0 must be finite")
        }
        B=B0
    }

    # Initial mixture weights
    if (any(!is.finite(w0)) | any(w0<0) | sum(w0)<=0) {
        stop("w0 must be nonnegative, finite, and have positive sum")
    }
    w0=w0/sum(w0)

    # Initial residual variances
    v=YY/n

    #### Coordinate algorithm
    if (update_scheme=="coordinate") {
        ELBO=numeric(maxIter)
        n_updates=0L

        if (acceleration!="none") {
            w0_residual_trace=numeric(maxIter)
            acceleration_accepted_trace=logical(maxIter)
        }

        if (acceleration=="squarem") {
            squarem_state=list(phase=0L, w00=NULL, w01=NULL)
            alpha_trace=rep(NA_real_, maxIter)
        }

        if (acceleration=="daarem") {
            daarem_state=list(w0=list(), residual=list(), since_restart=0L)
            daarem_order_trace=integer(maxIter)
            daarem_damping_trace=numeric(maxIter)
            daarem_restart_trace=logical(maxIter)
        }

        for (iter in 1:maxIter) {
            timeIn=proc.time()[3]

            w0_old=w0

            if (version=="Rcpp") {
                tmp=update_map_coordinate_rcpp(B=B, w0=w0, v=v, XX=XX, XX_diag=XX_diag_rcpp,
                                               XY=XY_rcpp, YY=YY, n=n, S0=S0_rcpp)
            } else {
                tmp=update_map_coordinate.rssx(B=B, w0=w0, v=v, XX=XX, XY=XY, YY=YY, n=n, S0=S0)
            }

            n_updates=n_updates+1L

            B=tmp$B
            W=tmp$W
            Sj=tmp$Sj
            KL=tmp$KL
            ERSS=tmp$ERSS

            w0_raw=tmp$w0
            v=tmp$v

            w0_residual=max(abs(w0_raw-w0_old))

            #### No acceleration
            if (acceleration=="none") {
                w0=w0_raw

            #### SQUAREM
            } else if (acceleration=="squarem") {
                tmp_acc=squarem_w0.rssx(w0=w0_old, w0_raw=w0_raw, state=squarem_state)

                w0=tmp_acc$w0
                squarem_state=tmp_acc$state
                alpha_trace[iter]=tmp_acc$alpha
                acceleration_accepted_trace[iter]=tmp_acc$accelerated

            #### DAAREM
            } else if (acceleration=="daarem") {
                tmp_acc=daarem_w0.rssx(w0=w0_old, w0_raw=w0_raw, state=daarem_state,
                                       order=daarem_order, damping=daarem_damping)

                w0=tmp_acc$w0
                daarem_state=tmp_acc$state
                daarem_order_trace[iter]=tmp_acc$order
                daarem_damping_trace[iter]=tmp_acc$damping
                daarem_restart_trace[iter]=tmp_acc$restarted
                acceleration_accepted_trace[iter]=tmp_acc$accelerated
            }

            if (acceleration!="none") {
                w0_residual_trace[iter]=w0_residual
            }

            ELBO[iter]=ELBO.rssx(v=v, ERSS=ERSS, KL=KL, n=n)

            if (iter==1) {
                delta_ELBO=NA_real_
            } else {
                delta_ELBO=ELBO[iter]-ELBO[iter-1]
            }

            timeOut=proc.time()[3]

            if (verbose) {
                cat("iteration =", iter,
                    "| ELBO =", round(ELBO[iter], 3),
                    "| delta =", round(delta_ELBO, 3),
                    "| time/iter =", round(timeOut-timeIn, 3), "s",
                    "\n")
            }

            if (iter>1 & abs(delta_ELBO)<tol) {break}
        }

        ELBO=ELBO[seq_len(iter)]

        res=list(B=B,
                 v=v,
                 w0=w0,
                 w1=W,
                 Sj=Sj,
                 ELBO=ELBO[length(ELBO)],
                 ELBO_trace=ELBO,
                 iterations=iter,
                 n_updates=n_updates,
                 update_scheme="coordinate",
                 acceleration=acceleration,
                 version=version)

        if (acceleration!="none") {
            res$w0_residual_trace=w0_residual_trace[seq_len(iter)]
            res$acceleration_accepted_trace=acceleration_accepted_trace[seq_len(iter)]
        }

        if (acceleration=="squarem") {
            res$alpha_trace=alpha_trace[seq_len(iter)]
        }

        if (acceleration=="daarem") {
            res$daarem_order_trace=daarem_order_trace[seq_len(iter)]
            res$daarem_damping_trace=daarem_damping_trace[seq_len(iter)]
            res$daarem_restart_trace=daarem_restart_trace[seq_len(iter)]
        }

    #### Blockwise algorithm
    } else if (update_scheme=="blockwise") {
        XX_diag=sapply(XX, diag)

        # SNP blocks
        block_id=ceiling(seq_len(p)/block_size)
        blocks=split(seq_len(p), block_id)
        n_blocks=length(blocks)

        if (length(negative_gain_patience)!=1L | !is.finite(negative_gain_patience) |
            negative_gain_patience<1 | negative_gain_patience!=as.integer(negative_gain_patience)) {
            stop("negative_gain_patience must be a positive integer")
        }

        negative_gain_patience=as.integer(negative_gain_patience)

        # Posterior storage
        Sj=vector("list", p)
        W=matrix(0, nrow=p, ncol=K)
        KL=numeric(p)

        # Initial ERSS
        ERSS=numeric(R)

        for (r in 1:R) {
            br=B[,r]
            ERSS[r]=YY[r]-2*sum(br*XY[[r]])+drop(t(br)%*%XX[[r]]%*%br)
        }

        ELBO=numeric(maxIter)

        # Block behavior
        block_gain_trace=matrix(NA_real_, nrow=n_blocks, ncol=maxIter)
        active_blocks_trace=integer(maxIter)
        negative_gain_streak=integer(n_blocks)
        block_evaluations_trace=matrix(0L, nrow=n_blocks, ncol=maxIter)

        if (acceleration!="none") {
            w0_residual_trace=numeric(maxIter)
            acceleration_accepted_trace=logical(maxIter)
        }

        if (acceleration=="squarem") {
            squarem_state=list(phase=0L, w00=NULL, w01=NULL)
            alpha_trace=rep(NA_real_, maxIter)
        }

        if (acceleration=="daarem") {
            daarem_state=list(w0=list(), residual=list(), since_restart=0L)
            daarem_order_trace=integer(maxIter)
            daarem_damping_trace=numeric(maxIter)
            daarem_restart_trace=logical(maxIter)
        }

        elbo_old=-Inf

        for (iter in 1:maxIter) {
            timeIn=proc.time()[3]

            w0_old=w0

            # Current residual cross-products
            XR_bar=vector("list", R)

            for (r in 1:R) {
                XR_bar[[r]]=XY[[r]]-XX[[r]]%*%B[,r]
            }

            # Active blocks
            if (iter<=2) {
                active_blocks=seq_len(n_blocks)
            } else {
                previous_gain=block_gain_trace[,iter-1]
                active_blocks=which(is.finite(previous_gain) &
                                    ((previous_gain>block_tol) |
                                     (previous_gain<0 & negative_gain_streak<negative_gain_patience)))
            }

            if (length(active_blocks)==0) {
                if (verbose) {
                    cat("No active blocks remain at iteration", iter, "\n")
                }
                iter=iter-1L
                break
            }

            # Visit active blocks
            for (b in active_blocks) {
                idx=blocks[[b]]

                B_block_start=B[idx,,drop=FALSE]
                Sj_block_start=Sj[idx]
                KL_block_start=KL[idx]

                XR_block_start=lapply(XR_bar, function(x) x[idx])
                XX_block=lapply(XX, function(x) x[idx,idx,drop=FALSE])
                XX_block_diag=XX_diag[idx,,drop=FALSE]

                block_theta=as.vector(B_block_start)

                if (version=="Rcpp") {
                    XX_block_rcpp=list2array.rssx(XX_block)
                    XR_block_start_rcpp=do.call(cbind, XR_block_start)
                    Sj_block_start_rcpp=Sj2array.rssx(Sj_block_start, R)
                }

                previous_gain=0
                block_gain_b=0

                # Repeated within-block sweeps
                for (sweep in 1:max_block_sweeps) {
                    if (version=="Rcpp") {
                        tmp_block=update_map_block_rcpp(block_theta=block_theta,
                                                        B_block_start=B_block_start,
                                                        Sj_block_start=Sj_block_start_rcpp,
                                                        KL_block_start=KL_block_start,
                                                        XR_block_start=XR_block_start_rcpp,
                                                        XX_block=XX_block_rcpp,
                                                        XX_block_diag=XX_block_diag,
                                                        v=v, S0=S0_rcpp, w0=w0)
                    } else {
                        tmp_block=update_map_block.rssx(block_theta=block_theta,
                                                        B_block_start=B_block_start,
                                                        Sj_block_start=Sj_block_start,
                                                        KL_block_start=KL_block_start,
                                                        XR_block_start=XR_block_start,
                                                        XX_block=XX_block,
                                                        XX_block_diag=XX_block_diag,
                                                        v=v, S0=S0, w0=w0)
                    }

                    block_theta=tmp_block$theta
                    block_gain_sweep=tmp_block$block_gain-previous_gain
                    previous_gain=tmp_block$block_gain
                    block_gain_b=tmp_block$block_gain

                    if (block_gain_sweep<=block_tol) {break}
                }

                block_evaluations_trace[b,iter]=sweep

                if (sweep>=max_block_sweeps & block_gain_sweep>block_tol) {
                    warning("Maximum number of sweeps reached for block ", b, " at iteration ", iter, ".")
                }

                B[idx,]=tmp_block$B
                Sj[idx]=tmp_block$Sj
                W[idx,]=tmp_block$W
                KL[idx]=tmp_block$KL

                delta_B=B[idx,,drop=FALSE]-B_block_start

                for (r in seq_len(R)) {
                    XR_bar[[r]]=XR_bar[[r]]-XX[[r]][,idx,drop=FALSE]%*%delta_B[,r]
                }

                ERSS=ERSS+tmp_block$delta_ERSS
                block_gain_trace[b,iter]=block_gain_b

                if (is.finite(block_gain_b) & block_gain_b<0) {
                    negative_gain_streak[b]=negative_gain_streak[b]+1L
                } else if (is.finite(block_gain_b)) {
                    negative_gain_streak[b]=0L
                }
            }

            active_blocks_trace[iter]=length(active_blocks)

            # Global updates
            w0_raw=colMeans(W)
            v=ERSS/n

            w0_residual=max(abs(w0_raw-w0_old))

            #### No acceleration
            if (acceleration=="none") {
                w0=w0_raw

            #### SQUAREM
            } else if (acceleration=="squarem") {
                tmp_acc=squarem_w0.rssx(w0=w0_old, w0_raw=w0_raw, state=squarem_state)

                w0=tmp_acc$w0
                squarem_state=tmp_acc$state
                alpha_trace[iter]=tmp_acc$alpha
                acceleration_accepted_trace[iter]=tmp_acc$accelerated

            #### DAAREM
            } else if (acceleration=="daarem") {
                tmp_acc=daarem_w0.rssx(w0=w0_old, w0_raw=w0_raw, state=daarem_state,
                                       order=daarem_order, damping=daarem_damping)

                w0=tmp_acc$w0
                daarem_state=tmp_acc$state
                daarem_order_trace[iter]=tmp_acc$order
                daarem_damping_trace[iter]=tmp_acc$damping
                daarem_restart_trace[iter]=tmp_acc$restarted
                acceleration_accepted_trace[iter]=tmp_acc$accelerated
            }

            if (acceleration!="none") {
                w0_residual_trace[iter]=w0_residual
            }

            ELBO[iter]=ELBO.rssx(v=v, ERSS=ERSS, KL=KL, n=n)

            if (iter==1) {
                delta_ELBO=NA_real_
            } else {
                delta_ELBO=ELBO[iter]-elbo_old
            }

            timeOut=proc.time()[3]

            if (verbose) {
                cat("iteration =", iter,
                    "| ELBO =", round(ELBO[iter], 3),
                    "| delta =", round(delta_ELBO, 3),
                    "| active blocks =", length(active_blocks),
                    "| time/iter =", round(timeOut-timeIn, 3), "s",
                    "\n")
            }

            if (iter>1 & abs(delta_ELBO)<tol) {break}

            elbo_old=ELBO[iter]
        }

        ELBO=ELBO[seq_len(iter)]
        block_gain_trace=block_gain_trace[,seq_len(iter),drop=FALSE]
        active_blocks_trace=active_blocks_trace[seq_len(iter)]
        block_evaluations_trace=block_evaluations_trace[,seq_len(iter),drop=FALSE]

        res=list(B=B,
                 v=v,
                 w0=w0,
                 w1=W,
                 Sj=Sj,
                 ELBO=ELBO[length(ELBO)],
                 ELBO_trace=ELBO,
                 iterations=iter,
                 update_scheme="blockwise",
                 acceleration=acceleration,
                 version=version,
                 block_size=block_size,
                 block_gain_trace=block_gain_trace,
                 active_blocks=active_blocks_trace,
                 block_evaluations_trace=block_evaluations_trace)

        if (acceleration!="none") {
            res$w0_residual_trace=w0_residual_trace[seq_len(iter)]
            res$acceleration_accepted_trace=acceleration_accepted_trace[seq_len(iter)]
        }

        if (acceleration=="squarem") {
            res$alpha_trace=alpha_trace[seq_len(iter)]
        }

        if (acceleration=="daarem") {
            res$daarem_order_trace=daarem_order_trace[seq_len(iter)]
            res$daarem_damping_trace=daarem_damping_trace[seq_len(iter)]
            res$daarem_restart_trace=daarem_restart_trace[seq_len(iter)]
        }
    }

    res
}
