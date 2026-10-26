squarem_w0.rssx=function(w0, w0_raw, state) {
    K=length(w0)

    w0=pmax(w0, 0)
    if (any(!is.finite(w0)) | sum(w0)<=0) {
        w0=rep(1/K, K)
    } else {
        w0=w0/sum(w0)
    }

    w0_raw=pmax(w0_raw, 0)
    if (any(!is.finite(w0_raw)) | sum(w0_raw)<=0) {
        w0_raw=rep(1/K, K)
    } else {
        w0_raw=w0_raw/sum(w0_raw)
    }

    alpha=NA_real_
    accelerated=FALSE

    if (state$phase==0L) {
        state$w00=w0
        state$w01=w0_raw
        state$phase=1L

        w0_next=w0_raw

    } else {
        r=state$w01-state$w00
        u=(w0_raw-w0)-r

        nr=sqrt(sum(r^2))
        nu=sqrt(sum(u^2))

        if (is.finite(nr) & is.finite(nu) & nu>0) {
            alpha=min(-1, -nr/nu)
            w0_square=state$w00-2*alpha*r+alpha^2*u
            w0_square=pmax(w0_square, 0)

            if (all(is.finite(w0_square)) & sum(w0_square)>0) {
                w0_next=w0_square/sum(w0_square)
                accelerated=TRUE
            } else {
                w0_next=w0_raw
            }
        } else {
            w0_next=w0_raw
        }

        state$w00=NULL
        state$w01=NULL
        state$phase=0L
    }

    list(w0=w0_next,
         state=state,
         alpha=alpha,
         accelerated=accelerated)
}

daarem_w0.rssx=function(w0, w0_raw, state, order=10, damping=0.5) {
    K=length(w0)

    w0=pmax(w0, 0)
    if (any(!is.finite(w0)) | sum(w0)<=0) {
        w0=rep(1/K, K)
    } else {
        w0=w0/sum(w0)
    }

    w0_raw=pmax(w0_raw, 0)
    if (any(!is.finite(w0_raw)) | sum(w0_raw)<=0) {
        w0_raw=rep(1/K, K)
    } else {
        w0_raw=w0_raw/sum(w0_raw)
    }

    residual=w0_raw-w0

    accelerated=FALSE
    restarted=FALSE
    order_used=0L
    damping_used=0

    if (state$since_restart>=2*order) {
        state$w0=list()
        state$residual=list()
        state$since_restart=0L
        restarted=TRUE
    }

    w0_history=c(state$w0, list(w0))
    residual_history=c(state$residual, list(residual))

    nh=length(residual_history)
    order_used=min(order, nh-1L)

    w0_next=w0_raw

    if (order_used>0) {
        idx=(nh-order_used):(nh-1L)

        X=matrix(0, nrow=K, ncol=order_used)
        Rmat=matrix(0, nrow=K, ncol=order_used)

        for (j in seq_len(order_used)) {
            i=idx[j]
            X[,j]=w0_history[[i+1L]]-w0_history[[i]]
            Rmat[,j]=residual_history[[i+1L]]-residual_history[[i]]
        }

        RtR=crossprod(Rmat)
        scale=mean(diag(RtR))

        if (!is.finite(scale) | scale<=0) {
            scale=1
        }

        ridge=1e-8*scale
        gamma=try(solve(RtR+diag(ridge, order_used), crossprod(Rmat, residual)), silent=TRUE)

        if (!inherits(gamma, "try-error") & all(is.finite(gamma))) {
            damping_used=damping+(1-damping)*min(1, state$since_restart/max(1, order))
            w0_daarem=w0_raw-damping_used*drop((X+Rmat)%*%gamma)
            w0_daarem=pmax(w0_daarem, 0)

            if (all(is.finite(w0_daarem)) & sum(w0_daarem)>0) {
                w0_next=w0_daarem/sum(w0_daarem)
                accelerated=TRUE
            } else {
                state$w0=list()
                state$residual=list()
                state$since_restart=0L
                restarted=TRUE
            }
        } else {
            state$w0=list()
            state$residual=list()
            state$since_restart=0L
            restarted=TRUE
        }
    }

    state$w0=c(state$w0, list(w0))
    state$residual=c(state$residual, list(residual))

    if (length(state$w0)>order+1L) {
        keep=(length(state$w0)-order):length(state$w0)
        state$w0=state$w0[keep]
        state$residual=state$residual[keep]
    }

    state$since_restart=state$since_restart+1L

    list(w0=w0_next,
         state=state,
         order=order_used,
         damping=damping_used,
         accelerated=accelerated,
         restarted=restarted)
}
