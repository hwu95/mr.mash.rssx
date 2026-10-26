update_map_block.rssx=function(block_theta, B_block_start, Sj_block_start, KL_block_start, XR_block_start,
                               XX_block, XX_block_diag, v, S0, w0) {
    m=nrow(B_block_start)
    R=ncol(B_block_start)
    K=length(S0)

    B_block=matrix(block_theta, nrow=m, ncol=R)

    XR_block=vector("list", R)

    for (r in seq_len(R)) {
        XR_block[[r]]=XR_block_start[[r]]+XX_block[[r]]%*%(B_block_start[,r]-B_block[,r])
    }

    Sj_block=vector("list", m)
    W_block=matrix(0, nrow=m, ncol=K)
    KL_block=numeric(m)

    for (j_local in seq_len(m)) {
        tmp=BMSR.rssx.blockwise(XX=XX_block, XX_diag=XX_block_diag, j=j_local, XR_bar=XR_block,
                                bj=B_block[j_local,], Sj_old=NULL, KL_old=0, v=v, S0=S0, w0=w0)

        B_block[j_local,]=tmp$b
        Sj_block[[j_local]]=tmp$S
        W_block[j_local,]=tmp$w
        KL_block[j_local]=tmp$KL
        XR_block=tmp$XR
    }

    delta_ERSS=block_delta_ERSS.rssx(XX_block=XX_block, XX_diag_block=XX_block_diag,
                                     XR_block_start=XR_block_start, B_start=B_block_start,
                                     B_end=B_block, Sj_start=Sj_block_start, Sj_end=Sj_block)

    delta_KL=sum(KL_block-KL_block_start)
    block_gain=-0.5*sum(delta_ERSS/v)-delta_KL

    list(theta=as.vector(B_block),
         B=B_block,
         Sj=Sj_block,
         W=W_block,
         KL=KL_block,
         XR=XR_block,
         delta_ERSS=delta_ERSS,
         delta_KL=delta_KL,
         block_gain=block_gain)
}

BMSR.rssx.blockwise=function(XX, XX_diag, j, XR_bar, bj, Sj_old=NULL, KL_old=0, v, S0, w0) {
    K=length(S0)
    R=length(bj)
    xx=XX_diag[j,]

    XRj_bar=numeric(R)

    for (r in 1:R) {
        XRj_bar[r]=XR_bar[[r]][j]+xx[r]*bj[r]
    }

    bj_hat=XRj_bar/xx

    precision=xx/v
    Sj_hat=diag(1/precision)

    S1_j=vector("list", K)
    b1_j=matrix(0, nrow=R, ncol=K)
    logBF=numeric(K)

    for (k in seq_len(K)) {
        if (max(abs(S0[[k]]))==0) {
            S1_j[[k]]=matrix(0, R, R)
            b1_j[,k]=0
            logBF[k]=0
        } else {
            S0_plus_Shat=S0[[k]]+Sj_hat
            S0_plus_Shat_chol=chol(S0_plus_Shat)

            S1_j[[k]]=S0[[k]]%*%backsolve(S0_plus_Shat_chol, forwardsolve(t(S0_plus_Shat_chol), Sj_hat))
            b1_j[,k]=S1_j[[k]]%*%(precision*bj_hat)

            log_det_S0_plus_Shat=2*sum(log(diag(S0_plus_Shat_chol)))

            logBF[k]=-0.5*sum(log(precision))-0.5*log_det_S0_plus_Shat+
                     0.5*sum(bj_hat^2*precision)-0.5*sum((forwardsolve(t(S0_plus_Shat_chol), bj_hat))^2)
        }
    }

    w0_log=log(w0)
    log_post=w0_log+logBF
    m=max(log_post)
    w1=exp(log_post-m)
    w1=w1/sum(w1)
    logBF_mix=m+log(sum(exp(log_post-m)))

    bj_bar=drop(b1_j%*%w1)

    Sj_j=matrix(0, nrow=R, ncol=R)

    for (k in 1:K) {
        bk=matrix(b1_j[,k], nrow=R)
        Sj_j=Sj_j+w1[k]*(S1_j[[k]]+tcrossprod(bk))
    }

    Sj_j=Sj_j-bj_bar%*%t(bj_bar)

    ERSS_j_new=-2*bj_bar*XRj_bar+xx*(diag(Sj_j)+bj_bar^2)

    if (is.null(Sj_old)) {
        Sj_old=matrix(0, nrow=R, ncol=R)
    }

    ERSS_j_old=-2*bj*XRj_bar+xx*(diag(Sj_old)+bj^2)
    delta_ERSS=ERSS_j_new-ERSS_j_old

    KL_j=-0.5*sum(ERSS_j_new/v)-logBF_mix
    delta_KL=KL_j-KL_old
    delta_ELBO=-0.5*sum(delta_ERSS/v)-delta_KL

    bj_old=bj
    XR_bar_new=vector("list", R)

    for (r in 1:R) {
        XR_bar_new[[r]]=XR_bar[[r]]+XX[[r]][,j]*(bj_old[r]-bj_bar[r])
    }

    list(b=bj_bar,
         S=Sj_j,
         w=w1,
         KL=KL_j,
         XR=XR_bar_new,
         bj_hat=bj_hat,
         Sj_hat=Sj_hat,
         logBF=logBF,
         logBF_mix=logBF_mix,
         ERSS_j=ERSS_j_new,
         delta_ERSS=delta_ERSS,
         delta_KL=delta_KL,
         delta_ELBO=delta_ELBO)
}

block_delta_ERSS.rssx=function(XX_block, XX_diag_block, XR_block_start, B_start, B_end, Sj_start, Sj_end) {
    R=ncol(B_end)
    delta_B=B_end-B_start
    delta_ERSS=numeric(R)

    for (r in seq_len(R)) {
        d=delta_B[,r]
        sj_start=vapply(Sj_start, function(s) {if (is.null(s)) 0 else s[r,r]}, numeric(1))
        sj_end=vapply(Sj_end, function(s) {if (is.null(s)) 0 else s[r,r]}, numeric(1))

        delta_ERSS[r]=-2*sum(d*XR_block_start[[r]])+
                      drop(crossprod(d, XX_block[[r]]%*%d))+
                      sum(XX_diag_block[,r]*(sj_end-sj_start))
    }

    delta_ERSS
}
