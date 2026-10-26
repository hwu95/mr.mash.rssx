update_map_coordinate.rssx=function(B, w0, v, XX, XY, YY, n, S0) {
    p=nrow(B)
    R=ncol(B)
    K=length(S0)

    Sj=vector("list", p)
    W=matrix(0, nrow=p, ncol=K)
    KL=numeric(p)
    XR_bar=vector("list", R)

    for (r in 1:R) {
        XR_bar[[r]]=XY[[r]]-XX[[r]]%*%B[,r]
    }

    tmp=update_coordinate.rssx(XX=XX, XR_bar=XR_bar, B=B, Sj=Sj, W=W, KL=KL, v=v, S0=S0, w0=w0)

    B=tmp$B
    Sj=tmp$Sj
    W=tmp$W
    KL=tmp$KL

    w0=colMeans(W)

    tmp_v=updateV.rssx.coordinate(XX=XX, XY=XY, YY=YY, B=B, Sj=Sj, n=n)
    v=tmp_v$v
    ERSS=tmp_v$ERSS

    list(B=B,
         v=v,
         w0=w0,
         W=W,
         Sj=Sj,
         KL=KL,
         ERSS=ERSS)
}

update_coordinate.rssx=function(XX, XR_bar, B, Sj, W, KL, v, S0, w0) {
    p=nrow(B)

    for (j in 1:p) {
        tmp=BMSR.rssx.coordinate(XX=XX, j=j, XR_bar=XR_bar, bj=B[j,], v=v, S0=S0, w0=w0)

        B[j,]=tmp$b
        Sj[[j]]=tmp$S
        W[j,]=tmp$w
        KL[j]=tmp$KL
        XR_bar=tmp$XR
    }

    list(B=B,
         Sj=Sj,
         W=W,
         KL=KL)
}

BMSR.rssx.coordinate=function(XX, j, XR_bar, bj, v, S0, w0) {
    K=length(S0)
    R=length(bj)
    xx=numeric(R)

    for (r in 1:R) {
        xx[r]=XX[[r]][j,j]
    }

    XRj_bar=numeric(R)

    for (r in 1:R) {
        XRj_bar[r]=XR_bar[[r]][j]+xx[r]*bj[r]
    }

    bj_hat=XRj_bar/xx
    Sj_hat=diag(v/xx, nrow=R, ncol=R)

    S1_j=vector("list", K)
    b1_j=matrix(0, nrow=R, ncol=K)
    logBF=numeric(K)

    log_det_Sj_hat=sum(log(v/xx))
    quad_Sj_hat=sum(bj_hat^2*(xx/v))

    for (k in 1:K) {
        if (max(abs(S0[[k]]))==0) {
            S1_j[[k]]=matrix(0, R, R)
            b1_j[,k]=0
            logBF[k]=0
        } else {
            S0_plus_Shat=S0[[k]]+Sj_hat
            S0_plus_Shat_chol=chol(S0_plus_Shat)

            S1_j[[k]]=S0[[k]]%*%backsolve(S0_plus_Shat_chol, forwardsolve(t(S0_plus_Shat_chol), Sj_hat))
            z=forwardsolve(t(S0_plus_Shat_chol), bj_hat)
            b1_j[,k]=S0[[k]]%*%backsolve(S0_plus_Shat_chol, z)

            log_det_S0_plus_Shat=2*sum(log(diag(S0_plus_Shat_chol)))

            logBF[k]=0.5*log_det_Sj_hat-0.5*log_det_S0_plus_Shat-
                     0.5*sum(z^2)+0.5*quad_Sj_hat
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

    Sj_j=Sj_j-tcrossprod(bj_bar)

    delta_ERSS=-2*bj_bar*XRj_bar+xx*(diag(Sj_j)+bj_bar^2)
    KL_j=-0.5*sum(delta_ERSS/v)-logBF_mix

    bj_old=bj
    XR_bar_new=vector("list", R)

    for (r in 1:R) {
        XR_bar_new[[r]]=XR_bar[[r]]+drop(XX[[r]][,j,drop=FALSE]*(bj_old[r]-bj_bar[r]))
    }

    list(b=bj_bar,
         S=Sj_j,
         w=w1,
         KL=KL_j,
         XR=XR_bar_new,
         bj_hat=bj_hat,
         Sj_hat=Sj_hat,
         logBF=logBF,
         logBF_mix=logBF_mix)
}

updateV.rssx.coordinate=function(XX, XY, YY, B, Sj, n) {
    p=nrow(B)
    R=ncol(B)

    ERSS=numeric(R)
    v=numeric(R)

    for (r in 1:R) {
        br=B[,r]
        ERSS[r]=YY[r]-2*sum(br*XY[[r]])+drop(t(br)%*%XX[[r]]%*%br)

        for (j in 1:p) {
            ERSS[r]=ERSS[r]+XX[[r]][j,j]*Sj[[j]][r,r]
        }

        v[r]=ERSS[r]/n[r]
    }

    list(v=v,
         ERSS=ERSS)
}
