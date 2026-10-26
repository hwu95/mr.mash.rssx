#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

struct BMSRres {
    arma::vec b;
    arma::mat S;
    arma::vec w;
    double logBF_mix;
};

std::vector<arma::mat> list2mat(const List& x, arma::uword nr, arma::uword nc) {
    std::vector<arma::mat> out;
    out.reserve(x.size());

    for (int i=0; i<x.size(); i++) {
        NumericMatrix xi=x[i];
        if ((arma::uword)xi.nrow()!=nr || (arma::uword)xi.ncol()!=nc) {
            stop("XX has wrong dimension") ;
        }
        out.emplace_back(xi.begin(), nr, nc, false, true);
    }

    return out;
}

List cube2list(const arma::cube& x) {
    List out(x.n_slices);

    for (arma::uword j=0; j<x.n_slices; j++) {
        out[j]=wrap(x.slice(j));
    }

    return out;
}

BMSRres BMSR_cpp(const arma::vec& XRj_bar,
                 const arma::vec& xx,
                 const arma::vec& v,
                 const arma::cube& S0,
                 const arma::vec& w0) {

    arma::uword R=XRj_bar.n_elem;
    arma::uword K=S0.n_slices;

    arma::vec bj_hat=XRj_bar/xx;
    arma::vec shat=v/xx;
    arma::mat Sj_hat=arma::diagmat(shat);

    arma::cube S1(R, R, K, arma::fill::zeros);
    arma::mat b1(R, K, arma::fill::zeros);
    arma::vec logBF(K, arma::fill::zeros);

    double log_det_Sj_hat=arma::accu(arma::log(shat));
    double quad_Sj_hat=arma::accu(arma::square(bj_hat)%(xx/v));

    for (arma::uword k=0; k<K; k++) {
        arma::mat S0k=S0.slice(k);

        if (arma::abs(S0k).max()==0) {
            continue;
        }

        arma::mat A=S0k+Sj_hat;
        arma::mat U;

        if (!arma::chol(U, A)) {
            stop("chol failed") ;
        }

        arma::mat tmp=arma::solve(arma::trimatl(U.t()), Sj_hat, arma::solve_opts::fast);
        arma::mat Ainv_Shat=arma::solve(arma::trimatu(U), tmp, arma::solve_opts::fast);
        S1.slice(k)=S0k*Ainv_Shat;

        arma::vec z=arma::solve(arma::trimatl(U.t()), bj_hat, arma::solve_opts::fast);
        arma::vec Ainv_b=arma::solve(arma::trimatu(U), z, arma::solve_opts::fast);
        b1.col(k)=S0k*Ainv_b;

        double log_det_A=2*arma::accu(arma::log(U.diag()));

        logBF(k)=0.5*log_det_Sj_hat-
                 0.5*log_det_A-
                 0.5*arma::dot(z,z)+
                 0.5*quad_Sj_hat;
    }

    arma::vec log_post=arma::log(w0)+logBF;
    double m=log_post.max();

    arma::vec w1=arma::exp(log_post-m);
    double sw=arma::accu(w1);
    w1/=sw;

    double logBF_mix=m+std::log(sw);

    arma::vec bj_bar=b1*w1;
    arma::mat Sj(R, R, arma::fill::zeros);

    for (arma::uword k=0; k<K; k++) {
        arma::vec bk=b1.col(k);
        Sj+=w1(k)*(S1.slice(k)+bk*bk.t());
    }

    Sj-=bj_bar*bj_bar.t();

    return {bj_bar, Sj, w1, logBF_mix};
}

// [[Rcpp::export]]
List update_map_coordinate_rcpp(const arma::mat& B,
                                const arma::vec& w0,
                                const arma::vec& v,
                                const List& XX,
                                const arma::mat& XX_diag,
                                const arma::mat& XY,
                                const arma::vec& YY,
                                const arma::vec& n,
                                const arma::cube& S0) {

    arma::uword p=B.n_rows;
    arma::uword R=B.n_cols;
    arma::uword K=S0.n_slices;

    std::vector<arma::mat> XX_mat=list2mat(XX, p, p);

    arma::mat B_new=B;
    arma::mat XR_bar(p, R, arma::fill::zeros);

    for (arma::uword r=0; r<R; r++) {
        XR_bar.col(r)=XY.col(r)-XX_mat[r]*B_new.col(r);
    }

    arma::cube Sj(R, R, p, arma::fill::zeros);
    arma::mat W(p, K, arma::fill::zeros);
    arma::vec KL(p, arma::fill::zeros);

    for (arma::uword j=0; j<p; j++) {
        arma::vec bj=B_new.row(j).t();
        arma::vec xx=XX_diag.row(j).t();
        arma::vec XRj_bar(R);

        for (arma::uword r=0; r<R; r++) {
            XRj_bar(r)=XR_bar(j,r)+xx(r)*bj(r);
        }

        BMSRres tmp=BMSR_cpp(XRj_bar, xx, v, S0, w0);

        B_new.row(j)=tmp.b.t();
        Sj.slice(j)=tmp.S;
        W.row(j)=tmp.w.t();

        arma::vec delta_ERSS=-2*(tmp.b%XRj_bar)+xx%(tmp.S.diag()+arma::square(tmp.b));
        KL(j)=-0.5*arma::accu(delta_ERSS/v)-tmp.logBF_mix;

        arma::vec db=bj-tmp.b;

        for (arma::uword r=0; r<R; r++) {
            XR_bar.col(r)+=XX_mat[r].col(j)*db(r);
        }
    }

    arma::vec w0_new=arma::mean(W,0).t();
    arma::vec ERSS(R, arma::fill::zeros);

    for (arma::uword r=0; r<R; r++) {
        arma::vec br=B_new.col(r);

        ERSS(r)=YY(r)-2*arma::dot(br,XY.col(r))+
                arma::as_scalar(br.t()*XX_mat[r]*br);

        for (arma::uword j=0; j<p; j++) {
            ERSS(r)+=XX_diag(j,r)*Sj(r,r,j);
        }
    }

    arma::vec v_new=ERSS/n;

    return List::create(_["B"]=B_new,
                        _["v"]=v_new,
                        _["w0"]=w0_new,
                        _["W"]=W,
                        _["Sj"]=cube2list(Sj),
                        _["KL"]=KL,
                        _["ERSS"]=ERSS);
}

// [[Rcpp::export]]
List update_map_block_rcpp(const arma::vec& block_theta,
                           const arma::mat& B_block_start,
                           const arma::cube& Sj_block_start,
                           const arma::vec& KL_block_start,
                           const arma::mat& XR_block_start,
                           const arma::cube& XX_block,
                           const arma::mat& XX_block_diag,
                           const arma::vec& v,
                           const arma::cube& S0,
                           const arma::vec& w0) {

    arma::uword m=B_block_start.n_rows;
    arma::uword R=B_block_start.n_cols;
    arma::uword K=S0.n_slices;

    arma::mat B=arma::reshape(block_theta, m, R);
    arma::mat XR_bar(m, R, arma::fill::zeros);

    for (arma::uword r=0; r<R; r++) {
        XR_bar.col(r)=XR_block_start.col(r)+
                      XX_block.slice(r)*(B_block_start.col(r)-B.col(r));
    }

    arma::cube Sj(R, R, m, arma::fill::zeros);
    arma::mat W(m, K, arma::fill::zeros);
    arma::vec KL(m, arma::fill::zeros);

    for (arma::uword j=0; j<m; j++) {
        arma::vec bj=B.row(j).t();
        arma::vec xx=XX_block_diag.row(j).t();
        arma::vec XRj_bar(R);

        for (arma::uword r=0; r<R; r++) {
            XRj_bar(r)=XR_bar(j,r)+xx(r)*bj(r);
        }

        BMSRres tmp=BMSR_cpp(XRj_bar, xx, v, S0, w0);

        B.row(j)=tmp.b.t();
        Sj.slice(j)=tmp.S;
        W.row(j)=tmp.w.t();

        arma::vec ERSS_j=-2*(tmp.b%XRj_bar)+xx%(tmp.S.diag()+arma::square(tmp.b));
        KL(j)=-0.5*arma::accu(ERSS_j/v)-tmp.logBF_mix;

        arma::vec db=bj-tmp.b;

        for (arma::uword r=0; r<R; r++) {
            XR_bar.col(r)+=XX_block.slice(r).col(j)*db(r);
        }
    }

    arma::mat delta_B=B-B_block_start;
    arma::vec delta_ERSS(R, arma::fill::zeros);

    for (arma::uword r=0; r<R; r++) {
        arma::vec d=delta_B.col(r);
        arma::vec delta_sj(m);

        for (arma::uword j=0; j<m; j++) {
            delta_sj(j)=Sj(r,r,j)-Sj_block_start(r,r,j);
        }

        delta_ERSS(r)=-2*arma::dot(d,XR_block_start.col(r))+
                      arma::dot(d,XX_block.slice(r)*d)+
                      arma::dot(XX_block_diag.col(r),delta_sj);
    }

    double delta_KL=arma::accu(KL-KL_block_start);
    double block_gain=-0.5*arma::accu(delta_ERSS/v)-delta_KL;

    return List::create(_["theta"]=arma::vectorise(B),
                        _["B"]=B,
                        _["Sj"]=cube2list(Sj),
                        _["W"]=W,
                        _["KL"]=KL,
                        _["XR"]=XR_bar,
                        _["delta_ERSS"]=delta_ERSS,
                        _["delta_KL"]=delta_KL,
                        _["block_gain"]=block_gain);
}
