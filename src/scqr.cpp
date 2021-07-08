# include <RcppArmadillo.h>
# include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
int sgn(const double x) {
  return (x > 0) - (x < 0);
}

// [[Rcpp::export]]
double mad(const arma::vec& x) {
  return 1.482602 * arma::median(arma::abs(x - arma::median(x)));
}

// [[Rcpp::export]]
arma::vec getH(const arma::vec& tauSeq) {
  return -arma::log(1 - tauSeq);
}

// [[Rcpp::export]]
void updateHuber(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const int n, const double tau, const double n1) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (std::abs(cur) <= tau) {
      der(i) = -cur;
    } else {
      der(i) = -tau * sgn(cur);
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
arma::vec huberReg(const arma::mat& Z, const arma::vec& Y, arma::vec& der, arma::vec& gradOld, arma::vec& gradNew, const int n, const int p, 
                   const double n1, const double tol = 0.0001, const double constTau = 1.345, const int iteMax = 5000) {
  double tau = constTau * mad(Y);
  updateHuber(Z, Y, der, gradOld, n, tau, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  tau = constTau * mad(res);
  updateHuber(Z, res, der, gradNew, n, tau, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    tau = constTau * mad(res);
    updateHuber(Z, res, der, gradNew, n, tau, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::mat standardize(arma::mat X, const arma::rowvec& mx, const arma::vec& sx, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - mx(i)) / sx(i);
  }
  return X;
}

// [[Rcpp::export]]
void updateGauss(const arma::mat& Z, const arma::uvec& censor, const arma::vec& res, const arma::vec& accu, arma::vec& der, arma::vec& grad, 
                 const double n1, const double h1) {
  der = censor % arma::normcdf(res * h1) - accu;
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
void updateUnif(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const int n, const double tau, const double h, 
                const double n1, const double h1) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (cur <= -h) {
      der(i) = 1 - tau;
    } else if (cur < h) {
      der(i) = 0.5 - tau - 0.5 * h1 * cur;
    } else {
      der(i) = -tau;
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
void updatePara(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const int n, const double tau, const double h, 
                const double n1, const double h1, const double h3) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (cur <= -h) {
      der(i) = 1 - tau;
    } else if (cur < h) {
      der(i) = 0.5 - tau - 0.75 * h1 * cur + 0.25 * h3 * cur * cur * cur;
    } else {
      der(i) = -tau;
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
void updateTrian(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const int n, const double tau, const double h, 
                 const double n1, const double h1, const double h2) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (cur <= -h) {
      der(i) = 1 - tau;
    } else if (cur < 0) {
      der(i) = 0.5 - tau - h1 * cur - 0.5 * h2 * cur * cur;
    } else if (cur < h) {
      der(i) = 0.5 - tau - h1 * cur + 0.5 * h2 * cur * cur;
    } else {
      der(i) = -tau;
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
arma::vec step0Gauss(const arma::mat& Z, const arma::vec& Y, const arma::uvec& censor, arma::vec& accu, const double tau, const double h, 
                     const int n, const int p, const double n1, const double h1, const double constTau = 1.345, const double tol = 0.0001, 
                     const int iteMax = 500) {
  arma::vec gradOld(p + 1), gradNew(p + 1);
  arma::vec der(n);
  arma::vec beta = huberReg(Z, Y, der, gradOld, gradNew, n, p, n1, tol, constTau, iteMax);
  arma::vec quant = {tau};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  arma::vec res = Z * beta - Y;
  accu = tau * arma::ones(n);
  updateGauss(Z, censor, res, accu, der, gradOld, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res += Z * betaDiff;
  updateGauss(Z, censor, res, accu, der, gradNew, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res += Z * betaDiff;
    updateGauss(Z, censor, res, accu, der, gradNew, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::vec step0GaussIni(const arma::mat& Z, const arma::vec& Y, const arma::uvec& censor, const arma::mat& betaMat, arma::vec& accu, const double tau, 
                        const double h, const int n, const int p, const double n1, const double h1, const double tol = 0.0001, const int iteMax = 500) {
  arma::vec gradOld(p + 1), gradNew(p + 1);
  arma::vec der(n);
  arma::vec beta = betaMat.col(0);
  arma::vec res = Z * beta - Y;
  accu = tau * arma::ones(n);
  updateGauss(Z, censor, res, accu, der, gradOld, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res += Z * betaDiff;
  updateGauss(Z, censor, res, accu, der, gradNew, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res += Z * betaDiff;
    updateGauss(Z, censor, res, accu, der, gradNew, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::vec stepkGauss(const arma::mat& Z, const arma::vec& Y, const arma::uvec& censor, arma::vec& beta, arma::vec& accu, const arma::vec& HSeq, 
                     const int k, const double h, const int n, const int p, const double n1, const double h1, const double tol = 0.0001, 
                     const int iteMax = 500) {
  arma::vec gradOld(p + 1), gradNew(p + 1);
  arma::vec der(n);
  arma::vec res = Z * beta - Y;
  accu += arma::normcdf(-res * h1) * (HSeq(k) - HSeq(k - 1));
  updateGauss(Z, censor, res, accu, der, gradOld, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res += Z * betaDiff;
  updateGauss(Z, censor, res, accu, der, gradNew, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res += Z * betaDiff;
    updateGauss(Z, censor, res, accu, der, gradNew, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::vec stepkGaussIni(const arma::mat& Z, const arma::vec& Y, const arma::uvec& censor, const arma::mat& betaMat, arma::vec& beta, arma::vec& accu, 
                        const arma::vec& HSeq, const int k, const double h, const int n, const int p, const double n1, const double h1, const double tol = 0.0001, 
                        const int iteMax = 500) {
  arma::vec gradOld(p + 1), gradNew(p + 1);
  arma::vec der(n);
  arma::vec res = Z * beta - Y;
  accu += arma::normcdf(-res * h1) * (HSeq(k) - HSeq(k - 1));
  beta = betaMat.col(k);
  res = Z * beta - Y;
  updateGauss(Z, censor, res, accu, der, gradOld, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res += Z * betaDiff;
  updateGauss(Z, censor, res, accu, der, gradNew, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res += Z * betaDiff;
    updateGauss(Z, censor, res, accu, der, gradNew, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
Rcpp::List scqrGauss(const arma::mat& X, arma::vec Y, const arma::uvec& censor, const arma::vec& tauSeq, double h = 0.0, const double constTau = 1.345, 
                     const double tol = 0.0001, const int iteMax = 500) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int m = tauSeq.size();
  if (h <= 0.05) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double n1 = 1.0 / n;
  const double h1 = 1.0 / h;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx = arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec accu(n);
  arma::mat betaProc(p + 1, m);
  arma::vec beta = step0Gauss(Z, Y, censor, accu, tauSeq(0), h, n, p, n1, h1, constTau, tol, iteMax);
  betaProc.col(0) = beta;
  arma::vec HSeq = getH(tauSeq);
  for (int k = 1; k < m; k++) {
    beta = stepkGauss(Z, Y, censor, beta, accu, HSeq, k, h, n, p, n1, h1, tol, iteMax);
    betaProc.col(k) =  beta;
  }
  betaProc.rows(1, p).each_col() /= sx;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return Rcpp::List::create(Rcpp::Named("coeff") = betaProc, Rcpp::Named("bandwidth") = h);
}

// [[Rcpp::export]]
arma::vec scqrGaussIni(const arma::mat& X, arma::vec Y, const arma::uvec& censor, const arma::mat& betaProc, const arma::vec& tauSeq, 
                        const arma::vec& HSeq, const double h, const int p, const int m, const double h1, const double tol = 0.0001, const int iteMax = 500) {
  const int n = X.n_rows;
  const double n1 = 1.0 / n;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx = arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec accu(n);
  arma::vec beta = step0GaussIni(Z, Y, censor, betaProc, accu, tauSeq(0), h, n, p, n1, h1, tol, iteMax);
  for (int k = 1; k < m; k++) {
    beta = stepkGaussIni(Z, Y, censor, betaProc, beta, accu, HSeq, k, h, n, p, n1, h1, tol, iteMax);
  }
  beta.rows(1, p) /= sx;
  beta(0) += my - arma::as_scalar(mx * beta.rows(1, p));
  return beta;
}

// [[Rcpp::export]]
Rcpp::List scqrGaussInf(const arma::mat& X, arma::vec Y, const arma::uvec& censor, const arma::vec& tauSeq, const int B = 1000, double h = 0.0, 
                       const double constTau = 1.345, const double tol = 0.0001, const int iteMax = 500) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int m = tauSeq.size();
  if (h <= 0.05) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double n1 = 1.0 / n;
  const double h1 = 1.0 / h;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx = arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec accu(n);
  arma::mat betaProc(p + 1, m);
  arma::vec beta = step0Gauss(Z, Y, censor, accu, tauSeq(0), h, n, p, n1, h1, constTau, tol, iteMax);
  betaProc.col(0) = beta;
  arma::vec HSeq = getH(tauSeq);
  for (int k = 1; k < m; k++) {
    beta = stepkGauss(Z, Y, censor, beta, accu, HSeq, k, h, n, p, n1, h1, tol, iteMax);
    betaProc.col(k) =  beta;
  }
  betaProc.rows(1, p).each_col() /= sx;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  Y += my;
  arma::mat rst(p + 1, B);
  for (int b = 0; b < B; b++) {
    arma::uvec idx = arma::find(arma::randi(n, arma::distr_param(0, 1)) == 1);
    arma::mat mbX = X.rows(idx);
    arma::vec mbY = Y.rows(idx);
    arma::uvec mbCensor = censor.rows(idx);
    rst.col(b) = scqrGaussIni(mbX, mbY, mbCensor, betaProc, tauSeq, HSeq, h, p, m, h1, tol, iteMax);
  }
  return Rcpp::List::create(Rcpp::Named("coeff") = betaProc, Rcpp::Named("boot") = rst, Rcpp::Named("bandwidth") = h);
}

