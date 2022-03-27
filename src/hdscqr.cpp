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
arma::mat standardize(arma::mat X, const arma::rowvec& mx, const arma::vec& sx1, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - mx(i)) * sx1(i);
  }
  return X;
}

// [[Rcpp::export]]
arma::vec softThresh(const arma::vec& x, const arma::vec& lambda, const int p) {
  return arma::sign(x) % arma::max(arma::abs(x) - lambda, arma::zeros(p + 1));
}

// [[Rcpp::export]]
arma::vec cmptLambdaLasso(const double lambda, const int p) {
  arma::vec rst = lambda * arma::ones(p + 1);
  rst(0) = 0;
  return rst;
}

// [[Rcpp::export]]
arma::vec cmptLambdaSCAD(const arma::vec& beta, const double lambda, const int p, const double para = 3.7) {
  arma::vec rst = arma::zeros(p + 1);
  for (int i = 1; i <= p; i++) {
    double abBeta = std::abs(beta(i));
    if (abBeta <= lambda) {
      rst(i) = lambda;
    } else if (abBeta <= para * lambda) {
      rst(i) = (para * lambda - abBeta) / (para - 1);
    } 
  }
  return rst;
}

// [[Rcpp::export]]
arma::vec cmptLambdaMCP(const arma::vec& beta, const double lambda, const int p, const double para = 3.0) {
  arma::vec rst = arma::zeros(p + 1);
  for (int i = 1; i <= p; i++) {
    double abBeta = std::abs(beta(i));
    if (abBeta <= para * lambda) {
      rst(i) = lambda - abBeta / para;
    }
  }
  return rst;
}

// Huberize least square regression for the lowest quantile tau_L
// [[Rcpp::export]]
double lossL2(const arma::mat& Z, const arma::vec& Y, const arma::vec& beta, const double n1, const double tau) {
  arma::vec res = Y - Z * beta;
  double rst = 0.0;
  for (int i = 0; i < Y.size(); i++) {
    rst += (res(i) > 0) ? (tau * res(i) * res(i)) : ((1 - tau) * res(i) * res(i));
  }
  return 0.5 * n1 * rst;
}

// [[Rcpp::export]]
double updateL2(const arma::mat& Z, const arma::vec& Y, const arma::vec& beta, arma::vec& grad, const double n1, const double tau) {
  arma::vec res = Y - Z * beta;
  double rst = 0.0;
  grad = arma::zeros(grad.size());
  for (int i = 0; i < Y.size(); i++) {
    double temp = res(i) > 0 ? tau : (1 - tau);
    grad -= temp * res(i) * Z.row(i).t();
    rst += temp * res(i) * res(i);
  }
  grad *= n1;
  return 0.5 * n1 * rst;
}

// [[Rcpp::export]]
arma::vec indicator(const arma::vec& x, const int n) {
  arma::vec rst(n);
  for (int i = 0; i < n; i++) {
    rst(i) = x(i) >= 0 ? 1.0 : 0.0;
  }
  return rst;
}

// [[Rcpp::export]]
double calRes(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const arma::mat& betaHat, const arma::vec& tauSeq, const arma::vec& HSeq,
              const int m, const int n) {
  arma::vec indi = indicator(Y - Z * betaHat.col(0), n);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  arma::vec res = censor % (1 - indi) - tauSeq(0);
  arma::vec dev = arma::sqrt(-2 * (res + censor % arma::log(censor - res)));
  double rst = arma::mean(dev);
  for (int i = 1; i < m; i++) {
    accu += indi * (HSeq(i) - HSeq(i - 1));
    indi = indicator(Y - Z * betaHat.col(i), n);
    res = censor % (1 - indi) - accu;
    dev = arma::sqrt(-2 * (res + censor % arma::log(censor - res)));
    rst += arma::mean(dev);
  }
  return rst;
}

// [[Rcpp::export]]
double lossGauss(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const arma::vec& beta, const double h, 
                 const double h1, const double h2) {
  arma::vec res = Y - Z * beta;
  arma::vec temp = 0.3989423 * h  * arma::exp(-0.5 * h2 * arma::square(res)) + 0.5 * res - res % arma::normcdf(-h1 * res);
  return arma::mean(censor % temp);
}

// [[Rcpp::export]]
double updateGauss(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const arma::vec& accu, const arma::vec& beta, arma::vec& grad, 
                   arma::vec& gradReal, const double n1, const double h, const double h1, const double h2) {
  arma::vec res = Y - Z * beta;
  arma::vec der = censor % arma::normcdf(-h1 * res) - accu;
  gradReal = n1 * Z.t() * der;
  der = censor % arma::normcdf(-h1 * res) - 0.5 * censor;
  grad = n1 * Z.t() * der;
  arma::vec temp = 0.3989423 * h  * arma::exp(-0.5 * h2 * arma::square(res)) + 0.5 * res - res % arma::normcdf(-h1 * res);
  return arma::mean(censor % temp);
}

// [[Rcpp::export]]
void updateGaussLow(const arma::mat& Z, const arma::uvec& censor, const arma::vec& res, const arma::vec& accu, arma::vec& der, arma::vec& grad, 
                    const double n1, const double h1) {
  der = censor % arma::normcdf(res * h1) - accu;
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
  updateGaussLow(Z, censor, res, accu, der, gradOld, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res += Z * betaDiff;
  updateGaussLow(Z, censor, res, accu, der, gradNew, n1, h1);
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
    updateGaussLow(Z, censor, res, accu, der, gradNew, n1, h1);
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
  updateGaussLow(Z, censor, res, accu, der, gradOld, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res += Z * betaDiff;
  updateGaussLow(Z, censor, res, accu, der, gradNew, n1, h1);
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
    updateGaussLow(Z, censor, res, accu, der, gradNew, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::mat scqrGauss(const arma::mat& X, arma::vec Y, const arma::uvec& censor, const arma::vec& tauSeq, double h = 0.0, const double constTau = 1.345, 
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
  return betaProc;
}

// LAMM, update beta, return phi
// [[Rcpp::export]]
double lammL2(const arma::mat& Z, const arma::vec& Y, const arma::vec& Lambda, arma::vec& beta, const double tau, const double phi, const double gamma, 
              const int p, const double n1) {
  double phiNew = phi;
  arma::vec betaNew(p + 1);
  arma::vec grad(p + 1);
  double loss = updateL2(Z, Y, beta, grad, n1, tau);
  while (true) {
    arma::vec first = beta - grad / phiNew;
    arma::vec second = Lambda / phiNew;
    betaNew = softThresh(first, second, p);
    double fVal = lossL2(Z, Y, betaNew, n1, tau);
    arma::vec diff = betaNew - beta;
    double psiVal = loss + arma::as_scalar(grad.t() * diff) + 0.5 * phiNew * arma::as_scalar(diff.t() * diff);
    if (fVal <= psiVal) {
      break;
    }
    phiNew *= gamma;
  }
  beta = betaNew;
  return phiNew;
}

// [[Rcpp::export]]
double lammSq(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const arma::vec& Lambda, const arma::vec& accu, arma::vec& beta, 
              const double phi, const double gamma, const int p, const double h, const double n1, const double h1, const double h2) {
  double phiNew = phi;
  arma::vec betaNew(p + 1);
  arma::vec grad(p + 1);
  arma::vec gradReal(p + 1);
  double loss = updateGauss(Z, censor, Y, accu, beta, grad, gradReal, n1, h, h1, h2);
  while (true) {
    arma::vec first = beta - gradReal / phiNew;
    arma::vec second = Lambda / phiNew;
    betaNew = softThresh(first, second, p);
    double fVal = lossGauss(Z, censor, Y, betaNew, h, h1, h2);
    arma::vec diff = betaNew - beta;
    double psiVal = loss + arma::as_scalar(grad.t() * diff) + 0.5 * phiNew * arma::as_scalar(diff.t() * diff);
    if (fVal <= psiVal) {
      break;
    }
    phiNew *= gamma;
  }
  beta = betaNew;
  return phiNew;
}

// [[Rcpp::export]]
arma::vec lasso(const arma::mat& Z, const arma::vec& Y, const double lambda, const double tau, const int p, const double n1, const double phi0 = 0.1, 
                const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammL2(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec scad(const arma::mat& Z, const arma::vec& Y, const double lambda, const double tau, const int p, const double n1, const double phi0 = 0.1, 
               const double gamma = 1.2, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammL2(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 0;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 50) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaSCAD(beta, lambda, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammL2(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
      phi = std::max(phi0, phi / gamma);
      if (arma::norm(betaNew - beta, "inf") <= epsilon) {
        break;
      }
      beta = betaNew;
    }
    if (arma::norm(betaNew - beta0, "inf") <= epsilon) {
      break;
    }
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec mcp(const arma::mat& Z, const arma::vec& Y, const double lambda, const double tau, const int p, const double n1, const double phi0 = 0.1, 
              const double gamma = 1.2, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammL2(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 0;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 50) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaMCP(beta, lambda, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammL2(Z, Y, Lambda, betaNew, tau, phi, gamma, p, n1);
      phi = std::max(phi0, phi / gamma);
      if (arma::norm(betaNew - beta, "inf") <= epsilon) {
        break;
      }
      beta = betaNew;
    }
    if (arma::norm(betaNew - beta0, "inf") <= epsilon) {
      break;
    }
  }
  return betaNew;
}

// SCQR-Lasso, SCAD and MCP with particular tau and lambda. This is NOT the function for the QR process.
// [[Rcpp::export]]
arma::vec sqr0Lasso(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const double tau, 
                    const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.1, 
                    const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  arma::vec beta = lasso(Z, Y, lambda, tau, p, n1, phi0, gamma, epsilon, iteMax);
  arma::vec quant = {tau};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  arma::vec betaNew = beta;
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec sqrkLasso(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, arma::vec beta, 
                    const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.1, 
                    const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  arma::vec betaNew = beta;
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec sqr0Scad(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const double tau, 
                   const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.1, 
                   const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  arma::vec beta = lasso(Z, Y, lambda, tau, p, n1, phi0, gamma, epsilon, iteMax);
  arma::vec quant = {tau};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  arma::vec betaNew = beta;
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 1;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 3) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaSCAD(beta, lambda, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
      phi = std::max(phi0, phi / gamma);
      if (arma::norm(betaNew - beta, "inf") <= epsilon) {
        break;
      }
      beta = betaNew;
    }
    if (arma::norm(betaNew - beta0, "inf") <= epsilon) {
      break;
    }
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec sqrkScad(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, arma::vec beta, 
                   const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.1, 
                   const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  arma::vec betaNew = beta;
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 1;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 3) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaSCAD(beta, lambda, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
      phi = std::max(phi0, phi / gamma);
      if (arma::norm(betaNew - beta, "inf") <= epsilon) {
        break;
      }
      beta = betaNew;
    }
    if (arma::norm(betaNew - beta0, "inf") <= epsilon) {
      break;
    }
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec sqr0Mcp(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const double tau, 
                  const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.1, 
                  const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  arma::vec beta = lasso(Z, Y, lambda, tau, p, n1, phi0, gamma, epsilon, iteMax);
  arma::vec quant = {tau};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  arma::vec betaNew = beta;
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 1;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 3) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaMCP(beta, lambda, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
      phi = std::max(phi0, phi / gamma);
      if (arma::norm(betaNew - beta, "inf") <= epsilon) {
        break;
      }
      beta = betaNew;
    }
    if (arma::norm(betaNew - beta0, "inf") <= epsilon) {
      break;
    }
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec sqrkMcp(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, arma::vec beta, 
                  const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.1, 
                  const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  arma::vec betaNew = beta;
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 1;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 3) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaMCP(beta, lambda, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, gamma, p, h, n1, h1, h2);
      phi = std::max(phi0, phi / gamma);
      if (arma::norm(betaNew - beta, "inf") <= epsilon) {
        break;
      }
      beta = betaNew;
    }
    if (arma::norm(betaNew - beta0, "inf") <= epsilon) {
      break;
    }
  }
  return betaNew;
}

// SCQR process with a dilating lambda
// [[Rcpp::export]]
arma::mat SqrLasso(const arma::mat& X, const arma::vec& censor, arma::vec Y, const double lambda, const arma::vec& tauSeq, const double h, 
                    const double phi0 = 0.1, const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols;
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::mat betaProc(p + 1, m);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  arma::vec betaHat = sqr0Lasso(Z, censor, Y, lambda, accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  arma::vec HSeq = getH(tauSeq);
  arma::vec dilate = 1 + arma::log((1 - tauSeq(0)) / (1 - tauSeq));
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkLasso(Z, censor, Y, lambda * dilate(k), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}

// [[Rcpp::export]]
Rcpp::List cvSqrLasso(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, const arma::vec& tauSeq, 
                      const int kfolds, const double h, const double phi0 = 0.1, const double gamma = 1.2, const double epsilon = 0.01, 
                      const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols, nlambda = lambdaSeq.size();
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::vec betaHat(p + 1);
  arma::mat betaProc(p + 1, m);
  arma::vec mse = arma::zeros(nlambda);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec HSeq = getH(tauSeq);
  mse = arma::zeros(nlambda);
  for (int j = 1; j <= kfolds; j++) {
    arma::uvec idx = arma::find(folds == j);
    arma::uvec idxComp = arma::find(folds != j);
    int nTrain = idxComp.size();
    int nTest = idx.size();
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * arma::ones(nTrain);
      betaHat = sqr0Lasso(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkLasso(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, betaHat, p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, HSeq, m, nTest);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  betaHat = sqr0Lasso(Z, censor, Y, lambdaSeq(cvIdx), accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkLasso(Z, censor, Y, lambdaSeq(cvIdx), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return Rcpp::List::create(Rcpp::Named("beta") = betaProc, Rcpp::Named("lambda0") = lambdaSeq(cvIdx));
}

// [[Rcpp::export]]
Rcpp::List cvSqrLassoGrow(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, 
                          const arma::vec& tauSeq, const int kfolds, const double h, const double phi0 = 0.1, const double gamma = 1.2, 
                          const double epsilon = 0.01, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols, nlambda = lambdaSeq.size();
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::vec betaHat(p + 1);
  arma::mat betaProc(p + 1, m);
  arma::vec mse = arma::zeros(nlambda);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec HSeq = getH(tauSeq);
  mse = arma::zeros(nlambda);
  arma::vec dilate = 1 + arma::log((1 - tauSeq(0)) / (1 - tauSeq));
  for (int j = 1; j <= kfolds; j++) {
    arma::uvec idx = arma::find(folds == j);
    arma::uvec idxComp = arma::find(folds != j);
    int nTrain = idxComp.size();
    int nTest = idx.size();
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * arma::ones(nTrain);
      double lambda = lambdaSeq(i);
      betaHat = sqr0Lasso(trainZ, trainCensor, trainY, lambda, trainAccu, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkLasso(trainZ, trainCensor, trainY, lambda * dilate(k), trainAccu, betaHat, p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, HSeq, m, nTest);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  double lambda = lambdaSeq(cvIdx);
  betaHat = sqr0Lasso(Z, censor, Y, lambda, accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkLasso(Z, censor, Y, lambda * dilate(k), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return Rcpp::List::create(Rcpp::Named("beta") = betaProc, Rcpp::Named("lambda0") = lambdaSeq(cvIdx));
}

// [[Rcpp::export]]
Rcpp::List cvSqrLassoIncr(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& incrSeq, const arma::vec& folds, 
                          const arma::vec& tauSeq, const int kfolds, const double h, const double phi0 = 0.1, const double gamma = 1.2, 
                          const double epsilon = 0.01, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols, nlambda = lambdaSeq.size(), nincr = incrSeq.size();
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::vec betaHat(p + 1);
  arma::mat betaProc(p + 1, m);
  arma::vec mse = arma::zeros(nlambda, nincr);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec HSeq = getH(tauSeq);
  mse = arma::zeros(nlambda);
  for (int j = 1; j <= kfolds; j++) {
    arma::uvec idx = arma::find(folds == j);
    arma::uvec idxComp = arma::find(folds != j);
    int nTrain = idxComp.size();
    int nTest = idx.size();
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      double lambda = lambdaSeq(i);
      for (int l = 0; l < nincr; l++) {
        double incr = incrSeq(l);
        arma::vec trainAccu = tauSeq(0) * arma::ones(nTrain);
        betaHat = sqr0Lasso(trainZ, trainCensor, trainY, lambda, trainAccu, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(0) = betaHat;
        for (int k = 1; k < m; k++) {
          arma::vec trainRes = trainY - trainZ * betaHat;
          trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
          betaHat = sqrkLasso(trainZ, trainCensor, trainY, lambda + incr * k, trainAccu, betaHat, p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
          betaProc.col(k) =  betaHat;
        }
        mse(i, l) += calRes(testZ, testCensor, testY, betaProc, tauSeq, HSeq, m, nTest);
      }
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  double lambda = lambdaSeq((int)cvIdx % nlambda);
  double incr = incrSeq((int)cvIdx / nlambda);
  betaHat = sqr0Lasso(Z, censor, Y, lambda, accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkLasso(Z, censor, Y, lambda + incr * k, accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return Rcpp::List::create(Rcpp::Named("beta") = betaProc, Rcpp::Named("lambda0") = lambda, cpp::Named("incr0") = incr);
}

// [[Rcpp::export]]
int exam(arma::mat A) {
  arma::uword cvIdx = A.index_min();
  return (int)cvIdx % 3;
}
  
// [[Rcpp::export]]
arma::mat SqrScad(const arma::mat& X, const arma::vec& censor, arma::vec Y, const double lambda, const arma::vec& tauSeq, const double h, 
                  const double phi0 = 0.1, const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols;
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::mat betaProc(p + 1, m);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  arma::vec betaHat = sqr0Scad(Z, censor, Y, lambda, accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  arma::vec HSeq = getH(tauSeq);
  arma::vec dilate = 1 + arma::log((1 - tauSeq(0)) / (1 - tauSeq));
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkScad(Z, censor, Y, lambda * dilate(k), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}


// [[Rcpp::export]]
arma::mat cvSqrScad(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, const arma::vec& tauSeq, 
                    const int kfolds, const double h, const double phi0 = 0.1, const double gamma = 1.2, const double epsilon = 0.01, 
                    const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols, nlambda = lambdaSeq.size();
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::vec betaHat(p + 1);
  arma::mat betaProc(p + 1, m);
  arma::vec mse = arma::zeros(nlambda);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec HSeq = getH(tauSeq);
  mse = arma::zeros(nlambda);
  for (int j = 1; j <= kfolds; j++) {
    arma::uvec idx = arma::find(folds == j);
    arma::uvec idxComp = arma::find(folds != j);
    int nTrain = idxComp.size();
    int nTest = idx.size();
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * arma::ones(nTrain);
      betaHat = sqr0Scad(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkScad(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, betaHat, p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, HSeq, m, nTest);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  betaHat = sqr0Scad(Z, censor, Y, lambdaSeq(cvIdx), accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkScad(Z, censor, Y, lambdaSeq(cvIdx), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}

// [[Rcpp::export]]
arma::mat cvSqrScadGrow(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, 
                        const arma::vec& tauSeq, const int kfolds, const double h, const double phi0 = 0.1, const double gamma = 1.2, 
                        const double epsilon = 0.01, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols, nlambda = lambdaSeq.size();
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::vec betaHat(p + 1);
  arma::mat betaProc(p + 1, m);
  arma::vec mse = arma::zeros(nlambda);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec HSeq = getH(tauSeq);
  mse = arma::zeros(nlambda);
  arma::vec dilate = 1 + arma::log((1 - tauSeq(0)) / (1 - tauSeq));
  for (int j = 1; j <= kfolds; j++) {
    arma::uvec idx = arma::find(folds == j);
    arma::uvec idxComp = arma::find(folds != j);
    int nTrain = idxComp.size();
    int nTest = idx.size();
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * arma::ones(nTrain);
      double lambda = lambdaSeq(i);
      betaHat = sqr0Scad(trainZ, trainCensor, trainY, lambda, trainAccu, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkScad(trainZ, trainCensor, trainY, lambda * dilate(k), trainAccu, betaHat, p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, HSeq, m, nTest);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  double lambda = lambdaSeq(cvIdx);
  betaHat = sqr0Scad(Z, censor, Y, lambda, accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkScad(Z, censor, Y, lambda * dilate(k), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}

// [[Rcpp::export]]
arma::mat SqrMcp(const arma::mat& X, const arma::vec& censor, arma::vec Y, const double lambda, const arma::vec& tauSeq, const double h, 
                 const double phi0 = 0.1, const double gamma = 1.2, const double epsilon = 0.01, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols;
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::mat betaProc(p + 1, m);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  arma::vec betaHat = sqr0Mcp(Z, censor, Y, lambda, accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  arma::vec HSeq = getH(tauSeq);
  arma::vec dilate = 1 + arma::log((1 - tauSeq(0)) / (1 - tauSeq));
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkMcp(Z, censor, Y, lambda * dilate(k), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}


// [[Rcpp::export]]
arma::mat cvSqrMcp(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, const arma::vec& tauSeq, 
                   const int kfolds, const double h, const double phi0 = 0.1, const double gamma = 1.2, const double epsilon = 0.01, 
                   const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols, nlambda = lambdaSeq.size();
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::vec betaHat(p + 1);
  arma::mat betaProc(p + 1, m);
  arma::vec mse = arma::zeros(nlambda);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec HSeq = getH(tauSeq);
  mse = arma::zeros(nlambda);
  for (int j = 1; j <= kfolds; j++) {
    arma::uvec idx = arma::find(folds == j);
    arma::uvec idxComp = arma::find(folds != j);
    int nTrain = idxComp.size();
    int nTest = idx.size();
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * arma::ones(nTrain);
      betaHat = sqr0Mcp(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkMcp(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, betaHat, p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, HSeq, m, nTest);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  betaHat = sqr0Mcp(Z, censor, Y, lambdaSeq(cvIdx), accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkMcp(Z, censor, Y, lambdaSeq(cvIdx), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}

// [[Rcpp::export]]
arma::mat cvSqrMcpGrow(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, 
                       const arma::vec& tauSeq, const int kfolds, const double h, const double phi0 = 0.1, const double gamma = 1.2, 
                       const double epsilon = 0.01, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols, nlambda = lambdaSeq.size();
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::vec betaHat(p + 1);
  arma::mat betaProc(p + 1, m);
  arma::vec mse = arma::zeros(nlambda);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec HSeq = getH(tauSeq);
  mse = arma::zeros(nlambda);
  arma::vec dilate = 1 + arma::log((1 - tauSeq(0)) / (1 - tauSeq));
  for (int j = 1; j <= kfolds; j++) {
    arma::uvec idx = arma::find(folds == j);
    arma::uvec idxComp = arma::find(folds != j);
    int nTrain = idxComp.size();
    int nTest = idx.size();
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * arma::ones(nTrain);
      double lambda = lambdaSeq(i);
      betaHat = sqr0Mcp(trainZ, trainCensor, trainY, lambda, trainAccu, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkMcp(trainZ, trainCensor, trainY, lambda * dilate(k), trainAccu, betaHat, p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, HSeq, m, nTest);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * arma::ones(n);
  double lambda = lambdaSeq(cvIdx);
  betaHat = sqr0Mcp(Z, censor, Y, lambda, accu, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkMcp(Z, censor, Y, lambda * dilate(k), accu, betaHat, p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}


