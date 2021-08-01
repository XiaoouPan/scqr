# include <RcppArmadillo.h>
# include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
int sgn(const double x) {
  return (x > 0) - (x < 0);
}

// [[Rcpp::export]]
arma::vec getH(const arma::vec& tauSeq) {
  return -arma::log(1 - tauSeq);
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
arma::vec cmptLambdaLasso(const double lambda, const arma::vec& sx1, const int p) {
  arma::vec rst = arma::zeros(p + 1);
  rst.rows(1, p) = lambda * sx1;
  return rst;
}

// [[Rcpp::export]]
arma::vec cmptLambdaSCAD(const arma::vec& beta, const double lambda, const arma::vec& sx1, const int p) {
  arma::vec rst(p + 1);
  rst(0) = 0;
  for (int i = 1; i <= p; i++) {
    double abBeta = std::abs(beta(i));
    double lambdaSx = lambda * sx1(i - 1);
    if (abBeta <= lambdaSx) {
      rst(i) = lambdaSx;
    } else if (abBeta <= 3.7 * lambdaSx) {
      rst(i) = 0.37037 * (3.7 * lambdaSx - abBeta);
    } else {
      rst(i) = 0;
    }
  }
  return rst;
}

// [[Rcpp::export]]
arma::vec cmptLambdaMCP(const arma::vec& beta, const double lambda, const arma::vec& sx1, const int p) {
  arma::vec rst(p + 1);
  rst(0) = 0;
  for (int i = 1; i <= p; i++) {
    double abBeta = std::abs(beta(i));
    double lambdaSx = lambda * sx1(i - 1);
    rst(i) = (abBeta <= 3 * lambdaSx) ? (lambdaSx - 0.33333 * abBeta) : 0;
  }
  return rst;
}

// [[Rcpp::export]]
double lossL2(const arma::mat& Z, const arma::vec& Y, const arma::vec& beta, const double n1, const double tau) {
  arma::vec res = Z * beta - Y;
  double rst = 0.0;
  for (int i = 0; i < Y.size(); i++) {
    rst += (res(i) > 0) ? (tau * res(i) * res(i)) : ((1 - tau) * res(i) * res(i));
  }
  return 0.5 * n1 * rst;
}

// [[Rcpp::export]]
double updateL2(const arma::mat& Z, const arma::vec& Y, const arma::vec& beta, arma::vec& grad, const double n1, const double tau) {
  arma::vec res = Z * beta - Y;
  double rst = 0.0;
  grad = arma::zeros(grad.size());
  for (int i = 0; i < Y.size(); i++) {
    double temp = res(i) > 0 ? tau : (1 - tau);
    grad += temp * res(i) * Z.row(i).t();
    rst += temp * res(i) * res(i);
  }
  grad *= n1;
  return 0.5 * n1 * rst;
}

// [[Rcpp::export]]
double calRes(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const arma::mat& betaHat, const arma::vec& tauSeq, const int m) {
  double res = 0.0;
  arma::uvec idx = arma::find(censor == 1);
  arma::mat Zcen = Z.rows(idx);
  arma::mat Ycen = Y.rows(idx);
  int ncen = Ycen.size();
  for (int k = 0; k < m; k++) {
    arma::vec resCen = Ycen - Zcen * betaHat.col(k);
    double temp = 0.0;
    for (int i = 0; i < ncen; i++) {
      temp += resCen(i) >= 0 ? (resCen(i) * tauSeq(k)) : (resCen(i) * (tauSeq(k) - 1));
    }
    res += temp / ncen;
  }
  return res / m;
}

// [[Rcpp::export]]
double lossGauss(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const arma::vec& accu, const arma::vec& beta, const double tau, 
                 const double h, const double h1, const double h2) {
  arma::vec res = Z * beta - Y;
  arma::vec temp = 0.39894 * h  * arma::exp(-0.5 * h2 * arma::square(res)) - tau * res + res % arma::normcdf(h1 * res);
  return arma::mean(censor % temp - accu % (Z * beta));
}

// [[Rcpp::export]]
double updateGauss(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const arma::vec& accu, const arma::vec& beta, arma::vec& grad, 
                   const double tau, const double n1, const double h, const double h1, const double h2) {
  arma::vec res = Z * beta - Y;
  arma::vec der = censor % arma::normcdf(res * h1) - accu;
  grad = n1 * Z.t() * der;
  arma::vec temp = 0.39894 * h  * arma::exp(-0.5 * h2 * arma::square(res)) - tau * res + res % arma::normcdf(h1 * res);
  return arma::mean(censor % temp - accu % (Z * beta));
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
              const double phi, const double tau, const double gamma, const int p, const double h, const double n1, const double h1, const double h2) {
  double phiNew = phi;
  arma::vec betaNew(p + 1);
  arma::vec grad(p + 1);
  double loss = updateGauss(Z, censor, Y, accu, beta, grad, tau, n1, h, h1, h2);
  while (true) {
    arma::vec first = beta - grad / phiNew;
    arma::vec second = Lambda / phiNew;
    betaNew = softThresh(first, second, p);
    double fVal = lossGauss(Z, censor, Y, accu, betaNew, tau, h, h1, h2);
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
arma::vec lasso(const arma::mat& Z, const arma::vec& Y, const double lambda, const double tau, const arma::vec& sx1, const int p, const double n1, 
                const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
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
arma::vec scad(const arma::mat& Z, const arma::vec& Y, const double lambda, const double tau, const arma::vec& sx1, const int p, const double n1, 
               const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
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
    Lambda = cmptLambdaSCAD(beta, lambda, sx1, p);
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
arma::vec mcp(const arma::mat& Z, const arma::vec& Y, const double lambda, const double tau, const arma::vec& sx1, const int p, const double n1, 
              const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
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
    Lambda = cmptLambdaMCP(beta, lambda, sx1, p);
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
arma::vec sqr0Lasso(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const arma::vec& sx1, 
                    const double tau, const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.01, 
                    const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = lasso(Z, Y, lambda, tau, sx1, p, n1, phi0, gamma, epsilon, iteMax);
  arma::vec quant = {tau};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  arma::vec betaNew = beta;
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec sqrkLasso(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const arma::vec& sx1, 
                    arma::vec beta, const double tau, const int p, const double n1, const double h, const double h1, const double h2, 
                    const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec betaNew = beta;
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  return betaNew;
}

// [[Rcpp::export]]
arma::vec sqr0Scad(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const arma::vec& sx1, 
                   const double tau, const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.01, 
                   const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = lasso(Z, Y, lambda, tau, sx1, p, n1, phi0, gamma, epsilon, iteMax);
  arma::vec quant = {tau};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  arma::vec betaNew = beta;
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 0;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 1) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaSCAD(beta, lambda, sx1, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
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
arma::vec sqrkScad(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const arma::vec& sx1, 
                   arma::vec beta, const double tau, const int p, const double n1, const double h, const double h1, const double h2, 
                   const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec betaNew = beta;
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 0;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 1) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaSCAD(beta, lambda, sx1, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
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
arma::vec sqr0Mcp(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const arma::vec& sx1, 
                  const double tau, const int p, const double n1, const double h, const double h1, const double h2, const double phi0 = 0.01, 
                  const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec beta = lasso(Z, Y, lambda, tau, sx1, p, n1, phi0, gamma, epsilon, iteMax);
  arma::vec quant = {tau};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  arma::vec betaNew = beta;
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 0;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 5) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaMCP(beta, lambda, sx1, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
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
arma::vec sqrkMcp(const arma::mat& Z, const arma::vec& censor, const arma::vec& Y, const double lambda, const arma::vec& accu, const arma::vec& sx1, 
                  arma::vec beta, const double tau, const int p, const double n1, const double h, const double h1, const double h2, 
                  const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  arma::vec betaNew = beta;
  // Contraction
  arma::vec Lambda = cmptLambdaLasso(lambda, sx1, p);
  double phi = phi0;
  int ite = 0;
  while (ite <= iteMax) {
    ite++;
    phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 0;
  // Tightening
  arma::vec beta0(p + 1);
  while (iteT <= 5) {
    iteT++;
    beta = betaNew;
    beta0 = betaNew;
    Lambda = cmptLambdaMCP(beta, lambda, sx1, p);
    phi = phi0;
    ite = 0;
    while (ite <= iteMax) {
      ite++;
      phi = lammSq(Z, censor, Y, Lambda, accu, betaNew, phi, tau, gamma, p, h, n1, h1, h2);
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

// SCQR process with a particular lambda
// [[Rcpp::export]]
arma::mat SqrLasso(const arma::mat& X, const arma::vec& censor, arma::vec Y, const double lambda, const arma::vec& tauSeq, const double h, 
                   const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols;
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::mat betaProc(p + 1, m);
  arma::vec accu = tauSeq(0) * (1 - censor);
  arma::vec betaHat = sqr0Lasso(Z, censor, Y, lambda, accu, sx1, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  arma::vec HSeq = getH(tauSeq);
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkLasso(Z, censor, Y, lambda, accu, sx1, betaHat, tauSeq(k), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}

// [[Rcpp::export]]
arma::mat cvSqrLasso(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, const arma::vec& tauSeq, 
                     const int kfolds, const double h, const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, 
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
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * (1 - trainCensor);
      betaHat = sqr0Lasso(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, sx1, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkLasso(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, sx1, betaHat, tauSeq(k), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, m);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * (1 - censor);
  betaHat = sqr0Lasso(Z, censor, Y, lambdaSeq(cvIdx), accu, sx1, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkLasso(Z, censor, Y, lambdaSeq(cvIdx), accu, sx1, betaHat, tauSeq(k), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}

// [[Rcpp::export]]
arma::mat SqrScad(const arma::mat& X, const arma::vec& censor, arma::vec Y, const double lambda, const arma::vec& tauSeq, const double h, 
                  const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols;
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::mat betaProc(p + 1, m);
  arma::vec accu = tauSeq(0) * (1 - censor);
  arma::vec betaHat = sqr0Scad(Z, censor, Y, lambda, accu, sx1, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  arma::vec HSeq = getH(tauSeq);
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkScad(Z, censor, Y, lambda, accu, sx1, betaHat, tauSeq(k), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}


// [[Rcpp::export]]
arma::mat cvSqrScad(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, const arma::vec& tauSeq, 
                    const int kfolds, const double h, const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, 
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
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * (1 - trainCensor);
      betaHat = sqr0Scad(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, sx1, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkScad(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, sx1, betaHat, tauSeq(k), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, m);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * (1 - censor);
  betaHat = sqr0Scad(Z, censor, Y, lambdaSeq(cvIdx), accu, sx1, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkScad(Z, censor, Y, lambdaSeq(cvIdx), accu, sx1, betaHat, tauSeq(k), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}

// [[Rcpp::export]]
arma::mat SqrMcp(const arma::mat& X, const arma::vec& censor, arma::vec Y, const double lambda, const arma::vec& tauSeq, const double h, 
                 const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, const int iteMax = 500) {
  const int n = X.n_rows, p = X.n_cols;
  const int m = tauSeq.size();
  const double h1 = 1.0 / h, h2 = 1.0 / (h * h);
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::mat betaProc(p + 1, m);
  arma::vec accu = tauSeq(0) * (1 - censor);
  arma::vec betaHat = sqr0Mcp(Z, censor, Y, lambda, accu, sx1, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  arma::vec HSeq = getH(tauSeq);
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkMcp(Z, censor, Y, lambda, accu, sx1, betaHat, tauSeq(k), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}


// [[Rcpp::export]]
arma::mat cvSqrMcp(const arma::mat& X, const arma::vec& censor, arma::vec Y, const arma::vec& lambdaSeq, const arma::vec& folds, const arma::vec& tauSeq, 
                   const int kfolds, const double h, const double phi0 = 0.01, const double gamma = 1.5, const double epsilon = 0.001, 
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
    double n1Train = 1.0 / idxComp.size();
    arma::mat trainZ = Z.rows(idxComp), testZ = Z.rows(idx);
    arma::vec trainY = Y.rows(idxComp), testY = Y.rows(idx);
    arma::vec trainCensor = censor.rows(idxComp), testCensor = censor.rows(idx);
    for (int i = 0; i < nlambda; i++) {
      arma::vec trainAccu = tauSeq(0) * (1 - trainCensor);
      betaHat = sqr0Mcp(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, sx1, tauSeq(0), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
      betaProc.col(0) = betaHat;
      for (int k = 1; k < m; k++) {
        arma::vec trainRes = trainY - trainZ * betaHat;
        trainAccu += arma::normcdf(trainRes * h1) * (HSeq(k) - HSeq(k - 1));
        betaHat = sqrkMcp(trainZ, trainCensor, trainY, lambdaSeq(i), trainAccu, sx1, betaHat, tauSeq(k), p, n1Train, h, h1, h2, phi0, gamma, epsilon, iteMax);
        betaProc.col(k) =  betaHat;
      }
      mse(i) += calRes(testZ, testCensor, testY, betaProc, tauSeq, m);
    }
  }
  arma::uword cvIdx = arma::index_min(mse);
  arma::vec accu = tauSeq(0) * (1 - censor);
  betaHat = sqr0Mcp(Z, censor, Y, lambdaSeq(cvIdx), accu, sx1, tauSeq(0), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
  betaProc.col(0) = betaHat;
  for (int k = 1; k < m; k++) {
    arma::vec res = Y - Z * betaHat;
    accu += arma::normcdf(res * h1) * (HSeq(k) - HSeq(k - 1));
    betaHat = sqrkMcp(Z, censor, Y, lambdaSeq(cvIdx), accu, sx1, betaHat, tauSeq(k), p, 1.0 / n, h, h1, h2, phi0, gamma, epsilon, iteMax);
    betaProc.col(k) =  betaHat;
  }
  betaProc.rows(1, p).each_col() %= sx1;
  betaProc.row(0) += my - mx * betaProc.rows(1, p);
  return betaProc;
}
