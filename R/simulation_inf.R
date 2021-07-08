###### Simulation code for scqr

rm(list = ls())
Rcpp::sourceCpp("src/scqr.cpp")

library(quantreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(survival)
library(ggplot2)
library(tikzDevice)

getSigma = function(p) {
  sig = diag(p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      sig[i, j] = sig[j, i] = 0.5^(j - i)
    }
  }
  return (sig)
}

estError = function(betahat, beta, tauSeq) {
  nTau = min(length(tauSeq), ncol(betahat) + 1)
  diff = betahat - beta
  err = sqrt(colSums((betahat - beta)^2))
  accu = 0
  for (k in 2:nTau) {
    accu = accu + err[k - 1] * (tauSeq[k] - tauSeq[k - 1])
  }
  return (accu)
}

getPivCI = function(est, estBoot, alpha) {
  q1 = rowQuantiles(estBoot, probs = alpha / 2)
  q2 = rowQuantiles(estBoot, probs = 1 - alpha / 2)
  perCI = cbind(q1, q2)
  pivCI = cbind(2 * est - q2, 2 * est - q1)
  colnames(perCI) = colnames(pivCI) = c("lower", "upper")
  return (list(perCI = perCI, pivCI = pivCI))
}

getNormCI = function(est, sd, z) {
  lower = est - z * sd
  upper = est + z * sd
  return (cbind(lower, upper))
}

getCoverage = function(ci, beta) {
  return (ci[-1, 1] <  beta & ci[-1, 2] > beta)
}

getWidth = function(ci) {
  return (ci[-1, 2] - ci[-1, 1])
}


getCoverPlot = function(cover1, cover2, cover3, p) {
  rst1 = colMeans(cover1)
  rst2 = colMeans(cover2)
  rst3 = colMeans(cover3)
  meth = c(rep("per", p), rep("piv", p), rep("norm", p))
  meth = factor(meth, levels = c("per", "piv", "norm"))
  rst = data.frame("cover" = c(rst1, rst2, rst3), "method" = meth)
  return (rst)
}

getWidthPlot = function(width1, width2, width3, j, M) {
  meth = c(rep("per", M), rep("piv", M), rep("norm", M))
  meth = factor(meth, levels = c("per", "piv", "norm"))
  rst = data.frame("width" = c(width1[, j], width2[, j], width3[, j]), "method" = meth)
  return (rst)
}


#### Fixed scale,
n = 2000
p = n / 50
M = 1000
tauSeq = seq(0.1, 0.5, by = 0.1)
grid = seq(0.1, 0.5, by = 0.1)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
alpha = 0.05
z = qnorm(1 - alpha / 2)
cover1 = cover2 = cover3 = matrix(0, M, p)
width1 = width2 = width3 = matrix(0, M, p)
time = prop = rep(0, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  #X = sqrt(12) * draw.d.variate.uniform(n, p, Sigma) - sqrt(3)
  Sigma = getSigma(p)
  X = mvrnorm(n, rep(0, p), Sigma)
  #Sigma = getSigma(45)
  #X = cbind(mvrnorm(n, rep(0, 45), Sigma), 4 * draw.d.variate.uniform(n, 45, Sigma) - 2, matrix(rbinom(10 * n, 1, c(0.5, 0.5)), n, 10))
  err = rt(n, 2)
  ## Homo
  beta = runif(p, -2, 2)
  betaMat = rbind(beta0, matrix(beta, p, nTau))
  logT = X %*% beta + err
  ## Hetero
  #X[, 1] = abs(X[, 1])
  #beta = runif(p - 1, -2, 2)
  #betaMat = rbind(rep(0, nTau), beta0, matrix(beta, p - 1, nTau))
  #logT = X[, 1] * err + X[, -1] %*% beta
  w = sample(1:3, n, prob = c(1/3, 1/3, 1/3), replace = TRUE)
  logC = (w == 1) * rnorm(n, 0, 4) + (w == 2) * rnorm(n, 5, 1) + (w == 3) * rnorm(n, 10, 0.5)
  censor = logT <= logC
  prop[i] = 1 - sum(censor) / n
  Y = pmin(logT, logC)
  response = Surv(Y, censor, type = "right")
  
  start = Sys.time()
  list = scqrGaussInf(X, Y, censor, tauSeq, B = 1000)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  beta.hat = list$coeff
  beta.boot = list$boot
  ci.list = getPivCI(beta.hat[, nTau], beta.boot, alpha)
  ci.per = ci.list$perCI
  ci.piv = ci.list$pivCI
  ci.norm = getNormCI(beta.hat[, nTau], rowSds(beta.boot), z)
  cover1[i, ] = getCoverage(ci.per, beta)
  cover2[i, ] = getCoverage(ci.piv, beta)
  cover3[i, ] = getCoverage(ci.norm, beta)
  width1[i, ] = getWidth(ci.per)
  width2[i, ] = getWidth(ci.piv)
  width3[i, ] = getWidth(ci.norm)
  
  setTxtProgressBar(pb, i / M)
}


mean(time)

write.csv(prop, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/prop.csv")
write.csv(time, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/time.csv")
write.csv(cover1, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/cover1.csv")
write.csv(cover2, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/cover2.csv")
write.csv(cover3, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/cover3.csv")
write.csv(width1, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/width1.csv")
write.csv(width2, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/width2.csv")
write.csv(width3, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/width3.csv")

summ = getCoverPlot(cover1, cover2, cover3, p)
setwd("~/Dropbox/Conquer/censoredQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(summ, aes(x = method, y = cover, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Coverage") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


summ = getWidthPlot(width1, width2, width3, 1, M)
setwd("~/Dropbox/Conquer/censoredQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(summ, aes(x = method, y = width, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("CI Width") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



