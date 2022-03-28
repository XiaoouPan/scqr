###### Cross validation for the increment in lambda sequence, suggested by a referee. This is for the comparison with Zheng, Peng and He (2018)
library(quantreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(survival)
library(caret)
library(rqPen)
library(tikzDevice)
library(ggplot2)

rm(list = ls())
Rcpp::sourceCpp("src/hdscqr.cpp")

exam = function(trueSig, selectSet, beta.est, beta) {
  m = ncol(beta)
  p = length(trueSig)
  beta.hat = matrix(0, p, m)
  beta.hat[selectSet, ] = beta.est[-1, ]
  selectSig = rep(0, p)
  selectSig[selectSet] = 1
  TPR = sum(trueSig != 0 & selectSig != 0) / sum(trueSig != 0)
  TNR = sum(trueSig == 0 & selectSig == 0) / sum(trueSig == 0)
  FDR = 0
  if (sum(selectSig != 0) > 0) {
    FDR = sum(trueSig == 0 & selectSig != 0) / sum(selectSig != 0)
  }
  err = 0
  for (i in 1:m) {
    err = err + norm(beta[, i] - beta.hat[, i], "2")
  }
  return (list("TPR" = TPR, "TNR" = TNR, "FDR" = FDR, "error" = err / m))
}

getSet = function(beta.hat, m) {
  active = 1 * (beta.hat[-1, ] != 0)
  uniActive = which(rowMaxs(active) != 0)
  voteActive = which(rowSums(active) > 0.5 * m)
  return (list("union" = uniActive, "vote" = voteActive))
}


#### Quantile process with fixed scale, hard to visualize
n = 200
p = 500
s = 5
M = 100
kfolds = 3
h = 0.5 * (log(p) / n)^(1/4)
tauSeq = seq(0.1, 0.7, by = 0.05)
m = length(tauSeq)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
Sigma = toeplitz(0.5^(0:(p - 1)))
lambdaSeq = exp(seq(log(0.01), log(0.2), length.out = 50))
incrSeq = seq(0.001, 0.02, length.out = 20)
trueSig = c(rep(1, s), rep(0, p - s))

TPR1 = TNR1 = FDR1 = error1 = time1 = rep(NA, M)
TPR2 = TNR2 = FDR2 = error2 = time2 = rep(NA, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  X = mvrnorm(n, rep(0, p), Sigma)
  err = rt(n, 2)
  ## Homo
  #beta = c(runif(s, 1, 1.5), rep(0, p - s))
  #betaMat = rbind(beta0, matrix(beta, p, nTau))
  #logT = X %*% beta + err
  ## Hetero
  X[, 1] = abs(X[, 1])
  beta = c(runif(s - 1, 1, 1.5), rep(0, p - s))
  betaMat = rbind(rep(0, nTau), beta0, matrix(beta, p - 1, nTau))
  logT = X[, 1] * err + X[, -1] %*% beta
  
  w = sample(1:3, n, prob = c(1/3, 1/3, 1/3), replace = TRUE)
  logC = (w == 1) * rnorm(n, 0, 4) + (w == 2) * rnorm(n, 5, 1) + (w == 3) * rnorm(n, 10, 0.5)
  censor = logT <= logC
  Y = pmin(logT, logC)
  folds = createFolds(censor, kfolds, FALSE)
  
  ## SCQR-Lasso
  start = Sys.time()
  fit = cvSqrLassoGrow(X, censor, Y, lambdaSeq, folds, tauSeq, kfolds, h)
  end = Sys.time()
  beta.lasso = fit$beta
  time1[i] = as.numeric(difftime(end, start, units = "secs"))
  activeSet = getSet(beta.lasso, m)
  uniSet = activeSet$union
  Xunion = X[, uniSet, drop = FALSE]
  ## scqr on the union set
  if (length(uniSet) > 0) {
    beta.union = scqrGauss(Xunion, Y, censor, tauSeq)
    test = exam(trueSig, uniSet, beta.union, betaMat[-1, ])
    TPR1[i] = test$TPR
    TNR1[i] = test$TNR
    FDR1[i] = test$FDR
    error1[i] = test$error
  }
  
  ## SCQR-Lasso with cv for increment
  start = Sys.time()
  fit = cvSqrLassoIncr(X, censor, Y, lambdaSeq, incrSeq, folds, tauSeq, kfolds, h)
  end = Sys.time()
  beta.lasso2 = fit$beta
  time2[i] = as.numeric(difftime(end, start, units = "secs"))
  activeSet = getSet(beta.lasso2, m)
  uniSet = activeSet$union
  Xunion = X[, uniSet, drop = FALSE]
  ## scqr on the union set
  if (length(uniSet) > 0) {
    beta.union = scqrGauss(Xunion, Y, censor, tauSeq)
    test = exam(trueSig, uniSet, beta.union, betaMat[-1, ])
    TPR2[i] = test$TPR
    TNR2[i] = test$TNR
    FDR2[i] = test$FDR
    error2[i] = test$error
  }

  setTxtProgressBar(pb, i / M)
}


result1 = rbind(TPR1, TNR1, FDR1, error1, time)
result2 = rbind(TPR2, TNR2, FDR2, error2, time)
write.csv(result1, "~/Dropbox/Conquer/SCQR/AOS_rev/Simulation/Highd/result1_hetero.csv")
write.csv(result2, "~/Dropbox/Conquer/SCQR/AOS_rev/Simulation/Highd/result2_hetero.csv")



setwd("~/Dropbox/Conquer/SCQR/Code/Simulation/highd/hetero")
tp.lasso = as.matrix(read.csv("tp_lasso.csv")[, -1])
tp.scad = as.matrix(read.csv("tp_scad.csv")[, -1])
tp.mcp = as.matrix(read.csv("tp_mcp.csv")[, -1])
mtc.lasso = as.matrix(read.csv("mtc_lasso.csv")[, -1])
mtc.scad = as.matrix(read.csv("mtc_scad.csv")[, -1])
mtc.mcp = as.matrix(read.csv("mtc_mcp.csv")[, -1])



### Data for plots
meth = c(rep("Dilating $\\lambda$", M), rep("CV for $\\lambda$ increment", M))
meth = factor(meth, levels = c("Dilating $\\lambda$", "CV for $\\lambda$ increment"))
time = c(time1, time2)
TPR = c(TPR1, TPR2)
TNR = c(TNR1, TNR2)
FDR = c(FDR1, FDR2)
error = c(error1, error2)
dat = data.frame("TPR" = TPR, "TNR" = TNR, "FDR" = FDR, "error" = error, "time" = time, "method" = meth)


### TPR
setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = TPR, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Estimation") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### TNR
setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = TNR, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("True negative rate") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### FDR
setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = FDR, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("False discovery rate") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### error
setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = error, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("$$Estimation error in $||\\cdot||_2$") + ylim(0.4, 2) + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### time
setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = time, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Elapsed time (in seconds)") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

