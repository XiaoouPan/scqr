###### Simulation code for scqr
library(quantreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(survival)
library(caret)
library(rqPen)
library(hqreg)
library(tikzDevice)
library(ggplot2)

rm(list = ls())
Rcpp::sourceCpp("src/hdscqr.cpp")

exam = function(trueSig, selectSet) {
  p = length(trueSig)
  selectSig = rep(0, p)
  selectSig[selectSet] = 1
  TPR = sum(trueSig != 0 & selectSig != 0) / sum(trueSig != 0)
  TNR = sum(trueSig == 0 & selectSig == 0) / sum(trueSig == 0)
  FDR = 0
  if (sum(selectSig != 0) > 0) {
    FDR = sum(trueSig == 0 & selectSig != 0) / sum(selectSig != 0)
  }
  return (list("TPR" = TPR, "TNR" = TNR, "FDR" = FDR))
}

getSet = function(beta.hat, m) {
  active = 1 * (beta.hat[-1, ] != 0)
  uniActive = which(rowMaxs(active) != 0)
  voteActive = which(rowSums(active) > 0.5 * m)
  return (list("union" = uniActive, "vote" = voteActive))
}

H = function(x) {
  return (-log(1 - x))
}

## High-dim CQR-Lasso, modified based on Zheng, Peng and He (2018)
quantproc = function(y, x, delta, JJ, lambda, tol = 1e-4) {
  pp = dim(x)[2]                                                       
  tmpbeta = matrix(0, length(JJ), pp)
  Hseq = H(JJ)
  dilate = 1 + log((1 - JJ[1]) / (1 - JJ))
  
  augy1 = y[which(delta == 1)]                                           
  augx1 = x[which(delta == 1), ]                                          
  augx2 = -colSums(augx1)                                                                 
  augy = c(augy1, 1000, 1000)
  
  sto_weights = matrix(0, length(JJ), dim(x)[1])                   
  rproc = matrix(0, length(JJ), dim(x)[1])                           
  sto_weights[1, ] = JJ[1] * 2                                            
  augx3 = sto_weights[1, ] %*% x                                       
  augx = rbind(augx1, augx2, augx3)
  beta0 = LASSO.fit(augy, augx, tau = 0.5, lambda = lambda, intercept = FALSE, coef.cutoff = 0.0001, weights = NULL)
  tmpbeta[1, ] = beta0
  rproc[1, ] = 1 * (y >= x %*% beta0)                                         
  for(s in 2:length(JJ)) {
    Hm = Hseq[s] - Hseq[s - 1]                                 
    sto_weights[s, ] = sto_weights[s - 1, ] + 2 * Hm * rproc[s - 1, ]           
    augx3 = sto_weights[s, ] %*% x                                    
    augx = rbind(augx1, augx2, augx3)
    tmpbeta[s,] = LASSO.fit(augy, augx, tau = 0.5, lambda = lambda * dilate[s], intercept = FALSE, coef.cutoff = 0.0001, weights = NULL)
    rproc[s, ] = 1*( y > x %*% tmpbeta[s, ])                              
  }
  return (t(tmpbeta))
}


#### High-dim quantile process with fixed scale
n = 40
p = 100
s = 2
M = 1
kfolds = 3
h = (log(p) / n)^(1/4)
tauSeq = seq(0.1, 0.7, by = 0.05)
grid = seq(0.1, 0.75, by = 0.05)
m = length(tauSeq)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
Sigma = toeplitz(0.5^(0:(p - 1)))
lambdaSeq = exp(seq(log(0.01), log(0.2), length.out = 50))
trueSig = c(rep(1, s), rep(0, p - s))

time1 = time2 = prop = rep(0, M)
TPR1 = TNR1 = FDR1 = error1 = rep(NA, M)
TPR2 = TNR2 = FDR2 = error2 = rep(NA, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  X = mvrnorm(n, rep(0, p), Sigma)
  err = rt(n, 2)
  ## Homo
  beta = c(runif(s, 1, 1.5), rep(0, p - s))
  betaMat = rbind(beta0, matrix(beta, p, nTau))
  logT = X %*% beta + err
  ## Hetero
  #X[, 1] = abs(X[, 1])
  #beta = c(runif(s - 1, 1, 1.5), rep(0, p - s))
  #betaMat = rbind(rep(0, nTau), beta0, matrix(beta, p - 1, nTau))
  #logT = X[, 1] * err + X[, -1] %*% beta
  
  w = sample(1:3, n, prob = c(1/3, 1/3, 1/3), replace = TRUE)
  logC = (w == 1) * rnorm(n, 0, 4) + (w == 2) * rnorm(n, 5, 1) + (w == 3) * rnorm(n, 10, 0.5)
  censor = logT <= logC
  prop[i] = 1 - sum(censor) / n
  Y = pmin(logT, logC)
  folds = createFolds(censor, kfolds, FALSE)
  
  ## SCQR-Lasso
  scqr.fit = cvSqrLasso(X, censor, Y, lambdaSeq, folds, tauSeq, kfolds, h)
  lambda0 = scqr.fit$lambda0
  start = Sys.time()
  beta.lasso = SqrLasso(X, censor, Y, lambda0, tauSeq, h)
  end = Sys.time()
  time1[i] = as.numeric(difftime(end, start, units = "secs"))
  

  ## HDCQR-Lasso using rqPen
  Z = cbind(1, X)
  start = Sys.time()
  beta.cqr = quantproc(Y, Z, censor, tauSeq, lambda0)
  end = Sys.time()
  time2[i] = as.numeric(difftime(end, start, units = "secs"))
  activeSet = getSet(beta.cqr, m)
  uniSet = activeSet$union
  voteSet = activeSet$vote
  ## scqr on the union set
  if (length(uniSet) > 0) {
    test = exam(trueSig, uniSet)
    TPR1[i] = test$TPR
    TNR1[i] = test$TNR
    FDR1[i] = test$FDR
  }
  ## scqr on the majority vote set
  if (length(voteSet) > 0) {
    test = exam(trueSig, voteSet)
    TPR2[i] = test$TPR
    TNR2[i] = test$TNR
    FDR2[i] = test$FDR
  }

  setTxtProgressBar(pb, i / M)
}


setwd("~/Dropbox/Conquer/SCQR/code/Simulation/highd/hetero/lasso")
tp.cqr = as.matrix(read.csv("tp_cqr.csv")[, -1])
mtc.cqr = as.matrix(read.csv("mtc_cqr.csv")[, -1])
setwd("~/Dropbox/Conquer/SCQR/code/Simulation/highd/hetero")
mtc.lasso = as.matrix(read.csv("mtc_lasso.csv")[, -1])


M = 500
### Data for plots
meth = c(rep("CV-SCQR-Lasso", M), rep("CQR-Lasso", M))
meth = factor(meth, levels = c("CV-SCQR-Lasso", "CQR-Lasso"))
time = c(tp.cqr[1, ], tp.cqr[2, ])
prop = c(tp.cqr[3, ], tp.cqr[3, ])
TPR = c(mtc.lasso[1, ], mtc.cqr[1, ])
TNR = c(mtc.lasso[2, ], mtc.cqr[2, ])
FDR = c(mtc.lasso[3, ], mtc.cqr[3, ])
dat = data.frame("TPR" = TPR, "TNR" = TNR, "FDR" = FDR, "time" = time, "prop" = prop, "method" = meth)



### TPR
setwd("~/Dropbox/Conquer/SCQR/code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = TPR, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Estimation") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



### FDR
setwd("~/Dropbox/Conquer/SCQR/code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = FDR, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("False discovery rate") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### error
setwd("~/Dropbox/Conquer/SCQR/code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = error, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("$$Estimation error in $||\\cdot||_2$") + ylim(0.4, 2) + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### time
setwd("~/Dropbox/Conquer/SCQR/code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = time, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.8, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Elapsed time (in seconds)") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



