###### Simulation code for scqr
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
n = 50
p = 200
s = 1
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

time = prop = rep(0, M)
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
  beta.lasso = scqr.fit$beta
  activeSet = getSet(beta.lasso, m)
  uniSet = activeSet$union
  voteSet = activeSet$vote
  Xunion = X[, uniSet, drop = FALSE]
  Xvote = X[, voteSet, drop = FALSE]
  ## scqr on the union set
  if (length(uniSet) > 0) {
    beta.union = scqrGauss(Xunion, Y, censor, tauSeq)
    test = exam(trueSig, uniSet, beta.union, betaMat[-1, ])
    TPR1[i] = test$TPR
    TNR1[i] = test$TNR
    FDR1[i] = test$FDR
    error1[i] = test$error
  }
  ## scqr on the majority vote set
  if (length(voteSet) > 0) {
    beta.vote = scqrGauss(Xvote, Y, censor, tauSeq)
    test = exam(trueSig, voteSet, beta.vote, betaMat[-1, ])
    TPR2[i] = test$TPR
    TNR2[i] = test$TNR
    FDR2[i] = test$FDR
    error2[i] = test$error
  }
  start = Sys.time()
  beta.lasso = SqrLasso(X, censor, Y, lambda0, tauSeq, h)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  

  ## HDCQR-Lasso using rqPen
  Z = cbind(1, X)
  start = Sys.time()
  beta.cqr = quantproc(Y, Z, censor, tauSeq, lambda0)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  activeSet = getSet(beta.cqr, m)
  uniSet = activeSet$union
  voteSet = activeSet$vote


  setTxtProgressBar(pb, i / M)
}




setwd("~/Dropbox/Conquer/SCQR/Code/Simulation/highd/hetero/lasso/")
mtc.lasso = as.matrix(read.csv("mtc_scqr.csv")[, -1])
mtc.scad = as.matrix(cbind(read.csv("mtc_scad1.csv")[, 2:101], 
                           read.csv("mtc_scad2.csv")[, 102:201],
                           read.csv("mtc_scad3.csv")[, 202:301],
                           read.csv("mtc_scad4.csv")[, 302:401],
                           read.csv("mtc_scad5.csv")[, 402:501]))
mtc.mcp = as.matrix(cbind(read.csv("mtc_mcp1.csv")[, 2:101],
                          read.csv("mtc_mcp2.csv")[, 102:201],
                          read.csv("mtc_mcp3.csv")[, 202:301],
                          read.csv("mtc_mcp4.csv")[, 302:401],
                          read.csv("mtc_mcp5.csv")[, 402:501]))

ind1 = 1:13
ind2 = 14:26
ind3 = 27:39
ind4 = 40:52
ind5 = 53:65

### Dataframe construction
TPR = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind1], rowMeans(mtc.scad, na.rm = TRUE)[ind1], rowMeans(mtc.mcp, na.rm = TRUE)[ind1])
TNR = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind2], rowMeans(mtc.scad, na.rm = TRUE)[ind2], rowMeans(mtc.mcp, na.rm = TRUE)[ind2])
PPV = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind3], rowMeans(mtc.scad, na.rm = TRUE)[ind3], rowMeans(mtc.mcp, na.rm = TRUE)[ind3])
FDR = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind4], rowMeans(mtc.scad, na.rm = TRUE)[ind4], rowMeans(mtc.mcp, na.rm = TRUE)[ind4])
error = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind5], rowMeans(mtc.scad, na.rm = TRUE)[ind5], rowMeans(mtc.mcp, na.rm = TRUE)[ind5])
dat = as.data.frame(cbind(TPR, TNR, PPV, FDR, error))
colnames(dat) = c("TPR", "TNR", "PPV", "FDR", "error")
dat$tau = rep(tauSeq, 3)
dat$type = c(rep("\\texttt{Lasso}", nTau), rep("\\texttt{SCAD}", nTau), rep("\\texttt{MCP}", nTau))
dat$type = factor(dat$type, levels = c("\\texttt{Lasso}", "\\texttt{SCAD}", "\\texttt{MCP}"))

### TPR
setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = tau, y = TPR)) +
  geom_line(aes(y = TPR, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Quantile level $\\tau$") + ylab("True positive rate") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.75), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

##TNR
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = tau, y = TNR)) +
  geom_line(aes(y = TNR, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Quantile level $\\tau$") + ylab("True negative rate") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.75), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

## PPV
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = tau, y = PPV)) +
  geom_line(aes(y = PPV, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Quantile level $\\tau$") + ylab("Precision") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.75), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


## FDR
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = tau, y = FDR)) +
  geom_line(aes(y = FDR, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Quantile level $\\tau$") + ylab("False discover rate") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.75), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


## error
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = tau, y = error)) +
  geom_line(aes(y = error, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Quantile level $\\tau$") + ylab("$\\ell_2$ error") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.75), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

## RE
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = tau, y = RE)) +
  geom_line(aes(y = RE, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Quantile level $\\tau$") + ylab("Relative error") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.75), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### Box plots
iii = 3
rst1 = c(mtc.lasso[ind4, ][iii, ], mtc.scad[ind4, ][iii, ], mtc.mcp[ind4, ][iii, ])
meth = c(rep("\\texttt{Lasso}", M), rep("\\texttt{SCAD}", M), rep("\\texttt{MCP}", M))
meth = factor(meth, levels = c("\\texttt{Lasso}", "\\texttt{SCAD}", "\\texttt{MCP}"))
dat = data.frame("est" = rst1, "method" = meth)

setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = method, y = est, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Estimation") + 
  #scale_y_continuous(breaks = c(5, 15, 25)) + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)





