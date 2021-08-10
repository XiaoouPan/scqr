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
library(glmnet)

rm(list = ls())
Rcpp::sourceCpp("src/hdscqr.cpp")

exam = function(beta, beta.hat) {
  m = ncol(beta)
  TPR = TNR = PPV = FDR = err = rep(0, m)
  for (i in 1:m) {
    TPR[i] = sum(beta[-1, i] != 0 & beta.hat[-1, i] != 0) / sum(beta[-1, i] != 0)
    TNR[i] = sum(beta[-1, i] == 0 & beta.hat[-1, i] == 0) / sum(beta[-1, i] == 0)
    PPV[i] = 0
    FDR[i] = 0
    if (sum(beta.hat[-1, i] != 0) > 0) {
      PPV[i] = sum(beta[-1, i] != 0 & beta.hat[-1, i] != 0) / sum(beta.hat[-1, i] != 0)
      FDR[i] = sum(beta[-1, i] == 0 & beta.hat[-1, i] != 0) / sum(beta.hat[-1, i] != 0)
    }
    err[i] = norm(beta[, i] - beta.hat[, i], "2")
  }
  return (list("TPR" = TPR, "TNR" = TNR, "PPV" = PPV, "FDR" = FDR, "error" = err))
}

calRes = function(Z, censor, Y, beta.hat, tauSeq, HSeq) {
  m = length(tauSeq)
  n = length(Y)
  rst = 0
  indi = 1 * (Y >= Z %*% beta.hat[, 1])
  accu = rep(tauSeq[1], n)
  res = censor * (1 - indi) - tauSeq[1]
  dev = sqrt(-2 * (res + censor * log(censor - res)))
  rst = rst + mean(dev)
  for (i in 2:m) {
    Hgap = HSeq[i] - HSeq[i - 1]
    accu = accu + indi * Hgap
    indi = 1 * (Y >= Z %*% beta.hat[, i])
    res = censor * (1 - indi) - accu
    dev = sqrt(-2 * (res + censor * log(censor - res)))
    rst = rst + mean(dev)
  }
  return (rst)
}

H = function(x) {
  return (-log(1 - x))
}

## Code from Zheng Peng and He (2018)
quantproc = function(y, x, delta, JJ, lambda, tol=1e-6){
  pp=dim(x)[2];                                                        #### the number of covariates
  tmpbeta=matrix(0,length(JJ),pp);                                     #### the coefficient matrix
  
  augy1=y[which(delta==1)];                                            #### the observed event time   
  augx1=x[which(delta==1),];                                           #### the corresponding covariates
  augx2= -colSums(augx1)                                               #### the 2nd part in the objective function                        
  augy=c(augy1, 5000, 5000);
  
  sto_weights=matrix(0,length(JJ),dim(x)[1]);                          #### stochastic weights
  rproc=matrix(0,length(JJ),dim(x)[1]);                                #### the indictor (y>=x%*%betahat);
  sto_weights[1,]=JJ[1]*2;                                             #### the initial weight 2tau0
  augx3=sto_weights[1,]%*%x;                                           #### the 3rd part in the objective function
  augx=rbind(augx1,augx2,augx3);
  tmpb = coef(rq(augy~0+augx,method="lasso",tau=0.5,lambda=lambda))
  #tmpb = coef(rq(y~0+x,method="lasso",tau=JJ[1],lambda=lambda))
  beta0 = tmpb*(abs(tmpb)>=tol)
  tmpbeta[1,]=beta0;
  rproc[1,]=1*(y>=x%*%beta0);                                          #### the initial indicator y>=x%*%betahat0;
  
  for(s in 2:length(JJ)){
    tuning=lambda
    Hm = H(JJ[s])-H(JJ[s-1]);                                        #### H(tau[s])-H(tau[s-1]) 
    sto_weights[s,]=sto_weights[s-1,]+2*Hm*rproc[s-1,];            #### update the stochastic weight for tau[s]
    augx3=sto_weights[s,]%*%x;                                     #### the 3rd part in the objective function
    augx=rbind(augx1,augx2,augx3);
    tmpb=coef(rq(augy~0+augx,method="lasso",tau=0.5,lambda=tuning));  #### quantile fit at tau[s];
    tmpbeta[s,]=tmpb*(abs(tmpb)>=tol);                             #### hard threshholding;
    rproc[s,] = 1*(y>x%*%tmpbeta[s,]);                               #### update the indicator y>=x%*%betahat at tau[s];	
  }
  return (t(tmpbeta));
}

cvCqr = function(X, censor, Y, lambdaSeq, tauSeq, K, folds) {
  l = length(lambdaSeq)
  mse = rep(0, l)
  Z = cbind(1, X)
  HSeq = H(tauSeq)
  for (k in 1:K) {
    indTest = which(folds == k)
    indTrain = which(folds != k)
    Ztrain = Z[indTrain, ]
    Ztest = Z[indTest, ]
    Ytrain = Y[indTrain]
    Ytest = Y[indTest]
    censorTrain = censor[indTrain]
    censorTest = censor[indTest]
    for (i in 1:l) {
      lambda = lambdaSeq[i]
      betaHat = quantproc(Ytrain, Ztrain, censorTrain, tauSeq, lambda)
      mse[i] = mse[i] + calRes(Ztest, censorTest, Ytest, betaHat, tauSeq, HSeq)
    }
  }
  index = which.min(mse)
  return (quantproc(Y, Z, censor, tauSeq, lambdaSeq[index]))
}


#### Quantile process with fixed scale, hard to visualize
n = 100
p = 200
s = 3
M = 1
kfolds = 3
h = (log(p) / n)^(1/4)
tauSeq = seq(0.1, 0.7, by = 0.05)
m = length(tauSeq)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
Sigma = toeplitz(0.5^(0:(p - 1)))
lambdaSeq = exp(seq(log(0.01), log(0.2), length.out = 50))

time = prop = rep(0, M)
TPR = TNR = PPV = FDR = error = matrix(0, m, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  #X = mvrnorm(n, rep(0, p), Sigma)
  X = matrix(runif(n * p, -1, 1), n, p)
  #err = rt(n, 2)
  err = rnorm(n)
  ## Homo
  #beta = c(runif(s, 1, 1.5), rep(0, p - s))
  beta = c(0.5, 1, 1.5, rep(0, p - s))
  betaMat = rbind(beta0, matrix(beta, p, nTau))
  logT = X %*% beta + err
  ## Hetero
  #X[, 1] = abs(X[, 1])
  #beta = c(runif(s - 1, 1, 1.5), rep(0, p - s))
  #betaMat = rbind(rep(0, nTau), beta0, matrix(beta, p - 1, nTau))
  #logT = X[, 1] * err + X[, -1] %*% beta
  
  #w = sample(1:3, n, prob = c(1/3, 1/3, 1/3), replace = TRUE)
  #logC = (w == 1) * rnorm(n, 0, 4) + (w == 2) * rnorm(n, 5, 1) + (w == 3) * rnorm(n, 10, 0.5)
  logC = rnorm(n, 3, sqrt(17.25))
  censor = logT <= logC
  prop[i] = 1 - sum(censor) / n
  Y = pmin(logT, logC)
  response = Surv(Y, censor, type = "right")
  folds = createFolds(censor, kfolds, FALSE)

  ## HDCQR-Lasso using quantreg
  beta.lasso = LASSO.fit(Y, X, tau = 0.1, lambda = lambdaSeq[20], intercept = TRUE, coef.cutoff = 0.0001)
  fit = cv.rq.pen(X, Y, tau =.1, lambda = lambdaSeq, penalty = "LASSO", nfolds = 3)
  
  
  
  start = Sys.time()
  beta.cqr = cvCqr(X, censor, Y, lambdaSeq, tauSeq, kfolds, folds)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  test = exam(betaMat, beta.cqr)
  TPR[, i] = test$TPR
  TNR[, i] = test$TNR
  PPV[, i] = test$PPV
  FDR[, i] = test$FDR
  error[, i] = test$error
  RE[, i] = test$RE
  
  ## SCQR-Lasso
  start = Sys.time()
  beta.lasso = cvSqrLasso(X, censor, Y, lambdaSeq, folds, tauSeq, kfolds, h)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  test = exam(betaMat, beta.lasso, beta.oracle)
  TPR[, i] = test$TPR
  TNR[, i] = test$TNR
  PPV[, i] = test$PPV
  FDR[, i] = test$FDR
  error[, i] = test$error
  RE[, i] = test$RE

  setTxtProgressBar(pb, i / M)
}




setwd("~/Dropbox/Conquer/SCQR/Code/Simulation/highd/homo")
mtc.lasso = as.matrix(read.csv("mtc_lasso.csv")[, -1])
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

ind1 = 1:11
ind2 = 12:22
ind3 = 23:33
ind4 = 34:44
ind5 = 45:55
ind6 = 56:66

### Dataframe construction
TPR = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind1], rowMeans(mtc.scad, na.rm = TRUE)[ind1], rowMeans(mtc.mcp, na.rm = TRUE)[ind1])
TNR = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind2], rowMeans(mtc.scad, na.rm = TRUE)[ind2], rowMeans(mtc.mcp, na.rm = TRUE)[ind2])
PPV = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind3], rowMeans(mtc.scad, na.rm = TRUE)[ind3], rowMeans(mtc.mcp, na.rm = TRUE)[ind3])
FDR = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind4], rowMeans(mtc.scad, na.rm = TRUE)[ind4], rowMeans(mtc.mcp, na.rm = TRUE)[ind4])
error = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind5], rowMeans(mtc.scad, na.rm = TRUE)[ind5], rowMeans(mtc.mcp, na.rm = TRUE)[ind5])
RE = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind6], rowMeans(mtc.scad, na.rm = TRUE)[ind6], rowMeans(mtc.mcp, na.rm = TRUE)[ind6])
dat = as.data.frame(cbind(TPR, TNR, PPV, FDR, error, RE))
colnames(dat) = c("TPR", "TNR", "PPV", "FDR", "error", "RE")
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





