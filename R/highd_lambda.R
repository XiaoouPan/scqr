###### Simulation code for scqr
library(quantreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(survival)
library(tikzDevice)
library(ggplot2)
library(glmnet)
library(caret)

rm(list = ls())
Rcpp::sourceCpp("src/hdscqr.cpp")


metric = function(beta, beta.hat) {
  m = ncol(beta)
  TPR = TNR = PPV = FDR = error = rep(0, m)
  for (i in 1:m) {
    TPR[i] = sum(beta[-1, i] != 0 & beta.hat[-1, i] != 0) / sum(beta[-1, i] != 0)
    TNR[i] = sum(beta[-1, i] == 0 & beta.hat[-1, i] == 0) / sum(beta[-1, i] == 0)
    PPV[i] = 0
    FDR[i] = 0
    if (sum(beta.hat[-1, i] != 0) > 0) {
      PPV[i] = sum(beta[-1, i] != 0 & beta.hat[-1, i] != 0) / sum(beta.hat[-1, i] != 0)
      FDR[i] = sum(beta[-1, i] == 0 & beta.hat[-1, i] != 0) / sum(beta.hat[-1, i] != 0)
    }
    error[i] = norm(beta[, i] - beta.hat[, i], "2")
  }
  return (list("TPR" = TPR, "TNR" = TNR, "PPV" = PPV, "FDR" = FDR, "error" = error))
}


exam = function(beta, beta.hat, HSeq) {
  m = ncol(beta)
  err = 0
  for (i in 1:(m - 1)) {
    err = err + norm(beta[, i] - beta.hat[, i], "2") * (HSeq[i + 1] - HSeq[i])
  }
  return (err)
}

calRes = function(X, censor, Y, beta.hat, tauSeq, HSeq) {
  m = ncol(beta.hat)
  res = censor * (Y <= (beta.hat[1, m] + X %*% beta.hat[-1, m])) - tauSeq[1]
  for (i in 1:(m - 1)) {
    res = res - (Y >= (beta.hat[1, i] + X %*% beta.hat[-1, i])) * (HSeq[i + 1] - HSeq[i])
  }
  return (mean(res^2))
}

calResSum = function(X, censor, Y, beta.hat, tauSeq) {
  m = ncol(beta.hat)
  res = 0
  X.cen = X[which(censor == 1), ]
  Y.cen = Y[which(censor == 1)]
  for (i in 1:m) {
    res.cen = Y.cen - beta.hat[1, i] - X.cen %*% beta.hat[-1, i]
    temp = res.cen * (tauSeq[i] - (res.cen < 0))
    res = res + mean(temp)
  }
  return (res / m)
}


#### Quantile process with fixed scale, hard to visualize
n = 400
p = 200
s = 10
M = 1
kfolds = 3
h = 0.2
tauSeq = seq(0.2, 0.7, by = 0.05)
m = length(tauSeq)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
lambdaSeq = exp(seq(log(0.01), log(0.2), length.out = 50))
HSeq = as.numeric(getH(tauSeq))
error1 = error2 = error3 = res1 = res2 = res3 = matrix(0, 50, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  Sigma = toeplitz(0.5^(0:(p - 1)))
  X = mvrnorm(n, rep(0, p), Sigma)
  err = rt(n, 2)
  ## Homo
  beta = c(runif(s, 1.5, 2), rep(0, p - s))
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
  Y = pmin(logT, logC)
  
  for (j in 1:50) {
    ## SCQR-Lasso
    beta.lasso = SqrLasso(X, censor, Y, lambdaSeq[j], tauSeq, h)
    error1[j, i] = exam(betaMat, beta.lasso, HSeq)
    res1[j, i] = calResSum(X, censor, Y, beta.lasso, tauSeq)

    ## SCQR-SCAD
    beta.scad = SqrScad(X, censor, Y, lambdaSeq[j], tauSeq, h)
    error2[j, i] = exam(betaMat, beta.scad, HSeq)
    res2[j, i] = calResSum(X, censor, Y, beta.scad, tauSeq)
    
    ## SCQR-MCP
    beta.mcp = SqrMcp(X, censor, Y, lambdaSeq[j], tauSeq, h)
    error3[j, i] = exam(betaMat, beta.mcp, HSeq)
    res3[j, i] = calResSum(X, censor, Y, beta.mcp, tauSeq)
    
    setTxtProgressBar(pb, (j + (i - 1) * 50) / (50 * M))
  }
}


cbind(rowMeans(error1), rowMeans(error2), rowMeans(error3))
cbind(rowMeans(res1), rowMeans(res2), rowMeans(res3))
plot(rowMeans(error1), type = "l")
lines(rowMeans(error2), type = "l", col = "red")
lines(rowMeans(error3), type = "l", col = "blue")
plot(rowMeans(res1), type = "l")
lines(rowMeans(res2), type = "l", col = "red")
lines(rowMeans(res3), type = "l", col = "blue")


setwd("~/Dropbox/Conquer/SCQR/Code/Simulation/highd/homo")
mtc.lasso = as.matrix(read.csv("lam_lasso.csv")[, -1])
mtc.scad = as.matrix(read.csv("lam_scad.csv")[, -1])
mtc.mcp = as.matrix(read.csv("lam_mcp.csv")[, -1])

lambdaSeq = exp(seq(log(0.02), log(0.2), length.out = 50))
ind1 = 1:50
ind2 = 51:100

### Dataframe construction
error = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind1], rowMeans(mtc.scad, na.rm = TRUE)[ind1], rowMeans(mtc.mcp, na.rm = TRUE)[ind1])
res = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind2], rowMeans(mtc.scad, na.rm = TRUE)[ind2], rowMeans(mtc.mcp, na.rm = TRUE)[ind2])
dat = as.data.frame(cbind(error, res))
colnames(dat) = c("error", "res")
dat$lambda = rep(lambdaSeq, 3)
dat$type = c(rep("\\texttt{Lasso}", 50), rep("\\texttt{SCAD}", 50), rep("\\texttt{MCP}", 50))
dat$type = factor(dat$type, levels = c("\\texttt{Lasso}", "\\texttt{SCAD}", "\\texttt{MCP}"))


## error
setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = lambda, y = error)) +
  geom_line(aes(y = error, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("$\\lambda$") + ylab("$Global error") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.75), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

## residual
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = lambda, y = res)) +
  geom_line(aes(y = res, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("$\\lambda$") + ylab("Martingale residual") +
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





