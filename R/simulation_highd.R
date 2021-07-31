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


getSigma = function(p) {
  sig = diag(p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      sig[i, j] = sig[j, i] = 0.5^(j - i)
    }
  }
  return (sig)
}

exam = function(beta, beta.hat, beta.oracle) {
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
  tt = ncol(beta.oracle)
  RE = rep(0, tt)
  for (i in 1:tt) {
    err.ora = norm(beta[, i] - beta.oracle[, i], "2")
    RE[i] = err[i] / err.ora
  }
  return (list("TPR" = TPR, "TNR" = TNR, "PPV" = PPV, "FDR" = FDR, "error" = err, "RE" = RE))
}


#### Quantile process with fixed scale, hard to visualize
n = 400
p = 1000
s = 10
M = 500
kfolds = 3
tauSeq = seq(0.2, 0.7, by = 0.05)
m = length(tauSeq)
grid = seq(0.2, 0.75, by = 0.05)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
Sigma = getSigma(p)
lambdaSeq = exp(seq(log(0.02), log(0.3), length.out = 50))

time = prop = rep(0, M)
TPR = TNR = PPV = FDR = error = RE = matrix(0, m, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  #X = sqrt(12) * draw.d.variate.uniform(n, p, Sigma) - sqrt(3)
  X = mvrnorm(n, rep(0, p), Sigma)
  #Sigma = getSigma(45)
  #X = cbind(mvrnorm(n, rep(0, 45), Sigma), 4 * draw.d.variate.uniform(n, 45, Sigma) - 2, matrix(rbinom(10 * n, 1, c(0.5, 0.5)), n, 10))
  err = rt(n, 2)
  ## Homo
  beta = c(runif(s, 2, 3), rep(0, p - s))
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
  folds = createFolds(censor, kfolds, FALSE)
  ##  Check if there are enough test samples
  if (sum(censor[folds == 1] == 1) <= 5 | sum(censor[folds == 2] == 1) <= 5 | sum(censor[folds == 3] == 1) <= 5) {
    setTxtProgressBar(pb, i / M)
    next
  }
  
  ## Peng and Huang on the oracle set
  list = crq(response ~ X[, 1:s], method = "PengHuang", grid = grid)
  beta.oracle = list$sol[2:(s + 2), ]
  tt = ncol(beta.oracle)
  beta.oracle = rbind(beta.oracle, matrix(0, p - s, tt))
  
  fit = cv.glmnet(X, Y, nlambda = 50)
  s.hat = sum(as.numeric(coef(fit, s = fit$lambda.min)) != 0)
  h = max(min((s.hat * sqrt(log(p) / n) + (s.hat * log(p) / n)^(0.25)) / 2, 1), 0.1)
  
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
  
  ## SCQR-SCAD
  start = Sys.time()
  beta.scad = cvSqrScad(X, censor, Y, lambdaSeq, folds, tauSeq, kfolds, h)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  test = exam(betaMat, beta.scad, beta.oracle)
  TPR[, i] = test$TPR
  TNR[, i] = test$TNR
  PPV[, i] = test$PPV
  FDR[, i] = test$FDR
  error[, i] = test$error
  RE[, i] = test$RE
  
  ## SCQR-MCP
  start = Sys.time()
  beta.mcp = cvSqrMcp(X, censor, Y, lambdaSeq, folds, tauSeq, kfolds, h)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  test = exam(betaMat, beta.mcp, beta.oracle)
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
mtc.scad = as.matrix(read.csv("mtc_scad.csv")[, -1])
mtc.mcp = as.matrix(read.csv("mtc_mcp.csv")[, -1])

ind1 = 1:13
ind2 = 14:26
ind3 = 27:39
ind4 = 40:52
ind5 = 53:65

### Dataframe construction
TPR = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind1], rowMeans(mtc.scad, na.rm = TRUE)[ind1], rowMeans(mtc.mcp, na.rm = TRUE)[ind1])
TNR = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind2], rowMeans(mtc.scad, na.rm = TRUE)[ind2], rowMeans(mtc.mcp, na.rm = TRUE)[ind2])
PPV = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind3], rowMeans(mtc.scad, na.rm = TRUE)[ind3], rowMeans(mtc.mcp, na.rm = TRUE)[ind3])
error = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind4], rowMeans(mtc.scad, na.rm = TRUE)[ind4], rowMeans(mtc.mcp, na.rm = TRUE)[ind4])
RE = c(rowMeans(mtc.lasso, na.rm = TRUE)[ind5], rowMeans(mtc.scad, na.rm = TRUE)[ind5], rowMeans(mtc.mcp, na.rm = TRUE)[ind5])
dat = as.data.frame(cbind(TPR, TNR, PPV, error, RE))
colnames(dat) = c("TPR", "TNR", "PPV", "error", "RE")
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

## error, not good
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

## RE, not good
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





