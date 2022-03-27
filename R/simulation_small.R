###### Simulation code for scqr with small sample size

rm(list = ls())
Rcpp::sourceCpp("src/scqr.cpp")

library(quantreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(survival)
library(tikzDevice)
library(ggplot2)


n = 400
p = 5
M = 500
df = 1.5
tauSeq = seq(0.05, 0.8, by = 0.05)
grid = seq(0.05, 0.85, by = 0.05)
nTau = length(tauSeq)
beta0 = qt(tauSeq, df)
coef1 = eff1 = matrix(0, M, nTau)
coef2 = eff2 = matrix(0, M, nTau)
time = matrix(0, 2, M)
prop = rep(0, M)


pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  Sigma = toeplitz(0.5^(0:(p - 1)))
  X = mvrnorm(n, rep(0, p), Sigma)
  err = rt(n, df)
  ## Homo
  beta = runif(p, -2, 2)
  betaMat = rbind(beta0, matrix(beta, p, nTau))
  logT = X %*% beta + err
  ## Hetero
  #X[, 1] = abs(X[, 1])
  #beta = runif(p - 1, -2, 2)
  #betaMat = rbind(rep(0, nTau), beta0, matrix(beta, p - 1, nTau))
  #logT = X[, 1] * err + X[, -1] %*% beta
  logC = rnorm(n, 2, 4) 
  censor = logT <= logC
  prop[i] = 1 - sum(censor) / n
  Y = pmin(logT, logC)
  response = Surv(Y, censor, type = "right")
  
  ## Smoothed CQR
  start = Sys.time()
  list = scqrGauss(X, Y, censor, tauSeq)
  end = Sys.time()
  time[1, i] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, ] = sqrt(colSums((list$coeff - betaMat)^2))
  eff1[i, ] = list$coeff[1, ]
  #eff1[i, ] = list$coeff[2, ]
  
  ## Peng and Huang
  start = Sys.time()
  list = crq(response ~ X, method = "PengHuang", grid = grid)
  end = Sys.time()
  time[2, i] = as.numeric(difftime(end, start, units = "secs"))
  tt = ncol(list$sol)
  coef2[i, 1:tt] = sqrt(colSums((list$sol[2:(p + 2), ] - betaMat[, 1:tt])^2))
  eff2[i, 1:tt] = list$sol[2, ]
  #eff2[i, 1:tt] = list$sol[3, ]
  
  setTxtProgressBar(pb, i / M)
}



### Estimation plots
mean1 = colMeans(coef1)
mean2 = colMeans(coef2)
dat = rbind(cbind(tauSeq, mean1), cbind(tauSeq, mean2))
dat = as.data.frame(dat)
colnames(dat) = c("quantile", "coef")
dat$type = c(rep("\\texttt{Our method}", nTau), rep("\\texttt{Peng} \\& \\texttt{Huang}", nTau))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Our method}"))

setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = quantile, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid")) +
  theme_bw() + xlab("Quantile level $\\tau$") + ylab("Estimation error in $||\\cdot||_2$") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.8), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(0.9, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



### Quantile effects plots
mean1 = colMeans(eff1)
mean2 = colMeans(eff2)
dat = rbind(cbind(tauSeq, mean1), cbind(tauSeq, mean2))
dat = rbind(dat, cbind(tauSeq, beta0))
dat = as.data.frame(dat)
colnames(dat) = c("quantile", "eff")
dat$type = c(rep("\\texttt{Our method}", nTau), rep("\\texttt{Peng} \\& \\texttt{Huang}", nTau), rep("\\texttt{True effects}", nTau))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Our method}", "\\texttt{True effects}"))

setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = quantile, y = eff, color = type)) +
  geom_line(aes(y = eff, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  theme_bw() + xlab("Quantile level $\\tau$") + 
  ylab("Estimated quantile effects") + 
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.75, 0.28), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(0.9, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)





### Running time plots
meth = c(rep("Our method", 500), rep("Peng \\& Huang", 500))
meth = factor(meth, levels = c("Our method", "Peng \\& Huang"))
rst = data.frame("time" = c(time[1, ], time[2, ]), "method" = meth)

setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(rst, aes(x = method, y = time, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Elapsed time (in seconds)") + theme_bw() + 
  #scale_y_continuous(breaks = seq(0, 125, 25)) + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

