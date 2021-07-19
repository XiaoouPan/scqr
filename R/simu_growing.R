###### Simulation code for scqr

rm(list = ls())
Rcpp::sourceCpp("src/scqr.cpp")

library(quantreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(survival)
library(tikzDevice)
library(ggplot2)

getSigma = function(p) {
  sig = diag(p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      sig[i, j] = sig[j, i] = 0.5^(j - i)
    }
  }
  return (sig)
}

accuError = function(betahat, beta, tauSeq) {
  nTau = min(length(tauSeq), ncol(betahat) + 1)
  diff = betahat - beta
  err = sqrt(colSums((betahat - beta)^2))
  accu = 0
  for (k in 2:nTau) {
    accu = accu + err[k - 1] * (tauSeq[k] - tauSeq[k - 1])
  }
  return (accu)
}


#### Growing dimension and sample size
nseq = seq(1000, 20000, by = 1000)
pseq = floor(nseq / 100)
l = length(nseq)
tauSeq = seq(0.05, 0.8, by = 0.05)
grid = seq(0.05, 0.85, by = 0.05)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
M = 500
coef1 = coef2 = coef3 = time1 = time2 = time3 = prop = matrix(0, M, l)
accu1 = accu2 = accu3 = matrix(0, M, l)

pb = txtProgressBar(style = 3)
for (j in 1:l) {
  n = nseq[j]
  p = pseq[j]
  for (i in 1:M) {
    set.seed((j - 1) * M + i)
    Sigma = getSigma(p)
    X = mvrnorm(n, rep(0, p), Sigma)
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
    prop[i, j] = 1 - sum(censor) / n
    Y = pmin(logT, logC)
    response = Surv(Y, censor, type = "right")
    
    ## Smoothed CQR
    start = Sys.time()
    list = scqrGauss(X, Y, censor, tauSeq)
    end = Sys.time()
    time1[i, j] = as.numeric(difftime(end, start, units = "secs"))
    accu1[i, j] = accuError(list$coeff, betaMat, tauSeq)
    coef1[i, j] = norm(list$coeff[, 14] - betaMat[, 14], "2")
    
    ## Peng and Huang
    start = Sys.time()
    list = crq(response ~ X, method = "PengHuang", grid = grid)
    end = Sys.time()
    time2[i, j] = as.numeric(difftime(end, start, units = "secs"))
    tt = ncol(list$sol)
    if (tt >= nTau - 1) {
      accu2[i, j] = accuError(list$sol[2:(p + 2), ], betaMat[, 1:tt], tauSeq)
    }
    coef2[i, j] = norm(list$sol[2:(p + 2), 14] - betaMat[, 14], "2")
    
    ## Portnoy
    start = Sys.time()
    list = crq(response ~ X, method = "Portnoy", grid = tauSeq)
    end = Sys.time()
    time3[i, j] = as.numeric(difftime(end, start, units = "secs"))
    accu3[i, j] = accuError(list$sol[2:(p + 2), 2:17], betaMat, tauSeq)
    coef3[i, j] = norm(list$sol[2:(p + 2), 15] - betaMat[, 14], "2")
    
    setTxtProgressBar(pb, ((j - 1) * M + i) / (l * M))
  }
}



#write.csv(prop, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Growing/prop.csv")
#write.csv(time1, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Growing/time1.csv")
#write.csv(time2, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Growing/time2.csv")
#write.csv(coef1, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Growing/coef1.csv")
#write.csv(coef2, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Growing/coef2.csv")
#time = as.matrix(read.csv("~/Dropbox/Conquer/censoredQR/Code/Simulation/time.csv"))[, -1]
#coef = as.matrix(read.csv("~/Dropbox/Conquer/censoredQR/Code/Simulation/coef.csv"))[, -1]


setwd("~/Dropbox/Conquer/SCQR/Code")
coef = as.matrix(read.csv("Simulation/growing/coef_hetero.csv")[, -1])
coef = cbind(coef, as.matrix(read.csv("Simulation/growing/coef_hetero_15.csv")[, -c(1, 2)]))
coef = cbind(coef, as.matrix(read.csv("Simulation/growing/coef_hetero_18.csv")[, -1]))
coef = cbind(coef, as.matrix(read.csv("Simulation/growing/coef_hetero_20.csv")[, -1]))
coef1 = coef[1:500, ]
coef2 = coef[501:1000, ]
accu = as.matrix(read.csv("Simulation/growing/accu_hetero.csv")[, -1])
accu = cbind(accu, as.matrix(read.csv("Simulation/growing/accu_hetero_15.csv")[, -c(1, 2)]))
accu = cbind(accu, as.matrix(read.csv("Simulation/growing/accu_hetero_18.csv")[, -1]))
accu = cbind(accu, as.matrix(read.csv("Simulation/growing/accu_hetero_20.csv")[, -1]))
accu1 = accu[1:500, ]
accu2 = accu[501:1000, ]
time = as.matrix(read.csv("Simulation/growing/time_hetero.csv")[, -1])
time = cbind(time, as.matrix(read.csv("Simulation/growing/time_hetero_15.csv")[, -c(1, 2)]))
time = cbind(time, as.matrix(read.csv("Simulation/growing/time_hetero_18.csv")[, -1]))
time = cbind(time, as.matrix(read.csv("Simulation/growing/time_hetero_20.csv")[, -1]))
time1 = time[1:500, ]
time2 = time[501:1000, ]


### Global estimation error 
mean1 = colMeans(accu1, na.rm = TRUE)
mean2 = colMeans(accu2, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Our method}", l), rep("\\texttt{Peng} \\& \\texttt{Huang}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Our method}"))

tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Global estimation error") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.7, 0.15), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### Estimation error at a certain quantile
mean1 = colMeans(coef1, na.rm = TRUE)
mean2 = colMeans(coef2, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Our method}", l), rep("\\texttt{Peng} \\& \\texttt{Huang}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Our method}"))

tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Estimation error at $\\tau = 0.7$") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  #theme(legend.position = c(0.7, 0.15), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
  #      legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
  #      axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



## Time plot
mean1 = colMeans(time1, na.rm = TRUE)
mean2 = colMeans(time2, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2))
dat = as.data.frame(dat)
colnames(dat) = c("size", "time")
dat$type = c(rep("\\texttt{Our method}", l), rep("\\texttt{Peng} \\& \\texttt{Huang}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Our method}"))

tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = time)) +
  geom_line(aes(y = time, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Elapsed time (in seconds)") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.35, 0.8), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



















### Old graphs without using ggplot
setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
plot(nseq, coef[2, ], type = "b", pch = 1, lwd = 5, cex = 1, col = "red", axes = FALSE, xlim = c(1000, 20000), 
     ylim = c(min(pretty(range(coef))), max(pretty(range(coef)))), xlab = "", ylab = "")
lines(nseq, coef[1, ], type = "b", pch = 2, lwd = 5, cex = 1, col = "blue")
axis(1, nseq[c(2, 4, 6, 8, 10)], nseq[c(2, 4, 6, 8, 10)] / 1000, tick = TRUE, line = 0, cex.axis = 1.5)
axis(2, pretty(range(coef)), line = 0, cex.axis = 1.5)
box()
abline(h = pretty(range(coef)), v = nseq[c(2, 4, 6, 8, 10)], col = "gray", lty = 2)
#color = c("red", "blue")
#labels = c("\\texttt{censored QR}", "\\texttt{proposed method}")
#pch = c(1, 2)
#legend(7, 0.014, labels, col = color, pch = pch, lwd = 5, cex = 2, box.lwd = 1, bg = "white")
title(xlab = "Sample size (in thousands)", line = 2.5, cex.lab = 1.8)
title(ylab = "Global estimation error", line = 2.5, cex.lab = 1.8)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

## Time plot
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
plot(nseq, time[2, ], type = "b", pch = 1, cex = 1, lwd = 5, col = "red", axes = FALSE, xlim = c(1000, 20000), 
     ylim = c(min(pretty(range(time))), max(pretty(range(time)))), xlab = "", ylab = "")
lines(nseq, time[1, ], type = "b", pch = 2, cex = 1, lwd = 5, col = "blue")
axis(1, nseq[c(2, 4, 6, 8, 10)], nseq[c(2, 4, 6, 8, 10)] / 1000, tick = TRUE, line = 0, cex.axis = 1.5)
axis(2, pretty(range(time)), line = 0, cex.axis = 1.5)
box()
abline(h = pretty(range(time)), v = nseq[c(2, 4, 6, 8, 10)], col = "gray", lty = 2)
color = c("red", "blue")
labels = c("CQR", "smoothed CQR")
pch = c(1, 2)
legend("topleft", labels, col = color, pch = pch, lwd = 3, cex = 1.5, box.lwd = 1, bg = "white")
title(xlab = "Sample size (in thousands)", line = 2.5, cex.lab = 1.8)
title(ylab = "Time elapsed (in seconds)", line = 2.5, cex.lab = 1.8)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

