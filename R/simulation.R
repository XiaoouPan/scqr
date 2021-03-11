###### Simulation code for scqr

rm(list = ls())
Rcpp::sourceCpp("src/scqr.cpp")

library(quantreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(survival)
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
  nTau = length(tauSeq)
  
}


#### Quantile process with fixed scale, hard to visualize
n = 5000
p = n / 50
M = 100
tauSeq = seq(0.05, 0.95, by = 0.05)
grid = seq(0.05, 0.95, by = 0.05)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
coef1 = itcp1 = matrix(0, M, nTau)
coef2 = itcp2 = matrix(0, M, nTau)
time = matrix(0, 2, M)
prop = rep(0, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  Sigma = getSigma(p)
  #X = sqrt(12) * draw.d.variate.uniform(n, p, Sigma) - sqrt(3)
  X = mvrnorm(n, rep(0, p), Sigma)
  beta = runif(p, -2, 2)
  betaMat = rbind(beta0, matrix(beta, p, nTau))
  err = rt(n, 2)
  logT = X %*% beta + err
  w = sample(1:3, n, prob = c(1/3, 1/3, 1/3), replace = TRUE)
  logC = (w == 1) * rnorm(n, 0, 4) + (w == 2) * rnorm(n, 5, 1) + (w == 3) * rnorm(n, 10, 0.5)
  censor = logT <= logC
  prop[i] = 1 - sum(censor) / n
  Y = pmin(logT, logC)
  response = Surv(Y, censor, type = "right")
  
  start = Sys.time()
  list = scqrGauss(X, Y, censor, tauSeq)
  end = Sys.time()
  time[1, i] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, ] = sqrt(colSums((list$coeff - betaMat)^2))
  itcp1[i, ] = list$coeff[1, ]
  
  start = Sys.time()
  list = crq(response ~ X, method = "PengHuang", grid = grid)
  end = Sys.time()
  time[2, i] = as.numeric(difftime(end, start, units = "secs"))
  tt = ncol(list$sol)
  coef2[i, 1:tt] = sqrt(colSums((list$sol[2:(p + 2), ] - betaMat[, 1:tt])^2))
  itcp2[i, 1:tt] = list$sol[2, ]
  
  setTxtProgressBar(pb, i / M)
}

time
sum(coef2[, nTau - 1] == 0)
est = rbind(colMeans(coef1), colMeans(coef2))
plot(tauSeq[-nTau], est[2, -nTau], type = "b", pch = 1, lwd = 5, cex = 1, col = "red", axes = FALSE, xlim = c(0, 1), 
     ylim = c(min(pretty(range(est))), max(pretty(range(est)))), xlab = "", ylab = "")
lines(tauSeq, est[1, ], type = "b", pch = 2, lwd = 5, cex = 1, col = "blue")
est = rbind(colMeans(itcp1), colMeans(itcp2))
plot(tauSeq[-nTau], est[2, -nTau], type = "b", pch = 1, lwd = 5, cex = 1, col = "red", axes = FALSE, xlim = c(0, 1), 
     ylim = c(min(pretty(range(est))), max(pretty(range(est)))), xlab = "", ylab = "")
lines(tauSeq, est[1, ], type = "b", pch = 2, lwd = 5, cex = 1, col = "blue")
lines(tauSeq, beta0, type = "l", lwd = 5, cex = 1)


#### Growing dimension and sample size
nseq = seq(2000, 6000, by = 2000)
pseq = floor(nseq / 40)
l = length(nseq)
M = 1
coef = time = matrix(0, 2, l)
prop = rep(0, l)
tau = 0.7
tauSeq = seq(0.1, tau, by = 0.1)
grid = seq(0.1, tau + 0.1, by = 0.1)
beta0 = 1

pb = txtProgressBar(style = 3)
for (j in 1:l) {
  n = nseq[j]
  p = pseq[j]
  beta = runif(p, 1, 2)
  for (i in 1:M) {
    set.seed((j - 1) * M + i)
    X = matrix(rnorm(n * p), n, p)
    err = rt(n, 2) - qt(tau, 2)
    logT = beta0 + X %*% beta + err
    logC = rnorm(n, 5, 4)
    censor = logT <= logC
    prop[j] = prop[j] + 1 - sum(censor) / n
    Y = pmin(logT, logC)
    response = Surv(Y, censor, type = "right")
    
    start = Sys.time()
    list = scqrGauss(X, Y, censor, tauSeq)
    end = Sys.time()
    time[1, j] = time[1, j] + as.numeric(difftime(end, start, units = "secs"))
    coef[1, j] = coef[1, j] + mean((list$coeff[-1] - beta)^2)
    
    start = Sys.time()
    list = crq(response ~ X, method = "PengHuang", grid = grid)
    end = Sys.time()
    time[2, j] = time[2, j] + as.numeric(difftime(end, start, units = "secs"))
    coef[2, j] = coef[2, j] + mean((as.numeric(list$sol[3:(p + 2), length(tauSeq)]) - beta)^2)
    
    setTxtProgressBar(pb, ((j - 1) * M + i) / (l * M))
  }
}

prop = prop / M
time = time / M
coef = coef / M
#write.csv(prop, "~/Dropbox/Conquer/censoredQR/Code/Simulation/prop.csv")
#write.csv(time, "~/Dropbox/Conquer/censoredQR/Code/Simulation/time.csv")
#write.csv(coef, "~/Dropbox/Conquer/censoredQR/Code/Simulation/coef.csv")
#time = as.matrix(read.csv("~/Dropbox/Conquer/censoredQR/Code/Simulation/time.csv"))[, -1]
#coef = as.matrix(read.csv("~/Dropbox/Conquer/censoredQR/Code/Simulation/coef.csv"))[, -1]

setwd("~/Dropbox/Conquer/censoredQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 6, height = 6)
plot(coef[2, ], type = "b", pch = 1, lwd = 5, cex = 1, col = "red", axes = FALSE, ylim = c(min(pretty(range(coef))), max(pretty(range(coef)))), xlab = "", ylab = "")
lines(coef[1, ], type = "b", pch = 2, lwd = 5, cex = 1, col = "blue")
grid(col = "grey", lty = 2)
color = c("red", "blue")
labels = c("\\texttt{censored QR}", "\\texttt{proposed method}")
pch = c(1, 2)
legend(7, 0.014, labels, col = color, pch = pch, lwd = 5, cex = 2, box.lwd = 1, bg = "white")
axis(1, c(1, seq(5, 20, by = 5)), nseq[c(1, seq(5, 20, by = 5))] / 1000, tick = TRUE, line = -1, cex.axis = 1.8)
axis(2, pretty(range(coef)), line = -0.7, cex.axis = 1.8)
title(xlab = "Sample Size (in thousands)", line = 1.5, cex.lab = 2)
title(ylab = "Estimation Error", line = 1.7, cex.lab = 2)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

## Time plot
tikz("plot.tex", standAlone = TRUE, width = 6, height = 6)
plot(time[2, ], type = "b", pch = 1, cex = 1, lwd = 5, col = "red", axes = FALSE, ylim = c(min(pretty(range(time))), max(pretty(range(time)))), xlab = "", ylab = "")
lines(time[1, ], type = "b", pch = 2, cex = 1, lwd = 5, col = "blue")
grid(col = "grey", lty = 2)
axis(1, c(1, seq(5, 20, by = 5)), nseq[c(1, seq(5, 20, by = 5))] / 1000, tick = TRUE, line = -0.6, cex.axis = 1.8)
axis(2, pretty(range(time)), line = -0.7, cex.axis = 1.8)
title(xlab = "Sample Size (in thousands)", line = 2, cex.lab = 2)
title(ylab = "Time Elapsed", line = 1.7, cex.lab = 2)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

