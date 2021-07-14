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

jack = function(X, Y, censor, n, grid, nTau, B = 1000) {
  mboot = 2 * ceiling(sqrt(n))
  rst = matrix(NA, p + 1, B)
  for (i in 1:B) {
    s = sample(1:n,  n - mboot)
    s = sort(s)
    yb = Y[s]
    xb = X[s, ]
    cb = censor[s]  
    resb = Surv(yb, cb, type = "right")
    list = crq(resb ~ xb, method = "PengHuang", grid = grid)
    tt = ncol(list$sol)
    if (tt >= nTau) {
      rst[, i] = list$sol[2:(p + 2), nTau]
    }
  }
  return (rst)
}

xyPair = function(X, Y, censor, n, grid, nTau, B = 1000) {
  rst = matrix(NA, p + 1, B)
  for (i in 1:B) {
    w = sample(1:n, n, replace = TRUE)
    s = sort(w)
    yb = Y[s]
    xb = X[s, ]
    cb = censor[s]  
    resb = Surv(yb, cb, type = "right")
    list = crq(resb ~ xb, method = "PengHuang", grid = grid)
    tt = ncol(list$sol)
    if (tt >= nTau) {
      rst[, i] = as.numeric(list$sol[2:(p + 2), nTau])
    }
  }
  return (rst)
}


getCoverage = function(ci, beta) {
  return (ci[-1, 1] <  beta & ci[-1, 2] > beta)
}

getWidth = function(ci) {
  return (ci[-1, 2] - ci[-1, 1])
}


getCoverPlot = function(cover1, cover2, cover3, p) {
  rst1 = colMeans(cover1, na.rm = TRUE)
  rst2 = colMeans(cover2, na.rm = TRUE)
  rst3 = colMeans(cover3, na.rm = TRUE)
  meth = c(rep("Multi", p), rep("Jack", p), rep("Pair", p))
  meth = factor(meth, levels = c("Multi", "Jack", "Pair"))
  rst = data.frame("cover" = c(rst1, rst2, rst3), "method" = meth)
  return (rst)
}

getWidthPlot = function(width1, width2, width3, j, M) {
  meth = c(rep("Multi", M), rep("Jack", M), rep("Pair", M))
  meth = factor(meth, levels = c("Multi", "Jack", "Pair"))
  rst = data.frame("width" = c(width1[, j], width2[, j], width3[, j]), "method" = meth)
  return (rst)
}

n = 1000
p = n / 50
M = 500
tauSeq = seq(0.1, 0.5, by = 0.1)
grid = seq(0.1, 0.6, by = 0.1)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
alpha = 0.05
z = qnorm(1 - alpha / 2)
cover1 = cover2 = cover3 = cover4 = cover5 = cover6 = cover7 = cover8 = cover9 = matrix(NA, M, p)
width1 = width2 = width3 = width4 = width5 = width6 = width7 = width8 = width9 = matrix(NA, M, p)
time1 = time2 = time3 = rep(NA, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  #Sigma = getSigma(p)
  #X = mvrnorm(n, rep(0, p), Sigma)
  Sigma = getSigma(8)
  X = cbind(mvrnorm(n, rep(0, 8), Sigma), 4 * draw.d.variate.uniform(n, 8, Sigma) - 2, matrix(rbinom(4 * n, 1, c(0.5, 0.5)), n, 4))
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
  Y = pmin(logT, logC)
  response = Surv(Y, censor, type = "right")
  
  ## Smoothed CQR
  start = Sys.time()
  list = scqrGaussInf(X, Y, censor, tauSeq, B = 1000)
  end = Sys.time()
  time1[i] = as.numeric(difftime(end, start, units = "secs"))
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
  
  ## Delete-d jackknife by Portnoy
  start = Sys.time()
  list = crq(response ~ X, method = "PengHuang", grid = grid)
  beta.boot = jack(X, Y, censor, n, grid, nTau, B = 1000)
  end = Sys.time()
  time2[i] = as.numeric(difftime(end, start, units = "secs"))
  beta.hat = list$sol
  tt = ncol(list$sol)
  if (tt >= nTau) {
    ci.list = getPivCI(beta.hat[2:(p + 2), nTau], beta.boot, alpha)
    ci.per = ci.list$perCI
    ci.piv = ci.list$pivCI
    ci.norm = getNormCI(beta.hat[2:(p + 2), nTau], rowSds(beta.boot), z)
    cover4[i, ] = getCoverage(ci.per, beta)
    cover5[i, ] = getCoverage(ci.piv, beta)
    cover6[i, ] = getCoverage(ci.norm, beta)
    width4[i, ] = getWidth(ci.per)
    width5[i, ] = getWidth(ci.piv)
    width6[i, ] = getWidth(ci.norm)
  }
  
  ## xy-pair resampling method
  start = Sys.time()
  list = crq(response ~ X, method = "PengHuang", grid = grid)
  beta.boot = xyPair(X, Y, censor, n, grid, nTau, B = 1000)
  end = Sys.time()
  time3[i] = as.numeric(difftime(end, start, units = "secs"))
  beta.hat = list$sol
  tt = ncol(list$sol)
  if (tt >= nTau) {
    ci.list = getPivCI(beta.hat[2:(p + 2), nTau], beta.boot, alpha)
    ci.per = ci.list$perCI
    ci.piv = ci.list$pivCI
    ci.norm = getNormCI(beta.hat[2:(p + 2), nTau], rowSds(beta.boot), z)
    cover7[i, ] = getCoverage(ci.per, beta)
    cover8[i, ] = getCoverage(ci.piv, beta)
    cover9[i, ] = getCoverage(ci.norm, beta)
    width7[i, ] = getWidth(ci.per)
    width8[i, ] = getWidth(ci.piv)
    width9[i, ] = getWidth(ci.norm)
  }
  
  setTxtProgressBar(pb, i / M)
}

rbind(time1, time2, time3)


#write.csv(time, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/time.csv")
#write.csv(cover1, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/cover1.csv")
#write.csv(cover2, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/cover2.csv")
#write.csv(cover3, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/cover3.csv")
#write.csv(width1, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/width1.csv")
#write.csv(width2, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/width2.csv")
#write.csv(width3, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Inference/width3.csv")




setwd("~/Dropbox/Conquer/SCQR/Code/Simulation/Inference/hetero")
time = as.matrix(read.csv("time.csv")[, -1])
cover_mb = as.matrix(read.csv("cover_mb.csv")[, -1])
cover_jack = as.matrix(read.csv("cover_jack.csv")[, -1])
cover_pair = as.matrix(read.csv("cover_pair.csv")[, -1])
width_mb = as.matrix(read.csv("width_mb.csv")[, -1])
width_jack = as.matrix(read.csv("width_jack.csv")[, -1])
width_pair = as.matrix(read.csv("width_pair.csv")[, -1])




### coverage plot
cover1 = cover_mb[1:500, ]
cover2 = cover_jack[1:500, ]
cover3 = cover_pair[1:500, ]

summ = getCoverPlot(cover1, cover2, cover3, p)
setwd("~/Dropbox/Conquer/censoredQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(summ, aes(x = method, y = cover, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Coverage") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



### width plot
width1 = width_mb[1:500, ]
width2 = width_jack[1:500, ]
width3 = width_pair[1:500, ]

summ = getWidthPlot(width1, width2, width3, 1, M)
setwd("~/Dropbox/Conquer/censoredQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(summ, aes(x = method, y = width, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("CI Width") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)




### time plot
meth = c(rep("Multi", M), rep("Jack", M), rep("Pair", M))
meth = factor(meth, levels = c("Multi", "Jack", "Pair"))
rst = data.frame("time" = c(time[1, ], time[2, ], time[3, ]), "method" = meth)
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(rst, aes(x = method, y = time, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Elapsed time (in seconds)") + theme_bw() + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



