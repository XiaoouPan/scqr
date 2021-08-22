##### real data section
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
Rcpp::sourceCpp("src/scqr.cpp")


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


#### Mayo data
X = cbind(pbc$age, pbc$edema, log(pbc$bili), log(pbc$albumin), log(pbc$protime)) 
idx = which(rowSums(is.na(X)) > 0)
X = X[-idx, ]
Y = log(pbc$time)[-idx]
censor = as.numeric(pbc$status == 2)[-idx]
n = nrow(X)
p = ncol(X)
1 - sum(censor) / n  ## censoring rate: 61.5%
tauSeq = seq(0.01, 0.9, by = 0.01)
grid = seq(0.01, 0.91, by = 0.01)
nTau = length(tauSeq)
alpha = 0.05
h = 0.5 * ((p + log(n)) / n)^(0.4)
lower = upper = NULL

list = scqrGauss(X, Y, censor, tauSeq, h = h)
est.h2 = list$coeff
for (i in 1:nTau) {
  list = scqrGaussInf(X, Y, censor, tauSeq[1:i], B = 1000)
  beta.boot = list$boot
  ci.list = getPivCI(est[, i], beta.boot, alpha)
  lower = cbind(lower, ci.list$perCI[, 1])
  upper = cbind(upper, ci.list$perCI[, 2])
}

response = Surv(Y, censor, type = "right")
list = crq(response ~ X, method = "PengHuang", grid = grid)
ph = list$sol[2:7, ]



### plots
setwd("~/Dropbox/Conquer/SCQR/Code")
est = as.matrix(read.csv("real/est.csv")[, -1])
low_up = as.matrix(read.csv("real/low_up.csv")[, -1])
beta.scqr = est[2:6, ]
beta.cqr = est[8:12, ]
lower = low_up[2:6, ]
upper = low_up[8:12, ]
beta.scqr2 = est.h2[2:6, ]

k = 2
dat = cbind(tauSeq, beta.scqr[k, ], lower[k, ], upper[k, ])
dat = rbind(dat, cbind(tauSeq, beta.scqr2[k, ], beta.scqr2[k, ], beta.scqr2[k, ]))
dat = rbind(dat, cbind(tauSeq, beta.cqr[k, ], beta.cqr[k, ], beta.cqr[k, ]))
dat = as.data.frame(dat)
colnames(dat) = c("quantile", "coeff", "lower", "upper")
dat$type = c(rep("\\texttt{Our method}", nTau), rep("\\texttt{Our method 2}", nTau), rep("\\texttt{Peng} \\& \\texttt{Huang}", nTau))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Our method}", "\\texttt{Our method 2}"))

tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = quantile, y = coeff, color = type)) +
  geom_line(aes(y = coeff, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid", "dashed")) +
  #geom_ribbon(aes(y = coeff, ymin = lower, ymax = upper, fill = type), alpha = 0.3) + 
  theme_bw() + xlab("Quantile level $\\tau$") + 
  ylab("Estimated coefficients") + 
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20)) 
  #theme(legend.position = c(0.65, 0.8), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
  #      legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
  #      axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)





plot(tauSeq, beta.scqr[1, ], type = "l", ylim = c(-0.035, 0.04))
lines(tauSeq, beta.cqr[1, ], type = "l", col = "red")
plot(tauSeq, beta.scqr[2, ], type = "l", ylim = c(-1.9, -0.25))
lines(tauSeq, beta.cqr[2, ], type = "l", col = "red")
plot(tauSeq, beta.scqr[3, ], type = "l", ylim = c(-0.65, -0.35))
lines(tauSeq, beta.cqr[3, ], type = "l", col = "red")
plot(tauSeq, beta.scqr[4, ], type = "l", ylim = c(0.8, 2.9))
lines(tauSeq, beta.cqr[4, ], type = "l", col = "red")
plot(tauSeq, beta.scqr[5, ], type = "l", ylim = c(-5.2, 2.2))
lines(tauSeq, beta.cqr[5, ], type = "l", col = "red")



## Highd data
Rcpp::sourceCpp("src/hdscqr.cpp")
dat = read.table("~/Dropbox/Conquer/SCQR/real_data/GSE68465.txt", header = FALSE)
index = which(is.na(dat[1, ]))
dat = dat[, -index]
X = t(as.matrix(dat[3:22285, 2:444]))
censor = as.numeric(dat[1, 2:444] == "vital_status: Dead")
Y = rep(NA, 443)
for (i in 1:443) {
  Y[i] = as.numeric(unlist(strsplit(dat[2, i + 1], " "))[2])
}
index = which(is.na(Y))
Y = Y[-index]
censor = censor[-index]
X = X[-index, ]
n = nrow(X)
p = ncol(X)
X = matrix(as.numeric(X), n, p)
1 - sum(censor) / n  ##censor rate 46.6%
rm(dat)
tauSeq = seq(0.1, 0.7, by = 0.01)
grid = seq(0.1, 0.71, by = 0.01)
nTau = length(tauSeq)
m = length(tauSeq)
lambdaSeq = exp(seq(log(0.14), log(0.06), length.out = 50))
h = 0.5 * (log(p) / n)^(1/4)
setRecord = matrix(0, p, 50)
time = rep(0, 50)

for (i in 1:50) {
  ### lasso
  start = Sys.time()
  beta.lasso = SqrLasso(X, censor, Y, lambdaSeq[i], tauSeq, h)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  activeSet = getSet(beta.lasso, m)
  uniSet = activeSet$union
  if (length(uniSet) > 0) {
    setRecord[uniSet, i] = 1
  }
  
  ### scad
  start = Sys.time()
  beta.scad = SqrScad(X, censor, Y, lambdaSeq[i], tauSeq, h)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  activeSet = getSet(beta.scad, m)
  uniSet = activeSet$union
  if (length(uniSet) > 0) {
    setRecord[uniSet, i] = 1
  }
  
  ### mcp
  start = Sys.time()
  beta.mcp = SqrMcp(X, censor, Y, lambdaSeq[i], tauSeq, h)
  end = Sys.time()
  time[i] = as.numeric(difftime(end, start, units = "secs"))
  activeSet = getSet(beta.mcp, m)
  uniSet = activeSet$union
  if (length(uniSet) > 0) {
    setRecord[uniSet, i] = 1
  }
}

#### cqr, for just one lambda
Z = cbind(1, X)
start = Sys.time()
beta.cqr = quantproc(Y, Z, censor, tauSeq, 0.1)
end = Sys.time()
time = as.numeric(difftime(end, start, units = "secs"))
activeSet = getSet(beta.cqr, m)
uniSet = activeSet$union



### read data
lambdaSeq = exp(seq(log(0.14), log(0.06), length.out = 50))
setwd("~/Dropbox/Conquer/SCQR/Code")
rec = as.matrix(read.csv("real/set_lasso.csv")[, -1])
set1 = rec[1:22283, ]
time1 = rec[22284, ]
plot(1:50, colSums(set1), type = "l")

rec = as.matrix(read.csv("real/set_scad.csv")[, -1])
set2 = rec[1:22283, ]
time2 = rec[22284, ]
plot(1:50, colSums(set2), type = "l")

rec = as.matrix(read.csv("real/set_mcp.csv")[, -1])
set3 = rec[1:22283, ]
time3 = rec[22284, ]
plot(1:50, colSums(set3), type = "l")

which(set1[, 25] == 1)
which(set2[, 25] == 1)
which(set3[, 25] == 1)

