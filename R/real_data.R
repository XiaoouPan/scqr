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

start = Sys.time()
#list = scqrGaussInf(X, Y, censor, tauSeq, B = 1000)
list = scqrGauss(X, Y, censor, tauSeq)
end = Sys.time()
as.numeric(difftime(end, start, units = "secs"))
beta.scqr = list$coeff[-1, ]

response = Surv(Y, censor, type = "right")
list = crq(response ~ X, method = "PengHuang", grid = grid)
beta.cqr = list$sol[3:7, ]

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

beta.boot = list$boot
ci.list = getPivCI(beta.hat[, nTau], beta.boot, alpha)
ci.per = ci.list$perCI
ci.piv = ci.list$pivCI
ci.norm = getNormCI(beta.hat[, nTau], rowSds(beta.boot), z)



## Highd data
Rcpp::sourceCpp("src/hdscqr.cpp")
dat = read.table("~/Dropbox/Conquer/SCQR/real_data/GSE68465.txt", header = FALSE)
index = which(is.na(dat[1, ]))
dat = dat[, -index]
dim(dat)
X = t(as.matrix(dat[3:22285, 2:444]))
censor = as.numeric(dat[1, 2:444] == "vital_status: Dead")
Y = rep(NA, n)
for (i in 1:n) {
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

### Run code



