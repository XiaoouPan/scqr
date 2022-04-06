###### Simulation with a fixed data scale

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

estError = function(betahat, beta, tauSeq) {
  nTau = min(length(tauSeq), ncol(betahat) + 1)
  diff = betahat - beta
  err = sqrt(colSums((betahat - beta)^2))
  accu = 0
  for (k in 2:nTau) {
    accu = accu + err[k - 1] * (tauSeq[k] - tauSeq[k - 1])
  }
  return (accu)
}


#### Quantile process with fixed scale, hard to visualize
n = 5000
p = n / 50
M = 500
tauSeq = seq(0.05, 0.8, by = 0.05)
grid = seq(0.05, 0.85, by = 0.05)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
coef1 = eff1 = matrix(0, M, nTau)
coef2 = eff2 = matrix(0, M, nTau)
coef3 = eff3 = matrix(0, M, nTau)
time = matrix(0, 3, M)
prop = rep(0, M)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  #X = sqrt(12) * draw.d.variate.uniform(n, p, Sigma) - sqrt(3)
  #Sigma = getSigma(p)
  #X = mvrnorm(n, rep(0, p), Sigma)
  Sigma = getSigma(45)
  X = cbind(mvrnorm(n, rep(0, 45), Sigma), 4 * draw.d.variate.uniform(n, 45, Sigma) - 2, matrix(rbinom(10 * n, 1, c(0.5, 0.5)), n, 10))
  err = rt(n, 2)
  ## Homo
  #beta = runif(p, -2, 2)
  #betaMat = rbind(beta0, matrix(beta, p, nTau))
  #logT = X %*% beta + err
  ## Hetero
  X[, 1] = abs(X[, 1])
  beta = runif(p - 1, -2, 2)
  betaMat = rbind(rep(0, nTau), beta0, matrix(beta, p - 1, nTau))
  logT = X[, 1] * err + X[, -1] %*% beta
  w = sample(1:3, n, prob = c(1/3, 1/3, 1/3), replace = TRUE)
  logC = (w == 1) * rnorm(n, 0, 4) + (w == 2) * rnorm(n, 5, 1) + (w == 3) * rnorm(n, 10, 0.5)
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
  #eff1[i, ] = list$coeff[1, ]
  eff1[i, ] = list$coeff[2, ]
  
  ## Peng and Huang
  start = Sys.time()
  list = crq(response ~ X, method = "PengHuang", grid = grid)
  end = Sys.time()
  time[2, i] = as.numeric(difftime(end, start, units = "secs"))
  tt = ncol(list$sol)
  coef2[i, 1:tt] = sqrt(colSums((list$sol[2:(p + 2), ] - betaMat[, 1:tt])^2))
  #eff2[i, 1:tt] = list$sol[2, ]
  eff2[i, 1:tt] = list$sol[3, ]
  
  ## Portnoy
  start = Sys.time()
  list = crq(response ~ X, method = "Portnoy", grid = tauSeq)
  end = Sys.time()
  time[3, i] = as.numeric(difftime(end, start, units = "secs"))
  coef3[i, ] = sqrt(colSums((list$sol[2:(p + 2), 2:17] - betaMat)^2))
  #eff3[i, ] = list$sol[2, 2:17]
  eff3[i, ] = list$sol[3, 2:17]
  
  setTxtProgressBar(pb, i / M)
}




setwd("~/Dropbox/Conquer/SCQR/Code")
coe.data = as.matrix(read.csv("Simulation/sensitivity/coef_homo.csv")[, -1])
coef1 = coe.data[1:500, ]
coef2 = coe.data[501:1000, ]
coef4 = coe.data[1001:1500, ]
coef5 = coe.data[1501:2000, ]
coe.data = as.matrix(read.csv("Simulation/coef_homo.csv")[, -1])
coef3 = coe.data[1:500, ]
coef6 = coe.data[501:1000, ]

### Estimation plots
#sd1 = colSds(coef1)
mean1 = colMeans(coef1)
#low = mean1 - sd1
#upp = mean1 + sd1
#dat = cbind(tauSeq, mean1, low, upp)
#sd2 = colSds(coef2)
mean2 = colMeans(coef2)
#low = mean2 - sd2
#upp = mean2 + sd2
mean4 = colMeans(coef4)
mean5 = colMeans(coef5)
mean3 = colMeans(coef3)
mean6 = colMeans(coef6)
dat = rbind(cbind(tauSeq, mean1), cbind(tauSeq, mean2), cbind(tauSeq, mean3), cbind(tauSeq, mean4), cbind(tauSeq, mean5), cbind(tauSeq, mean6))
#sd3 = colSds(coef3)
#mean3 = colMeans(coef3)
#low = mean3 - sd3
#upp = mean3 + sd3
#dat = rbind(dat, cbind(tauSeq, mean3, low, upp))
dat = as.data.frame(dat)
colnames(dat) = c("quantile", "coef")
dat$type = c(rep("\\texttt{Bandwidth} $h_1$", nTau), rep("\\texttt{Bandwidth} $h_2$", nTau), rep("\\texttt{Bandwidth} $h_3$", nTau), 
             rep("\\texttt{Bandwidth} $h_4$", nTau), rep("\\texttt{Bandwidth} $h_5$", nTau), rep("\\texttt{Peng} \\& \\texttt{Huang}", nTau))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Bandwidth} $h_1$", "\\texttt{Bandwidth} $h_2$",
                                       "\\texttt{Bandwidth} $h_3$", "\\texttt{Bandwidth} $h_4$", "\\texttt{Bandwidth} $h_5$"))

tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = quantile, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 1) + 
  scale_linetype_manual(values = c("solid", rep("twodash", 5))) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Quantile level $\\tau$") + ylab("Estimation error in $||\\cdot||_2$") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.65, 0.7), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(0.9, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



### Quantile effects plots
setwd("~/Dropbox/Conquer/SCQR/Code")
eff.data = as.matrix(read.csv("Simulation/sensitivity/eff_hetero.csv")[, -1])
eff1 = eff.data[1:500, ]
eff2 = eff.data[501:1000, ]
eff4 = eff.data[1001:1500, ]
eff5 = eff.data[1501:2000, ]
eff.data = as.matrix(read.csv("Simulation/eff_hetero.csv")[, -1])
eff3 = eff.data[1:500, ]
eff6 = eff.data[501:1000, ]


#sd1 = colSds(eff1)
mean1 = colMeans(eff1)
#low = mean1 - sd1
#upp = mean1 + sd1
#dat = cbind(tauSeq, mean1, low, upp)
#sd2 = colSds(eff2)
mean2 = colMeans(eff2)
#low = mean2 - sd2
#upp = mean2 + sd2
mean4 = colMeans(eff4)
mean5 = colMeans(eff5)
mean3 = colMeans(eff3)
mean6 = colMeans(eff6)
dat = rbind(cbind(tauSeq, mean1), cbind(tauSeq, mean2), cbind(tauSeq, mean3), cbind(tauSeq, mean4), cbind(tauSeq, mean5), cbind(tauSeq, mean6))
dat = rbind(dat, cbind(tauSeq, beta0))
dat = as.data.frame(dat)
colnames(dat) = c("quantile", "eff")
dat$type = c(rep("\\texttt{Bandwidth} $h_1$", nTau), rep("\\texttt{Bandwidth} $h_2$", nTau), rep("\\texttt{Bandwidth} $h_3$", nTau), 
             rep("\\texttt{Bandwidth} $h_4$", nTau), rep("\\texttt{Bandwidth} $h_5$", nTau), rep("\\texttt{Peng} \\& \\texttt{Huang}", nTau), 
             rep("\\texttt{True effects}", nTau))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Bandwidth} $h_1$", "\\texttt{Bandwidth} $h_2$",
                                       "\\texttt{Bandwidth} $h_3$", "\\texttt{Bandwidth} $h_4$", "\\texttt{Bandwidth} $h_5$", "\\texttt{True effects}"))

tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = quantile, y = eff, color = type)) +
  geom_line(aes(y = eff, color = type, linetype = type), size = 1) + 
  scale_linetype_manual(values = c("solid", rep("twodash", 5), "dashed")) +
  #geom_ribbon(aes(y = eff, ymin = low, ymax = upp, fill = type), alpha = 0.3) + 
  theme_bw() + xlab("Quantile level $\\tau$") + 
  ylab("Estimated quantile effects") + 
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  #theme(legend.position = c(0.78, 0.32), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(0.9, "cm"),
  #      legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
  #      axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)





### Running time plots
setwd("~/Dropbox/Conquer/SCQR/Code")
time = rbind(as.matrix(read.csv("Simulation/sensitivity/time_hetero.csv")[, -1]), as.matrix(read.csv("Simulation/time_hetero.csv")[1:2, -1]))
### For sensitivity analysis
meth = c(rep("$h_1$", M), rep("$h_2$", M), rep("$h_3$", M), rep("$h_4$", M), rep("$h_5$", M), rep("\\texttt{P} \\& \\texttt{H}", M))
method = factor(meth, levels = c("$h_1$", "$h_2$", "$h_3$", "$h_4$", "$h_5$", "\\texttt{P} \\& \\texttt{H}"))
rst = data.frame("time" = c(time[1, ], time[2, ], time[5, ], time[3, ], time[4, ], time[6, ]), "method" = method)
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(rst, aes(x = method, y = time, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Elapsed time (in seconds)") + theme_bw() + 
  #scale_y_continuous(breaks = seq(0, 125, 25)) + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### For main text
meth = c(rep("Our method", 500), rep("Peng \\& Huang", 500))
meth = factor(meth, levels = c("Our method", "Peng \\& Huang"))
rst = data.frame("time" = c(time[1, ], time[2, ]), "method" = meth)
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(rst, aes(x = method, y = time, fill = method)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Elapsed time (in seconds)") + theme_bw() + 
  #scale_y_continuous(breaks = seq(0, 125, 25)) + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

