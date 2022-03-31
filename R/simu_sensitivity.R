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

#### Growing dimension and sample size
n = 5000
pseq = seq(20, 100, by = 20)
l = length(pseq)
tauSeq = seq(0.05, 0.8, by = 0.05)
grid = seq(0.05, 0.85, by = 0.05)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
M = 100
coef1 = coef2 = coef3 =  matrix(NA, M, l)

## Quantil index of interest
index1 = 6  #tau = 0.3
index2 = 10 #tau = 0.5
index3 = 14 #tau = 0.7

pb = txtProgressBar(style = 3)
for (j in 1:l) {
  p = pseq[j]
  #h = 0.5 * ((p + log(n)) / n)^(0.4)
  h2 = max(0.5 * (log(p) / n)^(1/4), 0.05) ## bandwidth for high-d
  Sigma = toeplitz(0.5^(0:(p - 1)))
  for (i in 1:M) {
    set.seed((j - 1) * M + i)
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
    list = scqrGauss(X, Y, censor, tauSeq)
    coef1[i, j] = norm(list$coeff[, index2] - betaMat[, index2], "2")

    ## Peng and Huang
    list = crq(response ~ X, method = "PengHuang", grid = grid)
    tt = ncol(list$sol)
    if (tt >= index2) {
      coef2[i, j] = norm(list$sol[2:(p + 2), index2] - betaMat[, index2], "2")
    }

    ## Smoothed CQR with misspecified h
    list = scqrGauss(X, Y, censor, tauSeq, h = h2)
    coef3[i, j] = norm(list$coeff[, index2] - betaMat[, index2], "2")
    
    setTxtProgressBar(pb, ((j - 1) * M + i) / (l * M))
  }
}


#write.csv(prop, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Growing/prop.csv")
#write.csv(time1, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Growing/time1.csv")
#write.csv(time2, "~/Dropbox/Conquer/censoredQR/Code/Simulation/Growing/time2.csv")
write.csv(coef1, "~/Dropbox/Conquer/SCQR/AOS_rev/Simulation/Sensitivity/coef1_hetero.csv")
write.csv(coef2, "~/Dropbox/Conquer/SCQR/AOS_rev/Simulation/Sensitivity/coef2_hetero.csv")
write.csv(coef3, "~/Dropbox/Conquer/SCQR/AOS_rev/Simulation/Sensitivity/coef3_hetero.csv")
#write.csv(coef4, "~/Dropbox/Conquer/SCQR/AOS_rev/Simulation/Growing/coef4_hetero.csv")
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
diff = as.matrix(read.csv("Simulation/growing/diff_homo.csv")[, -1])


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
mean1 = colMeans(coef3, na.rm = TRUE)
mean2 = colMeans(coef4, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Our method}", l), rep("\\texttt{Peng} \\& \\texttt{Huang}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Our method}"))

setwd("~/Dropbox/Conquer/SCQR/Code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid")) +
  theme_bw() + xlab("Sample size") + ylab("Estimation error at $\\tau = 0.5$") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  #theme(legend.position = c(0.32, 0.88), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
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



###### Sensitivity analysis (for JASA rebuttal) with growing scale

setwd("~/Dropbox/Conquer/SCQR/Code")
coef = as.matrix(read.csv("Simulation/sensitivity/coef_grow_homo.csv")[, -1])
coef1 = coef[1:500, ]
coef2 = coef[501:1000, ]
coef3 = coef[1001:1500, ]
coef4 = coef[1501:2000, ]
coef5 = coef[2001:2500, ]
time = as.matrix(read.csv("Simulation/sensitivity/time_grow_hetero.csv")[, -1])
time1 = time[1:500, ]
time2 = time[501:1000, ]
time3 = time[1001:1500, ]
time4 = time[1501:2000, ]
time5 = time[2001:2500, ]


###### Sensitivity analysis (for AOS revision) with growing scale

coef1 = as.matrix(read.csv("~/Dropbox/Conquer/SCQR/AOS_rev/Simulation/sensitivity/coef3_hetero.csv")[, -1])

coef = as.matrix(read.csv("~/Dropbox/Conquer/SCQR/code/Simulation/growing/coef_hetero.csv")[, -1])
coef = cbind(coef, as.matrix(read.csv("~/Dropbox/Conquer/SCQR/code/Simulation/growing/coef_hetero_15.csv")[, -c(1, 2)]))
coef = cbind(coef, as.matrix(read.csv("~/Dropbox/Conquer/SCQR/code/Simulation/growing/coef_hetero_18.csv")[, -1]))
coef = cbind(coef, as.matrix(read.csv("~/Dropbox/Conquer/SCQR/code/Simulation/growing/coef_hetero_20.csv")[, -1]))
coef6 = coef[1:500, ]
coef7 = coef[501:1000, ]
time = as.matrix(read.csv("Simulation/growing/time_hetero.csv")[, -1])
time = cbind(time, as.matrix(read.csv("Simulation/growing/time_hetero_15.csv")[, -c(1, 2)]))
time = cbind(time, as.matrix(read.csv("Simulation/growing/time_hetero_18.csv")[, -1]))
time = cbind(time, as.matrix(read.csv("Simulation/growing/time_hetero_20.csv")[, -1]))
time6 = time[1:500, ]
time7 = time[501:1000, ]


### Estimation error at a certain quantile
mean1 = colMeans(coef1, na.rm = TRUE)
mean6 = colMeans(coef6, na.rm = TRUE)
mean7 = colMeans(coef7, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean6), cbind(nseq, mean7))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Misspecified bandwidth}", l), rep("\\texttt{Optimal bandwidth}", l), rep("\\texttt{Peng} \\& \\texttt{Huang}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Misspecified bandwidth}", "\\texttt{Optimal bandwidth}"))
#dat$color = c(rep("aquamarine", l), rep("bisque", l), rep("blueviolet", l), rep("cadetblue1", l), rep("coral", l), rep("brown1", l), rep("cyan", l))
#dat$color = factor(dat$color, levels = c("cyan", "aquamarine", "bisque", "blueviolet", "cadetblue1", "coral", "brown1"))

setwd("~/Dropbox/Conquer/SCQR/code")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("twodash", "dashed", "solid")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Estimation error at $\\tau = 0.7$") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  #theme(legend.position = c(0.7, 0.15), legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"),
  #    legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
  #    axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



## Time plot
mean1 = colMeans(time1, na.rm = TRUE)
mean2 = colMeans(time2, na.rm = TRUE)
mean3 = colMeans(time3, na.rm = TRUE)
mean4 = colMeans(time4, na.rm = TRUE)
mean5 = colMeans(time5, na.rm = TRUE)
mean6 = colMeans(time6, na.rm = TRUE)
mean7 = colMeans(time7, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2), cbind(nseq, mean3), cbind(nseq, mean4), cbind(nseq, mean5), cbind(nseq, mean6), cbind(nseq, mean7))
dat = as.data.frame(dat)
colnames(dat) = c("size", "time")
dat$type = c(rep("\\texttt{Bandwidth}$=0.05$", l), rep("\\texttt{Bandwidth}$=0.1$", l), rep("\\texttt{Bandwidth}$=0.15$", l),
             rep("\\texttt{Bandwidth}$=0.2$", l), rep("\\texttt{Bandwidth}$=0.25$", l), rep("\\texttt{Optimal bandwidth}", l), 
             rep("\\texttt{Peng} \\& \\texttt{Huang}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Peng} \\& \\texttt{Huang}", "\\texttt{Bandwidth}$=0.05$", "\\texttt{Bandwidth}$=0.1$",
                                       "\\texttt{Bandwidth}$=0.15$", "\\texttt{Bandwidth}$=0.2$", "\\texttt{Bandwidth}$=0.25$", 
                                       "\\texttt{Optimal bandwidth}"))

tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = time)) +
  geom_line(aes(y = time, color = type, linetype = type), size = 1.5) + 
  scale_linetype_manual(values = c("twodash", rep("solid", 5), "dashed")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Elapsed time (in seconds)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  #theme(legend.position = c(0.3, 0.68), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"),
  #      legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
  #      axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



