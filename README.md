# scqr

Smoothed Censored Quantile Regression Process

## Description

We approach the globally-concerned censored quantile regression process with a smoothing mechanism. 
In the low dimensional regime, the regression process is formulated as solving a sequence of smoothed estimating equations (SEE), which can be done via a quasi-Newton method.
Coordinatewise confidence intervals of coefficients can be constructed by multiplier bootstrap.
In the high dimensional regime, the sparse learning problem is solved by iteratively reweighted *&ell;<sub>1</sub>*-regularized regression.