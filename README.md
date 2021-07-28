# scqr

Smoothed Censored Quantile Regression Process

## Description

We approach the globally-concerned censored quantile regression process with a smoothing mechanism. 
In the low dimensional regime, the regression process is formulated as solving a sequence of smoothed estimating equations (SEE), which can be done via a quasi-Newton method.
Coordinatewise confidence intervals of coefficients can be constructed by multiplier bootstrap.
In the high dimensional regime, the sparse learning problem is solved by iteratively reweighted *&ell;<sub>1</sub>*-regularized regression, and each *&ell;<sub>1</sub>*-regularized regression is solved by a local MM algorithm.


## References

He, X., Pan, X., Tan, K. M., and Zhou, W.-X. (2020). Smoothed quantile regression with large-scale inference. [Paper](https://www.math.ucsd.edu/~xip024/Papers/sqr_main.pdf)

Koenker, R. and Bassett, G. (1978). Regression quantiles. *Econometrica* **46** 33-50. [Paper](https://www.jstor.org/stable/1913643?seq=1#metadata_info_tab_contents)

Pan, X., Sun, Q. and Zhou, W.-X. (2021). Iteratively reweighted *l<sub>1</sub>*-penalized robust regression. *Electron. J. Stat.* **15** 3287-3348. [Paper](https://doi.org/10.1214/21-EJS1862)

Peng, L. and Huang, Y. (2008). Survival analysis with quantile regression models. *J. Am. Stat. Assoc.* **103** 637–649. [Paper](https://doi.org/10.1198/016214508000000355)

Sanderson, C. and Curtin, R. (2016). Armadillo: A template-based C++ library for linear algebra. *J. Open Source Softw.* **1** 26. [Paper](https://joss.theoj.org/papers/10.21105/joss.00026.pdf)

