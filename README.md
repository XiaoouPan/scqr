# scqr

Smoothed Censored Quantile Regression Process

## Description

We approach the globally-concerned censored quantile regression process with a smoothing mechanism for efficient computation.
In the low dimensional regime, the regression process is formulated as solving a sequence of smoothed estimating equations (SEE), which can be done via a quasi-Newton method.
Coordinatewise confidence intervals of coefficients can be constructed by multiplier bootstrap.
In the high dimensional regime, the sparse learning problem is solved by iteratively reweighted *&ell;<sub>1</sub>*-regularized regression, and each *&ell;<sub>1</sub>*-regularized regression is solved by a local majorize-minimize algorithm.


## References

He, X., Pan, X., Tan, K. M., and Zhou, W.-X. (2022). Smoothed quantile regression with large-scale inference. *J. Econometrics*, to appear. [Paper](https://doi.org/10.1016/j.jeconom.2021.07.010)

Koenker, R. and Bassett, G. (1978). Regression quantiles. *Econometrica* **46** 33-50. [Paper](https://www.jstor.org/stable/1913643?seq=1#metadata_info_tab_contents)

Peng, L. and Huang, Y. (2008). Survival analysis with quantile regression models. *J. Am. Stat. Assoc.* **103** 637–649. [Paper](https://doi.org/10.1198/016214508000000355)

Zheng, Q., Peng, L. and He, X. (2018). High dimensional censored quantile regression. *Ann. Statist.* **46** 308-343. [Paper](https://doi.org/10.1214/17-AOS1551)

