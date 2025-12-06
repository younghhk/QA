[![Back to Hub](https://img.shields.io/badge/⬅️%20Back%20to%20Hub-2962FF?style=for-the-badge)](https://github.com/younghhk/NCI)

# Quantile-Adaptive Screening
This repository introduces a quantile-adaptive, model-free framework for nonlinear variable screening with high-dimensional heterogeneous data. The method is designed for settings where the set of important variables may differ across the outcome distribution and the number of candidate predictors is large.​

## Key features
* Quantile-specific effects: the active set of variables is allowed to vary across quantiles (for example, lower vs upper tail of the outcome).​

* Model-free: avoids specifying a full regression model in a high-dimensional covariate space.​

* Nonlinear screening: uses spline-based quantile regression to capture nonlinear marginal associations at a chosen quantile.​

* Theoretical support: under mild conditions, the procedure has a sure screening property in ultra-high dimensions.​

* Censored outcomes: extends naturally to right-censored survival data while retaining the same screening principle.​

## Example usage

```
source("QA.r")
## ------------- Uncensored example ----------------------
dat    <- simul_dat_example(N = 200, p = 1000, seed = 100)
x      <- dat$x
y      <- dat$y
active <- dat$active

out <- QaSIS(y = y, x = x, tau = 0.5)
rank(-out)[active]   # ranks of true active variables
[1] 2 1 3 4

## -------------- Survival example ------------------------
out_surv <- QaSIS.surv(
  x       = x,
  time    = time,
  delta   = delta,
  tau     = 0.5,
  w_trunc = 20,
  df_bs   = 3
)

rank(-out_surv)[active]   # ranks of true active variables
[1] 4 1 2 3
```


## Reference
He X, Wang L, Hong HG. [Quantile-adaptive model-free variable screening for high-dimensional heterogeneous data.](https://projecteuclid.org/journals/annals-of-statistics/volume-41/issue-1/Quantile-adaptive-model-free-variable-screening-for-high-dimensional-heterogeneous/10.1214/13-AOS1087.pdf) Annals of Statistics. 2013;41(1):342–369.​
