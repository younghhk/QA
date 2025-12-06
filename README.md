
# Quantile-adpative screening

```

# ----------------------------------------------------------------------
# Example usage (uncensored)
# ----------------------------------------------------------------------
source(Qasis.r)
## Generate data
dat    <- simul_dat_example(N = 200, p = 1000, seed = 100)
x      <- dat$x
y      <- dat$y
active <- dat$active

## Compute QaSIS statistics at tau = 0.5
out <- QaSIS(y = y, x = x, tau = 0.5)

## Ranks of the true active variables (should be near the top)
rank(-out)[active]
```

```
# ----------------------------------------------------------------------
# Example usage (censored)
# ----------------------------------------------------------------------
source(Qasissv.r)
## Generate data
dat    <- simul_dat_example(N = 200, p = 1000, seed = 100)
x      <- dat$x
y      <- dat$y
active <- dat$active

## Compute QaSIS statistics at tau = 0.5
out <- QaSIS(y = y, x = x, tau = 0.5)

## Ranks of the true active variables (should be near the top)
rank(-out)[active]
```

## References
