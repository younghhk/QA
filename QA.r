# ======================================================================
# Quantile-Adaptive Sure Independence Screening (QaSIS)
# He, X., Wang, L., & Hong, H.G. (2013), Annals of Statistics
#
# AUTHOR (original): Hyokyoung Grace Hong, Xuming He, Lan Wang
# This version: cleaned, documented, style-harmonized, with basic checks
#
# Inputs:
#   - y    : numeric response vector (uncensored)
#   - x    : numeric matrix of covariates, dimension n x p
#   - tau  : quantile of interest in (0,1)
#
# Output:
#   - numeric vector of length p: marginal QaSIS statistic for each covariate
# ======================================================================

library(quantreg)
library(splines)
library(MASS)
library(survival)

# ----------------------------------------------------------------------
# QaSIS for uncensored data
# ----------------------------------------------------------------------

QaSIS <- function(y, x, tau = 0.5) {
  # Basic checks
  if (!is.numeric(y)) stop("y must be numeric.")
  if (!is.matrix(x))  x <- as.matrix(x)
  if (!is.numeric(x)) stop("x must be numeric or coercible to numeric.")
  if (length(y) != nrow(x)) stop("Length of y must match nrow(x).")
  if (tau <= 0 || tau >= 1) stop("tau must be strictly between 0 and 1.")
  
  n <- length(y)
  p <- ncol(x)
  
  # Center outcome at its tau-quantile
  y_centered <- y - as.numeric(quantile(y, tau, na.rm = TRUE))
  
  # Standardize predictors (column-wise)
  x_scaled <- scale(x)
  
  fit_stat <- numeric(p)
  
  for (j in seq_len(p)) {
    xj <- x_scaled[, j]
    
    # Two interior knots at empirical 1/3 and 2/3 quantiles, degree = 1
    knots_j <- as.numeric(quantile(xj, probs = c(1/3, 2/3), na.rm = TRUE))
    # bs() will handle boundary knots automatically
    a_j <- bs(xj, knots = knots_j, degree = 1)
    
    # Marginal quantile regression
    # silence rq warnings (e.g. perfect fit) to avoid clutter
    b_j <- suppressWarnings(
      rq(y_centered ~ a_j, tau = tau)
    )
    
    # QaSIS statistic: average of squared fitted values
    fit_stat[j] <- mean(fitted(b_j)^2, na.rm = TRUE)
  }
  
  return(fit_stat)
}

# ----------------------------------------------------------------------
# QaSIS for right-censored survival data
#   time  : observed time = min(T, C)
#   delta : event indicator, 1 = event, 0 = censored
#   x     : covariate matrix
#   tau   : quantile of interest for failure time
#
# Uses inverse probability of censoring weights via Kaplan-Meier for C.
# ----------------------------------------------------------------------

QaSIS_surv <- function(x, time, delta, tau = 0.5) {
  # Basic checks
  if (!is.numeric(time))  stop("time must be numeric.")
  if (!is.numeric(delta)) stop("delta must be numeric (0/1).")
  if (!is.matrix(x))      x <- as.matrix(x)
  if (length(time) != nrow(x) || length(delta) != length(time)) {
    stop("Lengths of time, delta, and nrow(x) must match.")
  }
  if (tau <= 0 || tau >= 1) stop("tau must be strictly between 0 and 1.")
  
  N <- length(time)
  p <- ncol(x)
  
  # Kaplan–Meier for failure time to get tau-quantile (medy in original code)
  surv_y <- survfit(Surv(time, delta) ~ 1)
  # find smallest t with S(t) < tau
  idx_tau <- which(surv_y$surv < tau)
  if (length(idx_tau) == 0L) {
    stop("Cannot find a time point where survival falls below tau.")
  }
  med_y <- surv_y$time[min(idx_tau)]
  
  # Kaplan–Meier for censoring distribution (1 - delta)
  surv_c <- survfit(Surv(time, 1 - delta) ~ 1)
  
  # IPCW weights: w_i = delta_i / G_hat(time_i)
  w <- numeric(N)
  for (i in seq_len(N)) {
    # find largest time in KM for C <= observed time[i]
    idx_i <- max(which(round(surv_c$time, 4) <= round(time[i], 4)))
    G_hat  <- surv_c$surv[idx_i]
    w[i]   <- ifelse(G_hat > 0, delta[i] / G_hat, 0)
  }
  
  # Screening statistic for each covariate
  fit_stat <- numeric(p)
  for (k in seq_len(p)) {
    xk <- x[, k]
    
    # B-spline basis with df = 3 (as in original code)
    pix_k <- bs(xk, df = 3)
    
    # Weighted quantile regression of time on pix_k
    beta_k <- suppressWarnings(
      rq(time ~ pix_k, tau = tau, weights = w)
    )
    
    mu_k <- as.numeric(predict(beta_k))
    
    # QaSIS statistic: mean squared deviation from tau-quantile of T
    fit_stat[k] <- mean((mu_k - med_y)^2, na.rm = TRUE)
  }
  
  return(fit_stat)
}

# ----------------------------------------------------------------------
# Example usage (uncensored)
# ----------------------------------------------------------------------

simul_dat_example <- function(N, p, seed = 100) {
  if (!is.null(seed)) set.seed(seed)
  
  active <- 1:4
  
  # AR(1)-type covariance with rho = 0.8
  Sigma1 <- diag(p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i < j) {
        Sigma1[i, j] <- 0.8^abs(i - j)
        Sigma1[j, i] <- Sigma1[i, j]
      }
    }
  }
  
  # Nonlinear functions as in the AoS paper
  g1 <- function(x) (x)
  g2 <- function(x) ((2 * x - 1)^2)
  g3 <- function(x) (sin(2 * pi * x) / (2 - sin(2 * pi * x)))
  g4 <- function(x) (
    0.1 * sin(2 * pi * x) +
      0.2 * cos(2 * pi * x) +
      0.3 * sin(2 * pi * x)^2 +
      0.4 * cos(2 * pi * x)^3 +
      0.5 * sin(2 * pi * x)^3
  )
  
  X <- mvrnorm(N, mu = rep(0, p), Sigma = Sigma1)
  eps <- rnorm(N, mean = 0, sd = sqrt(1.74))
  
  Y <- 5 * g1(X[, 1]) +
    3 * g2(X[, 2]) +
    4 * g3(X[, 3]) +
    6 * g4(X[, 4]) +
    eps
  
  list(x = X, y = Y, active = active)
}

