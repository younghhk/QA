#######################################################################
# Quantile-Adaptive Sure Independence Screening (QaSIS)
# He, X., Wang, L., & Hong, H.G. (2013), Annals of Statistics
#
# AUTHOR (original): Hyokyoung Grace Hong, Xuming He, Lan Wang
# This version: cleaned and organized for easy reuse
#
# CONTENTS
#   1.  Packages
#   2.  Core functions
#         - QaSIS        (uncensored response)
#         - QaSIS.surv   (right-censored survival data)
#   3.  Simulation helpers
#         - simul_dat_example     (uncensored)
#         - simul_surv_example    (survival)
#   4.  Example usage
#######################################################################

#######################################################################
# 1. PACKAGES
#######################################################################

library(quantreg)
library(splines)
library(MASS)
library(survival)

#######################################################################
# 2. CORE FUNCTIONS
#######################################################################

# --------------------------------------------------------------------
# 2.1 QaSIS for uncensored data
# --------------------------------------------------------------------
# Inputs:
#   y    : numeric response vector (length n)
#   x    : numeric matrix of covariates, dimension n x p
#   tau  : quantile of interest in (0,1)
#
# Output:
#   numeric vector of length p (QaSIS statistic for each covariate)
# --------------------------------------------------------------------

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
    a_j     <- bs(xj, knots = knots_j, degree = 1)
    
    # Marginal quantile regression
    b_j <- suppressWarnings(
      rq(y_centered ~ a_j, tau = tau)
    )
    
    # QaSIS statistic: average of squared fitted values
    fit_stat[j] <- mean(fitted(b_j)^2, na.rm = TRUE)
  }
  
  fit_stat
}

# --------------------------------------------------------------------
# 2.2 QaSIS for right-censored survival data
# --------------------------------------------------------------------
# Inputs:
#   time   : observed time = min(T, C)
#   delta  : event indicator, 1 = event, 0 = censored
#   x      : covariate matrix (n x p)
#   tau    : quantile of interest for failure time T
#   w_trunc: truncation level for IPCW weights (numeric, default 20)
#   df_bs  : degrees of freedom for B-splines (default 3)
#
# Output:
#   numeric vector of length p (QaSIS–survival statistic)
#
# Uses inverse probability of censoring weights via Kaplan–Meier for C.
# --------------------------------------------------------------------

QaSIS.surv <- function(x, time, delta, tau = 0.5,
                       w_trunc = 20,
                       df_bs   = 3) {
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
  
  # KM for T to get tau-quantile on original time scale
  surv_y  <- survfit(Surv(time, delta) ~ 1)
  idx_tau <- which(surv_y$surv < tau)
  if (length(idx_tau) == 0L) {
    stop("Cannot find a time point where survival falls below tau.")
  }
  med_y <- surv_y$time[min(idx_tau)]
  
  # KM for censoring distribution G
  surv_c <- survfit(Surv(time, 1 - delta) ~ 1)
  t_c    <- surv_c$time
  G_hat  <- surv_c$surv
  
  # IPCW weights: w_i = delta_i / G_hat(time_i)
  w <- numeric(N)
  for (i in seq_len(N)) {
    idx_i <- which(t_c <= time[i])
    if (length(idx_i) == 0L) {
      g_i <- G_hat[1]
    } else {
      g_i <- G_hat[max(idx_i)]
    }
    if (is.na(g_i) || g_i <= 0) {
      w[i] <- 0
    } else {
      w[i] <- delta[i] / g_i
    }
  }
  
  # Truncate large weights for numerical stability
  w <- pmin(w, w_trunc)
  
  # QaSIS statistic for each covariate
  fit_stat <- numeric(p)
  for (k in seq_len(p)) {
    xk    <- x[, k]
    xk_sc <- scale(xk)             # standardize before spline
    
    pix_k <- bs(xk_sc, df = df_bs) # B-spline basis
    
    # Weighted quantile regression of time on pix_k
    beta_k <- suppressWarnings(
      rq(time ~ pix_k, tau = tau, weights = w)
    )
    
    mu_k <- as.numeric(predict(beta_k))
    
    # QaSIS statistic: mean squared deviation from tau-quantile of T
    fit_stat[k] <- mean((mu_k - med_y)^2, na.rm = TRUE)
  }
  
  fit_stat
}

#######################################################################
# 3. SIMULATION HELPERS
#######################################################################

# --------------------------------------------------------------------
# 3.1 Uncensored data (AoS-style example)
# --------------------------------------------------------------------

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
  
  X   <- mvrnorm(N, mu = rep(0, p), Sigma = Sigma1)
  eps <- rnorm(N, mean = 0, sd = sqrt(1.74))
  
  Y <- 5 * g1(X[, 1]) +
    3 * g2(X[, 2]) +
    4 * g3(X[, 3]) +
    6 * g4(X[, 4]) +
    eps
  
  list(x = X, y = Y, active = active)
}

# --------------------------------------------------------------------
# 3.2 Survival data example
# --------------------------------------------------------------------
# Linear signal on original time scale, AR(1) covariance, exponential
# censoring tuned to target censor_rate.
# --------------------------------------------------------------------

simul_surv_example <- function(N, p, censor_rate = 0.2, seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  
  active <- 1:4
  
  # AR(1) covariance, moderate correlation
  rho    <- 0.5
  Sigma1 <- diag(p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i < j) {
        Sigma1[i, j] <- rho^abs(i - j)
        Sigma1[j, i] <- Sigma1[i, j]
      }
    }
  }
  
  # Covariates
  X <- mvrnorm(N, mu = rep(0, p), Sigma = Sigma1)
  
  # Simple linear signal on original time scale
  g1 <- function(x) x
  g2 <- function(x) x
  g3 <- function(x) x
  g4 <- function(x) x
  
  mu_raw <- 3.5 * g1(X[, 1]) +
    4.0 * g2(X[, 2]) +
    3.8 * g3(X[, 3]) +
    4.5 * g4(X[, 4])
  
  # Add noise
  eps   <- rnorm(N, mean = 0, sd = 1.0)
  T_raw <- mu_raw + eps
  
  # Shift to make times positive
  shift  <- max(0.1 - min(T_raw), 0)
  T_true <- T_raw + shift
  
  # Exponential censoring tuned to target censor_rate
  lambda_c <- -log(1 - censor_rate) / median(T_true)
  C        <- rexp(N, rate = lambda_c)
  
  time  <- pmin(T_true, C)
  delta <- as.numeric(T_true <= C)
  
  list(
    x      = X,
    time   = time,
    delta  = delta,
    active = active,
    T_true = T_true,
    C      = C
  )
}

#######################################################################
# 4. EXAMPLE USAGE
#######################################################################
# You can comment out this section when using in a package.
#######################################################################

## ------------------- 4.1 Uncensored example -------------------------

# dat    <- simul_dat_example(N = 200, p = 1000, seed = 100)
# x      <- dat$x
# y      <- dat$y
# active <- dat$active
#
# out <- QaSIS(y = y, x = x, tau = 0.5)
# rank(-out)[active]   # ranks of true active variables

## ------------------- 4.2 Survival example ---------------------------

# sdat <- simul_surv_example(
#   N           = 2000,
#   p           = 500,
#   censor_rate = 0.2,
#   seed        = 2025
# )
#
# x      <- sdat$x
# time   <- sdat$time
# delta  <- sdat$delta
# active <- sdat$active
#
# out_surv <- QaSIS.surv(
#   x       = x,
#   time    = time,
#   delta   = delta,
#   tau     = 0.5,
#   w_trunc = 20,
#   df_bs   = 3
# )
#
# rank(-out_surv)[active]   # ranks of true active variables


