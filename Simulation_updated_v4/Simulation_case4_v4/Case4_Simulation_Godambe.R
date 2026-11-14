setwd('/Users/idongha/Library/CloudStorage/OneDrive-SimonFraserUniversity(1sfu)/SFU/Research/Deductible ratemaking/Submitted/2nd/sim_dh/261030_case_4')

load(file="Case4_2nd_round.Rdata")

# ==========================================================================
# CASE 4 — self-contained simulation, derived from the user's validated
# Case 1 workflow (October 30, 2026).
# rho_N = -0.6, rho_Y = -0.4. Conditional on observed D and X,
# ground-up claim count and individual severities are generated independently.
# This file fits all FEM/SEM/FDM/SDM models, assesses four premium metrics,
# reports parameter bias/RMSE, and optionally computes selected Godambe SEs.
# EDIT t, RUN_SANDWICH, SANDWICH_REPS, SANDWICH_TARGETS below.
# Requires: R + MASS. Run in a case-specific working directory.
# ==========================================================================

# [CHANGED 1] Only MASS is required; original Rcpp GB2 functions are implemented
# with mathematically equivalent, log-scale Lomax formulas in this same file.
library(MASS)

# Compatibility names preserve the original sim_loglikM_*_cpp() calls.
# The simulation fixes GB2 sigma=alpha1=1, so it is exactly Lomax.
# mu is log(scale), not the severity mean.
pGB2_cpp <- function(y, mu, sigma=1, alpha1=1, alpha2) {
  stopifnot(sigma == 1, alpha1 == 1, alpha2 > 0)
  -expm1(-alpha2 * log1p(y * exp(-as.numeric(mu))))
}

lomax_ll <- function(y, mu, alpha) {
  eta <- as.numeric(mu)
  log(alpha) - eta - (alpha + 1) * log1p(y * exp(-eta))
}

sim_loglikM_ground_GB2_cpp <- function(Claim, sim_covaS, parms) {
  k <- ncol(sim_covaS)
  eta <- drop(sim_covaS %*% parms[seq_len(k)])
  val <- lomax_ll(Claim, eta, parms[k+1L])
  if (any(!is.finite(val))) return(1e100)
  -sum(val)
}

sim_loglikM_censored_GB2_cpp <- function(LossAfterDeductible, Deduct,
                                          sim_covaS, parms) {
  k <- ncol(sim_covaS)
  eta <- drop(sim_covaS %*% parms[seq_len(k)])
  a <- parms[k+1L]
  logsurv <- -a * log1p(Deduct * exp(-eta))
  zeros <- LossAfterDeductible == 0
  ll <- numeric(length(Deduct))
  ll[zeros] <- log(-expm1(logsurv[zeros]))
  if (any(!zeros)) {
    ll[!zeros] <- lomax_ll(LossAfterDeductible[!zeros] + Deduct[!zeros],
                            eta[!zeros], a)
  }
  if (any(!is.finite(ll))) return(1e100)
  -sum(ll)
}

sim_loglikM_truncated_GB2_cpp <- function(Claim, Deduct,
                                           sim_covapos, parms) {
  k <- ncol(sim_covapos)
  eta <- drop(sim_covapos %*% parms[seq_len(k)])
  a <- parms[k+1L]
  logsurv <- -a * log1p(Deduct * exp(-eta))
  ll <- lomax_ll(Claim, eta, a) - logsurv
  if (any(!is.finite(ll))) return(1e100)
  -sum(ll)
}


# 1. Simulation Parameters
n <- 5000
t <- 500

# Covariates: X = (1, X1, X2)
# X1 ~ N(0,1), X2 ~ Bernoulli(0.6)

# Frequency: NB
r <- 1
beta_N <- c(-0.5, 0.5, 0.2)

# Deductible: proportional-odds logistic model
ded_levels <- c(500, 1000, 2500, 5000, 10000)
kappa <- c(0.1, 1.0, 1.8, 2.6)
beta_D <- c(-0.5, -0.5)

# Severity: Pareto (Lomax), GB2 special case
beta_Y <- c(8.5, 0.2, 0.5)
sigma <- 1
alpha1 <- 1
alpha2 <- 3

# Copula correlations for Case 4 (D--frequency, D--severity)
rho_N <- -0.6
rho_Y <- -0.4


# 2. Storage for subsequent estimation
D_table <- matrix(NA_real_, nrow = 6, ncol = t)

FEM1_table  <- matrix(NA_real_, nrow = 4, ncol = t)
FEM2G_table <- matrix(NA_real_, nrow = 4, ncol = t)
FEM2T_table <- matrix(NA_real_, nrow = 4, ncol = t)

SEM1_table <- matrix(NA_real_, nrow = 4, ncol = t)
SEM2_table <- matrix(NA_real_, nrow = 4, ncol = t)
SEM3_table <- matrix(NA_real_, nrow = 4, ncol = t)

FDM1_table <- matrix(NA_real_, nrow = 5, ncol = t)
FDM2_table <- matrix(NA_real_, nrow = 5, ncol = t)

SDM1_table <- matrix(NA_real_, nrow = 5, ncol = t)
SDM2_table <- matrix(NA_real_, nrow = 5, ncol = t)
SDM3_table <- matrix(NA_real_, nrow = 5, ncol = t)

# [ADDED] Copula convergence and optional Godambe outputs
copula_convergence <- matrix(NA_integer_, nrow=t, ncol=6,
  dimnames=list(NULL,c("D","SDM1","SDM2","SDM3","FDM1","FDM2")))
sandwich_rows <- list()

sim_dat.pol.list <- vector("list", t)
sim_dat.clm.list <- vector("list", t)
sim_dat.clmpos.list <- vector("list", t)

sim_px_G.list <- vector("list", t)
sim_px_T.list <- vector("list", t)

# [CHANGED] Omitted unused component-level MAE tables; S_metrics is
# created below for six aggregate insured loss model combinations.

# Same random seed as the uploaded Case 1 reference script.
set.seed(123)
# Clip latent CDF values only to avoid qnorm(0/1) and Lomax overflow.
clip_u <- function(u) pmin(1 - 1e-13, pmax(1e-13, u))

# 3. Data generation
for (s in 1:t) {
  
  # Same seed convention as professor's simulation
  
  # --------------------------------------------------------
  # Step 1: Generate covariates
  # --------------------------------------------------------
  X <- cbind(
    int = 1,
    X1 = rnorm(n, mean = 0, sd = 1),
    X2 = rbinom(n, size = 1, prob = 0.6)
  )
  
  # --------------------------------------------------------
  # Step 2: Generate discrete deductible
  # --------------------------------------------------------
  
  # Linear predictor without intercept
  eta_D <- drop(X[, 2:3, drop = FALSE] %*% beta_D)
  
  # Cumulative probabilities for the first 4 levels
  F_D <- sapply(kappa, function(k) plogis(k + eta_D))
  
  # F_D is n x 4:
  # P(D <= 500), P(D <= 1000), ..., P(D <= 5000)
  
  # Latent uniform for deductible selection
  U_D <- runif(n)
  
  # Inverse CDF sampling: deductible level 1, ..., 5
  D_level <- 1L + rowSums(
    U_D > F_D
  )
  
  D <- ded_levels[D_level]
  
  # --------------------------------------------------------
  # Step 3: Generate ground-up frequency N (legacy R: N_star)
  # Conditional Gaussian copula for a DISCRETE deductible.
  # The deductible atom is (F_D(d- | x), F_D(d | x)].
  # An independent latent U_D draw is used for N, conditional on D.
  # --------------------------------------------------------
  mu_N <- drop(exp(X %*% beta_N))

  atom_cdf <- cbind(0, F_D, 1)
  atom_lo <- atom_cdf[cbind(seq_len(n), D_level)]
  atom_hi <- atom_cdf[cbind(seq_len(n), D_level + 1L)]
  stopifnot(all(atom_hi > atom_lo), all(is.finite(atom_hi)))

  if (abs(rho_N) < 1e-14) {
    N_star <- rnbinom(n, size = r, mu = mu_N)
  } else {
    # Drawing within D's atom and THEN correlating latent normals
    # gives the correct D--N Gaussian copula conditional on the D atom.
    Ud_N <- atom_lo + (atom_hi - atom_lo) * runif(n)
    Zn <- rho_N * qnorm(clip_u(Ud_N)) +
      sqrt(1 - rho_N^2) * rnorm(n)
    N_star <- qnbinom(clip_u(pnorm(Zn)), size = r, mu = mu_N)
  }
  N_star <- as.integer(N_star)

  # --------------------------------------------------------
  # Step 4: Generate ground-up severity Y (legacy R: Y_star)
  # Resample an INDEPENDENT deductible-atom latent for each claim,
  # including independence from the frequency latent above.
  # This preserves N independent of individual severities given (D,X).
  # --------------------------------------------------------
  lambda_Y <- drop(exp(X %*% beta_Y))
  Y_star <- vector("list", n)

  for (i in seq_len(n)) {
    ni <- N_star[i]
    if (ni == 0L) {
      Y_star[[i]] <- numeric(0)
      next
    }
    if (abs(rho_Y) < 1e-14) {
      uy <- runif(ni)
    } else {
      Ud_Y <- atom_lo[i] + (atom_hi[i] - atom_lo[i]) * runif(ni)
      Zy <- rho_Y * qnorm(clip_u(Ud_Y)) +
        sqrt(1 - rho_Y^2) * rnorm(ni)
      uy <- clip_u(pnorm(Zy))
    }
    # F_Y^{-1}(u) for Lomax scale=lambda_Y and shape=alpha2.
    Y_star[[i]] <- lambda_Y[i] * expm1(-log1p(-uy) / alpha2)
  }
  Y_star_sum <- vapply(Y_star, sum, numeric(1))

  # --------------------------------------------------------
  # Step 5: Calculate number of positive payments N* (legacy R variable: N)
  # --------------------------------------------------------
  N <- integer(n)
  
  for (i in 1:n) {
    N[i] <- sum(Y_star[[i]] > D[i])
  }
  
  # --------------------------------------------------------
  # Step 6: Policy-level dataset
  # --------------------------------------------------------
  sim_dat.pol <- data.frame(
    policy_id = 1:n,
    X1 = X[, 2],
    X2 = X[, 3],
    Deduct = D,
    DeductLevel = D_level,
    Freq = N_star,
    pFreq = N,
    Sum.claim = Y_star_sum
  )
  
  # --------------------------------------------------------
  # Step 7: Claim-level dataset
  # --------------------------------------------------------
  
  policy_id <- rep(seq_len(n), times = N_star)
  claim_no <- unlist(
    lapply(N_star, seq_len),
    use.names = FALSE
  )
  
  sim_dat.clm <- data.frame(
    policy_id = policy_id,
    claim_no = claim_no,
    Claim = unlist(Y_star, use.names = FALSE),
    X1 = X[policy_id, 2],
    X2 = X[policy_id, 3],
    Deduct = D[policy_id],
    DeductLevel = D_level[policy_id],
    Freq = N_star[policy_id]
  )
  
  # --------------------------------------------------------
  # Step 8: Payment after deductible
  # --------------------------------------------------------
  sim_dat.clm$LossAfterDeductible <- pmax(
    sim_dat.clm$Claim - sim_dat.clm$Deduct,
    0
  )
  
  sim_dat.clmpos <- subset(
    sim_dat.clm,
    LossAfterDeductible > 0
  )
  
  # --------------------------------------------------------
  # Step 9: Total payment by policy
  # --------------------------------------------------------
  
  # [CHANGED 2] Empty claim datasets are possible in small pilot runs.
  sim_dat.pol$Sum.LossAfterDeductible <- 0
  if (nrow(sim_dat.clm) > 0L) {
    payment_sum <- rowsum(
      sim_dat.clm$LossAfterDeductible,
      group = sim_dat.clm$policy_id,
      reorder = TRUE
    )
    idx <- as.integer(rownames(payment_sum))
    sim_dat.pol$Sum.LossAfterDeductible[idx] <- as.numeric(payment_sum[, 1])
  }
  
  # --------------------------------------------------------
  # Step 10: Save the generated dataset in memory
  # --------------------------------------------------------
  sim_dat.pol.list[[s]] <- sim_dat.pol
  sim_dat.clm.list[[s]] <- sim_dat.clm
  sim_dat.clmpos.list[[s]] <- sim_dat.clmpos
  
  # DGP diagnostics are summarized across policies after estimation.
  # E.g. E[N|D=d,X] is altered by rho_N; E[Y|D=d,X] by rho_Y.

  # Progress
  if (s %% 10 == 0 || s == 1) {
    cat(sprintf(
      "Case 4: %d/%d | claims: %d | payments: %d\n",
      s, t,
      nrow(sim_dat.clm),
      nrow(sim_dat.clmpos)
    ))
  }
}




###############################################################################
# [ADDED 3] DISCRETE DEDUCTIBLE GAUSSIAN COPULA ROUTINES
# These functions operate on a thin adapter of YOUR original datasets.
# D atom is (F_D(d-), F_D(d)] rather than continuous F_D(d).
###############################################################################

DED_LEVELS <- ded_levels
KAPPA <- kappa
BETA_D <- beta_D
BETA_N <- beta_N
BETA_Y <- beta_Y
TRUE_R <- r
TRUE_ALPHA <- alpha2
RHO_N <- rho_N
RHO_Y <- rho_Y
Q_ATOM <- 24L        # Gauss-Legendre nodes inside each deductible atom
Q_SEV <- 80L         # Nodes for E[(Y-D)+ | D, X]
TAIL_EPS <- 1e-9     # Bound on unconditional omitted count tail
MAX_COUNT <- 700L   # If reached, mark likelihood evaluation invalid
COPULA_MAXIT <- 250L
# ---------------- USER CONTROLS: GODAMBE / COVERAGE --------------------
# Set TRUE to compute selected Godambe, stage-wise robust, and naive SEs.
# FALSE skips ALL standard-error calculations, but retains every model fit.
RUN_SANDWICH <- TRUE
# Which replications to assess? For a trial use 1:min(10,t), for coverage
# with 500 replications use seq_len(t) AFTER checking runtime.
SANDWICH_REPS <- seq_len(t)
# Any subset of: D, SEM1, SEM2, SEM3, FEM1, FEM2G, FEM2,
# SDM1, SDM2, SDM3, FDM1, FDM2. Prerequisite stages are automatic.
# Suggested Case 4 subset (customize to your manuscript or appendix):
SANDWICH_TARGETS <- c("SDM3", "FDM2")
# ----------------------------------------------------------------------
SANDWICH_TIMING <- data.frame(rep=integer(), model=character(), elapsed_seconds=double(), succeeded=logical())

# Input data adapter: preserves your names, adds paper's N,N*,Y,Y* aliases.
as_model_data <- function(pol, clm, pos) {
  pol$int <- 1
  pol$N <- pol$Freq
  pol$Nstar <- pol$pFreq
  pol$S <- pol$Sum.LossAfterDeductible
  clm$Y <- clm$Claim
  clm$Ystar <- clm$LossAfterDeductible
  pos$Y <- pos$Claim
  pos$Ystar <- pos$LossAfterDeductible
  list(pol = pol, clm = clm, pos = pos)
}

pn <- c("beta0", "beta1", "beta2")
pclip <- function(x) pmin(1 - 1e-13, pmax(1e-13, x))
safe_log <- function(x) log(pmax(x, 1e-300))

# Gauss-Legendre quadrature from the symmetric Jacobi matrix (base R only).
gl_nodes <- function(n) {
  j <- seq_len(n - 1L)
  a <- j / sqrt(4 * j^2 - 1)
  A <- matrix(0, n, n)
  A[cbind(j, j + 1L)] <- a
  A[cbind(j + 1L, j)] <- a
  eig <- eigen(A, symmetric = TRUE)
  ord <- order(eig$values)
  list(x = eig$values[ord], w = 2 * eig$vectors[1L, ord]^2)
}
GLA <- gl_nodes(Q_ATOM)
GLS <- gl_nodes(Q_SEV)

# All optimization is on unconstrained coordinates. The reported parameters
# are the original beta, alpha, r, rho, kappa scales.
d_unpack <- function(z) {
  k <- c(z[1L], z[1L] + cumsum(exp(z[2:4])))
  c(k, betaD1 = z[5L], betaD2 = z[6L])
}
d_pack <- function(k, b) c(k[1L], log(diff(k)), b)
s_unpack <- function(z, dep) {
  a <- c(setNames(z[1:3], pn), alpha = 1 + exp(z[4L]))
  if (dep) a <- c(a, rho = pmax(-0.995, pmin(0.995, tanh(z[5L]))))
  a
}
f_unpack <- function(z, dep) {
  a <- c(setNames(z[1:3], pn), r = exp(z[4L]))
  if (dep) a <- c(a, rho = pmax(-0.995, pmin(0.995, tanh(z[5L]))))
  a
}
true_par <- function(name) {
  if (name == "D") return(c(setNames(KAPPA, paste0("kappa", 1:4)), betaD1 = BETA_D[1], betaD2 = BETA_D[2]))
  if (startsWith(name, "S")) return(c(setNames(BETA_Y, pn), alpha = TRUE_ALPHA,
                                      if (startsWith(name, "SDM")) c(rho = RHO_Y)))
  c(setNames(BETA_N, pn), r = TRUE_R,
    if (startsWith(name, "FDM")) c(rho = RHO_N))
}

# Matrix columns: 0, F_D(d1), F_D(d2), F_D(d3), F_D(d4), 1
D_cdf_mat <- function(X, pars) {
  eta <- drop(X[, 2:3, drop = FALSE] %*% pars[5:6])
  cbind(0, sapply(pars[1:4], function(k) plogis(k + eta)), 1)
}
atom_bounds <- function(pol, Dtheta = NULL) {
  if (is.null(Dtheta)) return(list(lo = rep(0, nrow(pol)), hi = rep(1, nrow(pol))))
  Fd <- D_cdf_mat(as.matrix(pol[, c("int", "X1", "X2")]), d_unpack(Dtheta))
  ii <- seq_len(nrow(pol))
  list(lo = Fd[cbind(ii, pol$DeductLevel)], hi = Fd[cbind(ii, pol$DeductLevel + 1L)])
}

# Normal interval probabilities: avoid catastrophic cancellation near 1.
norm_between <- function(a, b) {
  left <- pnorm(b) - pnorm(a)
  right <- pnorm(a, lower.tail = FALSE) - pnorm(b, lower.tail = FALSE)
  pmax(0, ifelse(a > 0, right, left))
}

# Integral of a conditional Gaussian-normal CDF over uniform u_D in its atom.
# F(Z_other <= z | D=d, X=x) = integral Phi((z-rho*Phi^-1(u))/s) du / (hi-lo).
atom_normal_cdf <- function(z, lo, hi, rho) {
  if (abs(rho) < 1e-14) return(pnorm(z))
  n <- length(z)
  ans <- numeric(n)
  for (j in seq_along(GLA$x)) {
    u <- lo + (hi - lo) * (GLA$x[j] + 1) / 2
    ans <- ans + (GLA$w[j] / 2) * pnorm((z - rho * qnorm(pclip(u))) / sqrt(1 - rho^2))
  }
  ans
}

# Gaussian copula weight of a conditional density of Y* given deductible atom.
sev_weight <- function(y, scale, alpha, rho, lo, hi) {
  if (abs(rho) < 1e-14) return(rep(1, length(y)))
  sy <- (1 + y / scale)^(-alpha)
  zy <- qnorm(pclip(1 - sy))
  a <- (qnorm(pclip(lo)) - rho * zy) / sqrt(1 - rho^2)
  b <- (qnorm(pclip(hi)) - rho * zy) / sqrt(1 - rho^2)
  norm_between(a, b) / pmax(hi - lo, 1e-15)
}

# P(Y* > d | D,X), evaluated directly by conditioning on U_D.
retention <- function(Sp, pol, atom) {
  X <- as.matrix(pol[, c("int", "X1", "X2")])
  lam <- drop(exp(X %*% Sp[1:3])); alpha <- unname(Sp[4])
  rho <- if (length(Sp) >= 5) unname(Sp[5]) else 0
  survival <- (1 + pol$Deduct / lam)^(-alpha)
  if (abs(rho) < 1e-14) return(pmax(1e-12, survival))
  zd <- qnorm(pclip(1 - survival))
  ans <- numeric(nrow(pol))
  for (j in seq_along(GLA$x)) {
    u <- atom$lo + (atom$hi - atom$lo) * (GLA$x[j] + 1) / 2
    ans <- ans + GLA$w[j] / 2 * pnorm((rho * qnorm(pclip(u)) - zd) / sqrt(1 - rho^2))
  }
  pmin(1 - 1e-12, pmax(1e-12, ans))
}

# P(N*=k | D,X) using the appropriate Gaussian-copula rectangle.
count_pmf <- function(k, mu, r, rho, atom) {
  if (abs(rho) < 1e-14) return(dnbinom(k, size = r, mu = mu))
  zhi <- qnorm(pnbinom(k, size = r, mu = mu))
  zlo <- qnorm(pnbinom(k - 1L, size = r, mu = mu))
  ans <- numeric(length(mu))
  for (j in seq_along(GLA$x)) {
    u <- atom$lo + (atom$hi - atom$lo) * (GLA$x[j] + 1) / 2
    z <- rho * qnorm(pclip(u))
    s <- sqrt(1 - rho^2)
    ans <- ans + GLA$w[j] / 2 * norm_between((zlo - z) / s, (zhi - z) / s)
  }
  pmax(ans, 0)
}

# Exact thinning convolution given the discrete-deductible Gaussian copula:
# P(Nstar=m|D,X) = sum_{k>=m} Binomial(m;k,p_ret) P(N=k|D,X).
# The finite sum is controlled by TAIL_EPS. NB offset is NOT exact for FDM2.
count_thinned <- function(m, mu, r, rho, atom, pret) {
  if (abs(rho) < 1e-14) return(dnbinom(m, size = r, mu = mu * pret))
  minwidth <- min(atom$hi - atom$lo)
  K <- max(max(m), qnbinom(1 - TAIL_EPS * minwidth, size = r, mu = max(mu)))
  if (!is.finite(K) || K > MAX_COUNT) return(rep(NA_real_, length(mu)))
  ans <- numeric(length(mu))
  for (k in 0:as.integer(K)) {
    valid <- which(m <= k)
    if (!length(valid)) next
    av <- list(lo = atom$lo[valid], hi = atom$hi[valid])
    pk <- count_pmf(rep.int(k, length(valid)), mu[valid], r, rho, av)
    ans[valid] <- ans[valid] + pk * dbinom(m[valid], k, pret[valid])
  }
  pmax(ans, 1e-300)
}

# Convert claim-level contributions into independent policy-cluster contributions.
by_policy <- function(v, id, n) {
  ans <- numeric(n)
  if (!length(v)) return(ans)
  g <- rowsum(v, id, reorder = FALSE)
  ans[as.integer(rownames(g))] <- g[, 1L]
  ans
}

# Stage graph: all inference is sequential; predecessors are estimated, not fixed at truth.
DEPS <- list(D = character(), SEM1 = character(), SEM2 = character(), SEM3 = character(),
             FEM1 = character(), FEM2G = "SEM1", FEM2 = "SEM3",
             SDM1 = "D", SDM2 = "D", SDM3 = "D", FDM1 = "D",
             FDM2 = c("D", "SDM3"))
ORDER <- c("D", "SEM1", "SEM2", "SEM3", "FEM1", "FEM2G", "FEM2",
           "SDM1", "SDM2", "SDM3", "FDM1", "FDM2")

# fitz is a named list of unconstrained stage coordinates.
ll_policy <- function(name, z, fitz, dat) {
  pol <- dat$pol; clm <- dat$clm; pos <- dat$pos; n <- nrow(pol)
  Xp <- as.matrix(pol[, c("int", "X1", "X2")])
  if (name == "D") {
    Fd <- D_cdf_mat(Xp, d_unpack(z))
    lev <- pol$DeductLevel
    prob <- Fd[cbind(seq_len(n), lev + 1L)] - Fd[cbind(seq_len(n), lev)]
    return(safe_log(prob))
  }
  dep <- name %in% c("SDM1", "SDM2", "SDM3", "FDM1", "FDM2")
  atom <- atom_bounds(pol, if (dep) fitz$D else NULL)

  if (startsWith(name, "S")) {
    Sp <- s_unpack(z, dep)
    cc <- if (name %in% c("SEM3", "SDM3")) pos else clm
    if (!nrow(cc)) return(rep(0, n))
    ii <- cc$policy_id
    yy <- if (name %in% c("SEM1", "SDM1")) cc$Y else cc$Ystar + cc$Deduct
    lam <- drop(exp(Xp[ii, , drop = FALSE] %*% Sp[1:3]))
    a <- unname(Sp[4]); rho <- if (dep) unname(Sp[5]) else 0
    lden <- log(a) - log(lam) - (a + 1) * log1p(yy / lam)
    if (dep) {
      w <- sev_weight(yy, lam, a, rho, atom$lo[ii], atom$hi[ii])
      lden <- lden + safe_log(w)
    }
    if (name %in% c("SEM2", "SDM2", "SEM3", "SDM3")) {
      p <- retention(Sp, pol, atom)
      if (name %in% c("SEM3", "SDM3")) {
        lden <- lden - safe_log(p[ii])
      } else {
        lden[cc$Ystar == 0] <- log1p(-p[ii[cc$Ystar == 0]])
      }
    }
    return(by_policy(lden, ii, n))
  }

  Fp <- f_unpack(z, dep)
  mu <- drop(exp(Xp %*% Fp[1:3])); r <- unname(Fp[4])
  rho <- if (dep) unname(Fp[5]) else 0
  if (name %in% c("FEM1", "FDM1")) {
    return(safe_log(count_pmf(pol$N, mu, r, rho, atom)))
  }
  spec <- if (name == "FEM2") "SEM3" else if (name == "FEM2G") "SEM1" else "SDM3"
  Sp <- s_unpack(fitz[[spec]], spec == "SDM3")
  atS <- atom_bounds(pol, if (spec == "SDM3") fitz$D else NULL)
  pret <- retention(Sp, pol, atS)
  safe_log(count_thinned(pol$Nstar, mu, r, rho, atom, pret))
}

stage_nll <- function(name, z, fitz, dat) {
  ll <- ll_policy(name, z, fitz, dat)
  if (any(!is.finite(ll))) return(1e200)
  -sum(ll)
}

# Conditional expected number of loss events, E[N|D,X], using tail-sum identity.
freq_mean <- function(Fp, pol, atom) {
  X <- as.matrix(pol[, c("int", "X1", "X2")])
  mu <- drop(exp(X %*% Fp[1:3])); r <- unname(Fp[4])
  rho <- if (length(Fp) > 4) unname(Fp[5]) else 0
  if (abs(rho) < 1e-14) return(mu)
  K <- qnbinom(1 - TAIL_EPS * min(atom$hi - atom$lo), size = r, mu = max(mu))
  if (!is.finite(K) || K > MAX_COUNT) stop("Increase MAX_COUNT for frequency mean")
  EN <- numeric(nrow(pol))
  for (k in 0:as.integer(K)) {
    q <- qnorm(pnbinom(k, size = r, mu = mu))
    EN <- EN + pmax(0, 1 - atom_normal_cdf(q, atom$lo, atom$hi, rho))
  }
  EN
}

# Conditional insured severity mean E[(Y-D)+|D,X].
# For rho=0, use closed form. For a Gaussian copula with the deductible,
# evaluate only the insured-mean integral (no unneeded ground-up/positive means).
sev_paid_mean <- function(Sp, pol, atom) {
  X <- as.matrix(pol[, c("int", "X1", "X2")])
  lam <- drop(exp(X %*% Sp[1:3])); alpha <- unname(Sp[4])
  rho <- if (length(Sp) > 4) unname(Sp[5]) else 0
  d <- pol$Deduct
  if (abs(rho) < 1e-14) {
    return(lam / (alpha - 1) * (1 + d / lam)^(1 - alpha))
  }
  zlower <- qnorm(pclip(1 - (1 + d / lam)^(-alpha)))
  # Normal tail beyond 10 is negligible for this design; test Q_SEV convergence.
  zupper <- pmax(10, zlower + 0.02)
  width <- zupper - zlower
  val <- numeric(length(lam))
  a <- qnorm(pclip(atom$lo)); b <- qnorm(pclip(atom$hi))
  ss <- sqrt(1 - rho^2)
  for (j in seq_along(GLS$x)) {
    z <- zlower + width * (GLS$x[j] + 1) / 2
    tail <- pmax(pnorm(z, lower.tail = FALSE), 1e-290)
    yy <- lam * expm1(-log(tail) / alpha)
    w <- norm_between((a - rho * z) / ss, (b - rho * z) / ss) /
      pmax(atom$hi - atom$lo, 1e-15)
    val <- val + GLS$w[j] / 2 * width * pmax(0, yy - d) * dnorm(z) * w
  }
  val
}

# Base-R numeric differentiation of per-policy stage scores.
policy_scores <- function(name, z, fitz, dat) {
  h <- 1e-5 * (1 + abs(z))
  out <- matrix(NA_real_, nrow(dat$pol), length(z))
  for (j in seq_along(z)) {
    zp <- zm <- z; zp[j] <- zp[j] + h[j]; zm[j] <- zm[j] - h[j]
    out[, j] <- (ll_policy(name, zp, fitz, dat) -
                  ll_policy(name, zm, fitz, dat)) / (2 * h[j])
  }
  out
}

# Mixed Hessian of stage's negative log-likelihood with respect to
# current theta (rows) and estimated earlier theta (columns).
cross_hessian <- function(name, own, earlier_name, fitz, dat) {
  base <- fitz[[earlier_name]]
  h1 <- 1e-4 * (1 + abs(own)); h2 <- 1e-4 * (1 + abs(base))
  ans <- matrix(NA_real_, length(own), length(base))
  for (i in seq_along(own)) for (j in seq_along(base)) {
    pp <- pm <- mp <- mm <- own
    pp[i] <- pp[i] + h1[i]; pm[i] <- pm[i] + h1[i]
    mp[i] <- mp[i] - h1[i]; mm[i] <- mm[i] - h1[i]
    fplus <- fminus <- fitz
    fplus[[earlier_name]][j] <- base[j] + h2[j]
    fminus[[earlier_name]][j] <- base[j] - h2[j]
    ans[i, j] <- (stage_nll(name, pp, fplus, dat) -
                  stage_nll(name, pm, fminus, dat) -
                  stage_nll(name, mp, fplus, dat) +
                  stage_nll(name, mm, fminus, dat)) / (4 * h1[i] * h2[j])
  }
  ans
}

# Generalized sequential Godambe sandwich for every estimable model.
# The stage-wise robust sandwich intentionally ignores upstream estimation.
# Per-replication cache avoids repeating identical derivative calculations
# when requested models share a preceding stage (e.g., SDM3 and FDM2).
GODAMBE_CHAINS <- list(
  D = c("D"),
  SEM1 = c("SEM1"), SEM2 = c("SEM2"), SEM3 = c("SEM3"),
  FEM1 = c("FEM1"),
  FEM2G = c("SEM1", "FEM2G"),
  FEM2 = c("SEM3", "FEM2"),
  SDM1 = c("D", "SDM1"),
  SDM2 = c("D", "SDM2"),
  SDM3 = c("D", "SDM3"),
  FDM1 = c("D", "FDM1"),
  FDM2 = c("D", "SDM3", "FDM2")
)

cache_compute <- function(cache, key, fun) {
  if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  value <- fun()
  if (!is.null(cache)) assign(key, value, envir = cache)
  value
}

# Jacobian of the map from optimizer coordinates to reported parameters.
# Crucially, ordered D cutpoints kappa_2..kappa_4 share previous increments:
# a diagonal-only delta method would give incorrect cutpoint standard errors.
native_jacobian <- function(name, z) {
  q <- length(z)
  A <- diag(q)
  if (name == "D") {
    A <- matrix(0, q, q)
    A[1:4, 1L] <- 1
    for (j in 2:4) A[j:4, j] <- exp(z[j])
    A[5L, 5L] <- 1
    A[6L, 6L] <- 1
  } else {
    A[4L, 4L] <- if (startsWith(name, "S")) exp(z[4L]) else exp(z[4L])
    if (q == 5L) A[5L, 5L] <- 1 - tanh(z[5L])^2
  }
  A
}

reported_params <- function(name, z) {
  if (name == "D") {
    a <- d_unpack(z)
    names(a) <- c(paste0("kappa", 1:4), "betaD1", "betaD2")
    return(a)
  }
  if (startsWith(name, "S")) return(s_unpack(z, startsWith(name, "SDM")))
  f_unpack(z, startsWith(name, "FDM"))
}

godambe <- function(target, fits, dat, cache = NULL) {
  chain <- GODAMBE_CHAINS[[target]]
  if (is.null(chain)) stop("Unrecognized Godambe target: ", target)
  fitz <- lapply(fits, `[[`, "z")
  if (!all(chain %in% names(fitz))) stop("Missing fitted stages in chain: ", target)
  sizes <- vapply(chain, function(nm) length(fitz[[nm]]), integer(1))
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  idx <- Map(seq.int, starts, ends)
  q <- sum(sizes); n_pol <- nrow(dat$pol)
  H <- matrix(0, q, q)
  U <- matrix(0, n_pol, q)

  for (j in seq_along(chain)) {
    name <- chain[j]
    z <- fitz[[name]]
    U[, idx[[j]]] <- cache_compute(cache, paste0("score_", name), function() {
      policy_scores(name, z, fitz, dat)
    })
    H[idx[[j]], idx[[j]]] <- cache_compute(cache, paste0("hessian_", name), function() {
      optimHess(z, function(v) stage_nll(name, v, fitz, dat))
    })
    if (j > 1L) for (k in seq_len(j - 1L)) {
      before <- chain[k]
      if (before %in% DEPS[[name]]) {
        H[idx[[j]], idx[[k]]] <- cache_compute(
          cache, paste0("cross_", name, "_", before), function() {
            cross_hessian(name, z, before, fitz, dat)
          }
        )
      }
    }
  }
  if (any(!is.finite(H)) || any(!is.finite(U)))
    stop("Non-finite score or Hessian for ", target)

  J <- crossprod(U)
  inv <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(inv)) stop("Godambe Hessian singular for ", target)
  V <- inv %*% J %*% t(inv)
  last <- idx[[length(chain)]]
  hlast <- H[last, last, drop = FALSE]
  jlast <- J[last, last, drop = FALSE]
  Vinv <- tryCatch(solve(hlast), error = function(e) NULL)
  if (is.null(Vinv)) stop("Stage Hessian singular for ", target)
  Vstage <- Vinv %*% jlast %*% t(Vinv)

  z <- fitz[[target]]
  estimate <- reported_params(target, z)
  A <- native_jacobian(target, z)
  native_se <- function(M) {
    W <- A %*% M %*% t(A)
    v <- diag(W)
    if (any(!is.finite(v)) || any(v < -1e-8))
      stop("Negative or non-finite variance for ", target)
    sqrt(pmax(0, v))
  }
  ans <- data.frame(
    model = target,
    parameter = names(estimate),
    estimate = unname(estimate),
    true = unname(true_par(target)),
    naive_se = native_se(Vinv),
    stage_robust_se = native_se(Vstage),
    godambe_se = native_se(V[last, last, drop = FALSE])
  )
  if (any(!is.finite(as.matrix(ans[, c("estimate", "true", "naive_se",
                                        "stage_robust_se", "godambe_se")]))))
    stop("Invalid Godambe output for ", target)
  ans
}

# [ADDED 4] Optimize one conditional log-likelihood using the preceding stages.
# z = (beta_0,beta_1,beta_2, log(alpha-1) or log(r), atanh(rho)).
fit_copula_stage <- function(stage, z0, fitz, dat) {
  obj <- function(z) stage_nll(stage, z, fitz, dat)
  fit <- tryCatch(optim(z0, obj, method = "BFGS",
                        control = list(maxit=COPULA_MAXIT, reltol=2e-8)),
                  error=function(e) NULL)
  if (is.null(fit) || !is.finite(fit$value) || fit$convergence != 0L) {
    alt <- tryCatch(nlminb(z0, obj,
                           control=list(iter.max=COPULA_MAXIT,
                                        eval.max=3L*COPULA_MAXIT)),
                    error=function(e) NULL)
    if (!is.null(alt) && is.finite(alt$objective) &&
        (is.null(fit) || !is.finite(fit$value) || alt$objective < fit$value)) {
      fit <- list(par=alt$par, value=alt$objective, convergence=alt$convergence)
    }
  }
  if (is.null(fit) || !is.finite(fit$value))
    stop("Copula stage failed: ",stage)
  if (fit$convergence != 0)
    warning(sprintf("%s optimization convergence=%s",stage,fit$convergence))
  fit
}

# H=total derivative of scores; J=sum of policy outer products.
# The total (not average) versions already include 1/n; do not rescale.
if (!is.logical(RUN_SANDWICH) || length(RUN_SANDWICH) != 1L)
  stop("RUN_SANDWICH must be TRUE or FALSE")
if (any(!SANDWICH_TARGETS %in% names(GODAMBE_CHAINS)))
  stop("Unknown SANDWICH_TARGETS: ",
       paste(setdiff(SANDWICH_TARGETS, names(GODAMBE_CHAINS)), collapse = ", "))
if (any(!SANDWICH_REPS %in% seq_len(t)))
  stop("SANDWICH_REPS must be between 1 and t")

for (s in 1:t) {
  if (s %% max(1L, floor(t/10)) == 0) {
    cat(sprintf("%.0f%% completed (%d/%d)\n",
                100*s/t, s, t))
  }
  
  sim_dat.clm <- sim_dat.clm.list[[s]]
  sim_dat.clmpos <- sim_dat.clmpos.list[[s]]
  sim_dat.pol <- sim_dat.pol.list[[s]]
  
  sim_covaS <- model.matrix(~ X1 + X2, data=sim_dat.clm)
  sim_covapos <- model.matrix(~ X1 + X2, data=sim_dat.clmpos)
  
  
  ## Severity
  
  
  # loglikelihood GB2
  sim_loglikM_ground_GB2 <- function(parms) {
    sim_loglikM_ground_GB2_cpp(sim_dat.clm$Claim, sim_covaS, parms)
  }
  
  sim_loglikM_censored_GB2 <- function(parms) {
    sim_loglikM_censored_GB2_cpp(sim_dat.clm$LossAfterDeductible, sim_dat.clm$Deduct, sim_covaS, parms)
  }
  
  sim_loglikM_truncated_GB2 <- function(parms) {
    sim_loglikM_truncated_GB2_cpp(sim_dat.clmpos$Claim, sim_dat.clmpos$Deduct, sim_covapos, parms)
  }
  
  # SEM1: Full observations
  sim_fit_ground <- lm(log(Claim) ~ X1 + X2,data = sim_dat.clm)
  
  # SEM2: Positive observed payments for initialization
  sim_fit_censored <- lm(log(LossAfterDeductible) ~ X1 + X2,data = sim_dat.clmpos)
  
  # SEM3: Truncated observations
  sim_fit_truncated <- lm(log(LossAfterDeductible) ~ X1 + X2,data = sim_dat.clmpos)
  
  sim_init_ground_GB2 <- c(coef(sim_fit_ground), 2)
  
  sim_init_censored_GB2 <- c(coef(sim_fit_censored), 2)
  
  sim_init_truncated_GB2 <- c(coef(sim_fit_truncated), 2)
  
  # model fitting
  #system.time(
  sim_zop_ground_GB2      <- nlminb(sim_init_ground_GB2, sim_loglikM_ground_GB2, lower =c(rep(-Inf,ncol(sim_covaS))  ,1.001))
  #)
  
  #system.time(
  sim_zop_censored_GB2      <- nlminb(sim_init_censored_GB2, sim_loglikM_censored_GB2, lower =c(rep(-Inf,ncol(sim_covaS)) ,1.001))
  #)
  
  #system.time(
  sim_zop_truncated_GB2    <- nlminb(sim_init_truncated_GB2, sim_loglikM_truncated_GB2 , lower =c(rep(-Inf,ncol(sim_covapos)),1.001))
  #)
  
  # Estimates from Severity
  sim_betaest_ground   <- sim_zop_ground_GB2$par[1:ncol(sim_covaS)]
  sim_sigmaest_ground  <- 1
  sim_alpha1est_ground <- 1
  sim_alpha2est_ground <- sim_zop_ground_GB2$par[ncol(sim_covaS)+1]
  
  SEM1_table[,s] <- c(sim_betaest_ground, sim_alpha2est_ground)
  
  
  sim_betaest_censored   <- sim_zop_censored_GB2$par[1:ncol(sim_covaS)]
  sim_sigmaest_censored  <- 1
  sim_alpha1est_censored <- 1
  sim_alpha2est_censored <- sim_zop_censored_GB2$par[ncol(sim_covaS)+1]
  
  SEM2_table[,s] <- c(sim_betaest_censored, sim_alpha2est_censored)
  
  sim_betaest_truncated   <- sim_zop_truncated_GB2$par[1:ncol(sim_covapos)]
  sim_sigmaest_truncated  <- 1
  sim_alpha1est_truncated <- 1
  sim_alpha2est_truncated <- sim_zop_truncated_GB2$par[ncol(sim_covapos)+1]
  
  SEM3_table[,s] <- c(sim_betaest_truncated, sim_alpha2est_truncated)
  
  ## Frequency test
  
  sim_cov <- model.matrix(~ X1 + X2, data=sim_dat.pol)
  
  # weights for adjustment for truncated frequency
  sim_px_G.list[[s]]    <- 1- pGB2_cpp(sim_dat.pol$Deduct, sim_cov %*% sim_betaest_ground,
                                       sim_sigmaest_ground, sim_alpha1est_ground,sim_alpha2est_ground)
  
  sim_px_T.list[[s]]    <- 1- pGB2_cpp(sim_dat.pol$Deduct, sim_cov %*% sim_betaest_truncated,
                                       sim_sigmaest_truncated,
                                       sim_alpha1est_truncated,
                                       sim_alpha2est_truncated)
  sim_px_G <- sim_px_G.list[[s]]
  sim_px_T <- sim_px_T.list[[s]]
  
  sim_nop_censored   <- glm.nb(Freq~X1+X2, data=sim_dat.pol)
  sim_nop_truncatedG <- glm.nb(pFreq~X1+X2+offset(log(sim_px_G)), data=sim_dat.pol)
  sim_nop_truncatedT <- glm.nb(pFreq~X1+X2+offset(log(sim_px_T)), data=sim_dat.pol)
  
  FEM1_table[,s] <- c(coef(sim_nop_censored) ,
                      summary(sim_nop_censored)$theta)
  FEM2G_table[,s] <- c(coef(sim_nop_truncatedG),
                       summary(sim_nop_truncatedG)$theta)
  FEM2T_table[,s] <- c(coef(sim_nop_truncatedT),
                       summary(sim_nop_truncatedT)$theta)

  ###########################################################################
  # [ADDED 5] Ordinal deductible -> SDM1/2/3 -> FDM1/2 for replication s
  ###########################################################################
  dat <- as_model_data(sim_dat.pol, sim_dat.clm, sim_dat.clmpos)
  proportions <- cumsum(tabulate(dat$pol$DeductLevel, nbins=5L))[1:4] / n
  D0 <- d_pack(qlogis(pclip(proportions)), c(0,0))
  D_fit <- fit_copula_stage("D", D0, list(), dat)
  D_table[,s] <- unname(d_unpack(D_fit$par))

  # Copula stages start from your fitted SEM and FEM parameters.
  # log(alpha-1), log(r), atanh(rho) are unconstrained optimizer coordinates.
  fitted_z <- list(D=D_fit$par)
  stage_fits <- list()
  for (model in c("SDM1","SDM2","SDM3")) {
    base <- switch(model,
                   SDM1=SEM1_table[,s],
                   SDM2=SEM2_table[,s],
                   SDM3=SEM3_table[,s])
    z0 <- c(base[1:3], log(base[4]-1), 0)
    f <- fit_copula_stage(model, z0, fitted_z, dat)
    fitted_z[[model]] <- f$par
    stage_fits[[model]] <- f
    theta <- s_unpack(f$par, TRUE)
    if (model=="SDM1") SDM1_table[,s] <- unname(theta)
    if (model=="SDM2") SDM2_table[,s] <- unname(theta)
    if (model=="SDM3") SDM3_table[,s] <- unname(theta)
  }

  for (model in c("FDM1", "FDM2")) {
    base <- if (model == "FDM1") FEM1_table[,s] else FEM2T_table[,s]
    z0 <- c(base[1:3], log(base[4]), 0)
    f <- fit_copula_stage(model, z0, fitted_z, dat)
    fitted_z[[model]] <- f$par
    stage_fits[[model]] <- f
    theta <- f_unpack(f$par, TRUE)
    if (model == "FDM1") FDM1_table[,s] <- unname(theta)
    if (model == "FDM2") FDM2_table[,s] <- unname(theta)
  }
  copula_convergence[s,] <- c(D_fit$convergence,
       vapply(c("SDM1","SDM2","SDM3","FDM1","FDM2"),
              function(nm) stage_fits[[nm]]$convergence, numeric(1)))

  # Optional selectable Godambe computations: any estimation stage.
  # The score and Hessian cache is shared across targets within a replication.
  if (RUN_SANDWICH && s %in% SANDWICH_REPS && length(SANDWICH_TARGETS)) {
    fitz <- fitted_z
    for (nm in c("SEM1", "SEM2", "SEM3")) {
      vv <- switch(nm, SEM1 = SEM1_table[, s], SEM2 = SEM2_table[, s],
                   SEM3 = SEM3_table[, s])
      fitz[[nm]] <- c(vv[1:3], log(vv[4] - 1))
    }
    for (nm in c("FEM1", "FEM2G", "FEM2")) {
      vv <- switch(nm, FEM1 = FEM1_table[, s], FEM2G = FEM2G_table[, s],
                   FEM2 = FEM2T_table[, s])
      fitz[[nm]] <- c(vv[1:3], log(vv[4]))
    }
    estimate_obj <- lapply(fitz, function(z) list(z = z))
    cache <- new.env(parent = emptyenv())
    for (target in SANDWICH_TARGETS) {
      sandwich_start <- proc.time()[["elapsed"]]
      rr <- tryCatch(godambe(target, estimate_obj, dat, cache = cache),
                     error = function(e) {
                       warning(sprintf("Godambe %s rep %d: %s", target, s,
                                       conditionMessage(e)))
                       NULL
                     })
      SANDWICH_TIMING <- rbind(SANDWICH_TIMING,
        data.frame(rep = s, model = target,
                   elapsed_seconds = proc.time()[["elapsed"]] - sandwich_start,
                   succeeded = !is.null(rr)))
      cat(sprintf("Godambe rep %d %s: %.1f sec, success=%s\n",
                  s, target, tail(SANDWICH_TIMING$elapsed_seconds, 1L),
                  !is.null(rr)))
      if (!is.null(rr)) {
        rr$rep <- s
        sandwich_rows[[length(sandwich_rows) + 1L]] <- rr
      }
    }
  }

  cat(sprintf("Case 4 rep %d: D/SDM/FDM convergence = %s\n", s,
              paste(copula_convergence[s,], collapse=",")))
}
# Monte Carlo simulation estimation loop finished; proceed to premium assessment.

#######################################################################################################################################################################


###############################################################################
# [ADDED 6] Six aggregate insured-loss premium comparisons
# Paired models by reporting scheme; all estimated in your original loops.
###############################################################################

# The prior helper freq_mean() uses a tail-sum for conditional NB frequencies;
# sev_paid_mean() uses the conditional Lomax mean with Gaussian atom weighting.
# Under rho=0, these reduce to E[N]=exp(X beta_N) and the Lomax formula.

PREMIUM_PAIRS <- data.frame(
  scenario=c("Full","Full","Censored","Censored","Truncated","Truncated"),
  freq_model=c("FEM1","FDM1","FEM1","FDM1","FEM2T","FDM2"),
  sev_model =c("SEM1","SDM1","SEM2","SDM2","SEM3","SDM3"),
  stringsAsFactors=FALSE
)
PREMIUM_PAIRS$model <- paste(PREMIUM_PAIRS$freq_model,
                             PREMIUM_PAIRS$sev_model, sep="+")

get_par <- function(name, s) {
  switch(name,
    FEM1=FEM1_table[,s], FEM2T=FEM2T_table[,s],
    SEM1=SEM1_table[,s], SEM2=SEM2_table[,s], SEM3=SEM3_table[,s],
    FDM1=FDM1_table[,s], FDM2=FDM2_table[,s],
    SDM1=SDM1_table[,s], SDM2=SDM2_table[,s], SDM3=SDM3_table[,s],
    stop("Unknown model: ",name))
}

# UPDATED (METRICS): four criteria only.
# Observed realized aggregate loss: MAE_obs.
# True conditional premium: MAE_premium, RMSRE_premium, Relative_bias_premium.
# Relative bias is the total predicted premium minus total true premium,
# divided by total true premium (not mean of policy-level relative errors).
premium_metrics <- function(observed, truth, predicted) {
  keys <- c("MAE_obs", "MAE_premium", "RMSRE_premium",
            "Relative_bias_premium")
  fail <- setNames(rep(NA_real_, length(keys)), keys)
  if (!length(observed) || length(truth) != length(observed) ||
      length(predicted) != length(observed) ||
      any(!is.finite(c(observed, truth, predicted))) || any(truth <= 0))
    return(fail)

  rel_true <- (predicted - truth) / truth
  c(MAE_obs = mean(abs(predicted - observed)),
    MAE_premium = mean(abs(predicted - truth)),
    RMSRE_premium = sqrt(mean(rel_true^2)),
    Relative_bias_premium = (sum(predicted) - sum(truth)) / sum(truth))
}

model_names <- PREMIUM_PAIRS$model
metric_names <- c("MAE_obs", "MAE_premium", "RMSRE_premium",
                  "Relative_bias_premium")
S_metrics <- array(NA_real_, c(t, length(model_names),length(metric_names)),
                   dimnames=list(NULL,model_names,metric_names))
S_metrics_by_D <- array(NA_real_,
  c(t,length(ded_levels),length(model_names),length(metric_names)),
  dimnames=list(NULL,as.character(ded_levels),model_names,metric_names))
example_predictions <- NULL

for (s in seq_len(t)) {
  pol <- sim_dat.pol.list[[s]]
  dat <- as_model_data(pol, sim_dat.clm.list[[s]],sim_dat.clmpos.list[[s]])
  pol <- dat$pol
  X <- as.matrix(pol[,c("int","X1","X2")])
  atomD <- atom_bounds(pol, d_pack(D_table[1:4,s],D_table[5:6,s]))
  atom0 <- atom_bounds(pol, NULL)

  # DGP-matched oracle E[S|D,X] using TRUE parameters and the
  # TRUE discrete deductible atom, not the estimated D parameters.
  # Conditional independence N _|_ (Y_j) | (D,X) is built into the DGP.
  atom_true <- atom_bounds(pol, d_pack(kappa, beta_D))
  true_freq <- c(beta_N, r, rho_N)
  true_sev <- c(beta_Y, alpha2, rho_Y)
  oracle <- freq_mean(true_freq, pol, atom_true) *
    sev_paid_mean(true_sev, pol, atom_true)
  if (any(!is.finite(oracle)) || any(oracle <= 0))
    stop("Invalid oracle E[S|D,X] in replication ", s)

  # Recompute frequencies and severities from FITTED values, never from the
  # mu_N / lambda_Y left over from the most recent simulation iteration.
  f_models <- unique(PREMIUM_PAIRS$freq_model)
  s_models <- unique(PREMIUM_PAIRS$sev_model)
  freq_pred <- setNames(lapply(f_models,function(name) {
    at <- if (startsWith(name,"FDM")) atomD else atom0
    freq_mean(get_par(name,s), pol, at)
  }),f_models)
  sev_pred <- setNames(lapply(s_models,function(name) {
    at <- if (startsWith(name,"SDM")) atomD else atom0
    sev_paid_mean(get_par(name,s),pol,at)
  }),s_models)

  if (s==1L) example_predictions <- data.frame(
    policy_id=pol$policy_id,X1=pol$X1,X2=pol$X2,
    Deduct=pol$Deduct,S_observed=pol$S,S_true=oracle)

  for (j in seq_len(nrow(PREMIUM_PAIRS))) {
    fm <- PREMIUM_PAIRS$freq_model[j]; sm <- PREMIUM_PAIRS$sev_model[j]
    model <- PREMIUM_PAIRS$model[j]
    premium <- as.numeric(freq_pred[[fm]] * sev_pred[[sm]])
    S_metrics[s,model,] <- premium_metrics(pol$S,oracle,premium)
    if (s==1L) example_predictions[[model]] <- premium
    for (k in seq_along(ded_levels)) {
      idx <- pol$Deduct==ded_levels[k]
      if (any(idx))
        S_metrics_by_D[s,k,model,] <- premium_metrics(
          pol$S[idx],oracle[idx],premium[idx])
    }
  }
  if (s %% max(1L,floor(t/10))==0L || s==1L)
    cat(sprintf("Premium assessment: %d/%d replications\n",s,t))
}

# Summary: equal-weight replication means; MC_SD is SD of metric across reps.
overall_results <- expand.grid(Model=model_names,Metric=metric_names,
                                stringsAsFactors=FALSE)
overall_results$Mean <- NA_real_
overall_results$MC_SD <- NA_real_
overall_results$Valid_reps <- 0L
for (j in seq_len(nrow(overall_results))) {
  x <- S_metrics[,overall_results$Model[j],overall_results$Metric[j]]
  ok <- is.finite(x)
  overall_results$Valid_reps[j] <- sum(ok)
  if (any(ok)) overall_results$Mean[j] <- mean(x[ok])
  if (sum(ok)>1L) overall_results$MC_SD[j] <- sd(x[ok])
}

by_deductible_results <- expand.grid(
  Deductible=ded_levels,Model=model_names,Metric=metric_names,
  stringsAsFactors=FALSE)
by_deductible_results$Mean <- NA_real_
by_deductible_results$MC_SD <- NA_real_
by_deductible_results$Valid_reps <- 0L
for (j in seq_len(nrow(by_deductible_results))) {
  x <- S_metrics_by_D[,
    as.character(by_deductible_results$Deductible[j]),
    by_deductible_results$Model[j],by_deductible_results$Metric[j]]
  ok <- is.finite(x)
  by_deductible_results$Valid_reps[j] <- sum(ok)
  if (any(ok)) by_deductible_results$Mean[j] <- mean(x[ok])
  if (sum(ok)>1L) by_deductible_results$MC_SD[j] <- sd(x[ok])
}
print(overall_results)
write.csv(overall_results,"Case4_Aggregate_Overall.csv",row.names=FALSE)
write.csv(by_deductible_results,"Case4_Aggregate_ByDeductible.csv",row.names=FALSE)
write.csv(example_predictions,"Case4_Example_Predictions.csv",row.names=FALSE)
write.csv(copula_convergence,"Case4_Copula_Convergence.csv",row.names=FALSE)

# Graphs: same overview as your original, now all six lines.
plot_aggregate_by_deductible <- function(metric) {
  mat <- sapply(model_names, function(m) {
    vapply(ded_levels,function(d) {
      by_deductible_results$Mean[by_deductible_results$Deductible==d &
        by_deductible_results$Model==m & by_deductible_results$Metric==metric]
    },numeric(1))
  })
  yrange <- range(c(0,mat),finite=TRUE)
  if (!all(is.finite(yrange)) || diff(yrange)<1e-12) yrange <- c(0,1)
  matplot(ded_levels,mat,type="b",pch=rep(c(16,1),3),lty=rep(c(1,2),3),
          col=rep(1:3,each=2),log="x",xaxt="n",ylim=yrange,
          xlab="Per-occurrence deductible",ylab=metric,main=paste("Case 4",metric))
  axis(1,at=ded_levels,labels=ded_levels)
  legend("topright",legend=model_names,pch=rep(c(16,1),3),
         lty=rep(c(1,2),3),col=rep(1:3,each=2),cex=.68,bty="n")
}
pdf("Case4_Aggregate_ByDeductible.pdf",width=10,height=6)
for (metric in metric_names) plot_aggregate_by_deductible(metric)
dev.off()

# Parameter boxplots compare marginal vs copula estimates.
plot_parameter_box <- function(models, indices, pname, file) {
  pdf(file,width=11,height=5.8,onefile=TRUE)
  for (i in seq_along(pname)) {
    z <- lapply(seq_along(models),function(j) {
      M <- models[[j]]
      if (indices[j,i]<1L) return(numeric())
      v <- M[indices[j,i],]
      v[is.finite(v)]
    })
    names(z) <- names(models)
    z <- z[lengths(z)>0L]
    if (!length(z)) next
    boxplot(z,main=paste("Case 4:",pname[i]),ylab="Estimated parameter",
            las=2,col="grey92")
    truth <- switch(pname[i],beta0=if (grepl("severity",file)) beta_Y[1] else beta_N[1],
      beta1=if (grepl("severity",file)) beta_Y[2] else beta_N[2],
      beta2=if (grepl("severity",file)) beta_Y[3] else beta_N[3],
      alpha=alpha2,r=r,
      rho=if (grepl("severity",file)) rho_Y else rho_N)
    abline(h=truth,col="red",lty=2,lwd=2)
  }
  dev.off()
}
sev_models <- list(SEM1=SEM1_table,SEM2=SEM2_table,SEM3=SEM3_table,
                   SDM1=SDM1_table,SDM2=SDM2_table,SDM3=SDM3_table)
freq_models <- list(FEM1=FEM1_table,FEM2=FEM2T_table,
                    FDM1=FDM1_table,FDM2=FDM2_table)
si <- matrix(rep(1:5,each=6),nrow=6)
si[1:3,5] <- 0
fi <- matrix(rep(1:5,each=4),nrow=4)
fi[1:2,5] <- 0
plot_parameter_box(sev_models,si,c("beta0","beta1","beta2","alpha","rho"),
                   "Case4_Boxplots_Severity.pdf")
plot_parameter_box(freq_models,fi,c("beta0","beta1","beta2","r","rho"),
                   "Case4_Boxplots_Frequency.pdf")

# Directly inspect the Monte Carlo parameter estimates in the original tables.
saveRDS(list(D_table=D_table, FEM1_table=FEM1_table,FEM2G_table=FEM2G_table,
             FEM2T_table=FEM2T_table,SEM1_table=SEM1_table,SEM2_table=SEM2_table,
             SEM3_table=SEM3_table,FDM1_table=FDM1_table,FDM2_table=FDM2_table,
             SDM1_table=SDM1_table,SDM2_table=SDM2_table,SDM3_table=SDM3_table),
        "Case4_Parameter_Tables.rds")

# MC parameter summary: mean, empirical SD, bias, RMSE by model and parameter.
all_tables <- list(D=D_table,FEM1=FEM1_table,FEM2G=FEM2G_table,
  FEM2T=FEM2T_table,SEM1=SEM1_table,SEM2=SEM2_table,SEM3=SEM3_table,
  SDM1=SDM1_table,SDM2=SDM2_table,SDM3=SDM3_table,
  FDM1=FDM1_table,FDM2=FDM2_table)
par_summary_rows <- list()
for (name in names(all_tables)) {
  M <- all_tables[[name]]
  tr <- if (name=="D") c(kappa,beta_D) else if (startsWith(name,"S"))
    c(beta_Y,alpha2,if (startsWith(name,"SDM")) rho_Y) else
    c(beta_N,r,if (startsWith(name,"FDM")) rho_N)
  labels <- if (name=="D") c(paste0("kappa",1:4),"betaD1","betaD2") else
    c("beta0","beta1","beta2",if (startsWith(name,"S")) "alpha" else "r",
      if (name %in% c("FDM1","FDM2","SDM1","SDM2","SDM3")) "rho")
  for (j in seq_len(nrow(M))) {
    estimates <- M[j,]
    estimates <- estimates[is.finite(estimates)]
    par_summary_rows[[length(par_summary_rows)+1L]] <- data.frame(
      Model=name,Parameter=labels[j],True=tr[j],
      Mean=if (length(estimates)) mean(estimates) else NA_real_,
      Empirical_SD=if (length(estimates)>1) sd(estimates) else NA_real_,
      Bias=if (length(estimates)) mean(estimates)-tr[j] else NA_real_,
      RMSE=if (length(estimates)) sqrt(mean((estimates-tr[j])^2)) else NA_real_,
      Valid_reps=length(estimates))
  }
}
parameter_summary <- do.call(rbind,par_summary_rows)
write.csv(parameter_summary,"Case4_Parameter_Summary.csv",row.names=FALSE)

if (nrow(SANDWICH_TIMING)) {
  write.csv(SANDWICH_TIMING, "Case4_Godambe_Runtime.csv", row.names=FALSE)
}
if (length(sandwich_rows)) {
  se_results <- do.call(rbind,sandwich_rows)
  write.csv(se_results,"Case4_Godambe_SE.csv",row.names=FALSE)
}


###############################################################################
# Monte Carlo assessment of sequential standard errors
# Naive vs stage-wise robust vs full Godambe
###############################################################################

if (length(sandwich_rows) > 0L) {
  
  se_results <- do.call(rbind, sandwich_rows)
  
  models <- unique(SANDWICH_TARGETS)
  
  results <- list()
  row_id <- 0L
  
  for (mod in models) {
    
    dat_mod <- se_results[
      se_results$model == mod, ,
      drop = FALSE
    ]
    
    if (nrow(dat_mod) == 0L) next
    
    for (par in unique(dat_mod$parameter)) {
      
      z <- dat_mod[
        dat_mod$parameter == par, ,
        drop = FALSE
      ]
      
      # Only complete and valid replications
      valid <- is.finite(z$estimate) &
        is.finite(z$true) &
        is.finite(z$naive_se) &
        is.finite(z$stage_robust_se) &
        is.finite(z$godambe_se) &
        z$naive_se > 0 &
        z$stage_robust_se > 0 &
        z$godambe_se > 0
      
      z <- z[valid, , drop = FALSE]
      
      R_valid <- nrow(z)
      
      if (R_valid == 0L) next
      
      # Empirical Monte Carlo SD
      empirical_sd <- if (R_valid >= 2L) {
        sd(z$estimate)
      } else {
        NA_real_
      }
      
      # Mean estimated standard errors
      mean_naive <- mean(z$naive_se)
      mean_robust <- mean(z$stage_robust_se)
      mean_godambe <- mean(z$godambe_se)
      
      # Empirical coverage of nominal 95% intervals
      q <- qnorm(0.975)
      
      coverage_naive <- mean(
        abs(z$estimate - z$true) <=
          q * z$naive_se
      )
      
      coverage_robust <- mean(
        abs(z$estimate - z$true) <=
          q * z$stage_robust_se
      )
      
      coverage_godambe <- mean(
        abs(z$estimate - z$true) <=
          q * z$godambe_se
      )
      
      row_id <- row_id + 1L
      
      results[[row_id]] <- data.frame(
        Model = mod,
        Parameter = par,
        True = z$true[1],
        Empirical_SD = empirical_sd,
        Mean_SE_Naive = mean_naive,
        Mean_SE_Robust = mean_robust,
        Mean_SE_Godambe = mean_godambe,
        Coverage_Naive = coverage_naive,
        Coverage_Robust = coverage_robust,
        Coverage_Godambe = coverage_godambe,
        Valid_reps = R_valid
      )
    }
  }
  
  if (length(results) > 0L) {
    
    coverage_summary <- do.call(rbind, results)
    
    print(coverage_summary, row.names = FALSE)
    
    write.csv(
      coverage_summary,
      "Case4_Godambe_Coverage_Summary.csv",
      row.names = FALSE
    )
  }
}


# [ADDED 7] Algebraic identities evaluated at the INDEPENDENCE SUBMODEL
# (rho=0). These checks are valid even when the simulated DGP has nonzero
# rho_N or rho_Y; they are not tests that the actual data are independent.
pol <- as_model_data(sim_dat.pol.list[[1]],sim_dat.clm.list[[1]],
                     sim_dat.clmpos.list[[1]])$pol
xx <- as.matrix(pol[,c("int","X1","X2")])
mu <- drop(exp(xx %*% beta_N))
lam <- drop(exp(xx %*% beta_Y))
p <- (1+pol$Deduct/lam)^(-alpha2)
from_censored <- mu*lam/(alpha2-1)*(1+pol$Deduct/lam)^(1-alpha2)
from_truncated <- (mu*p)*(lam+pol$Deduct)/(alpha2-1)
cat("Maximum premium formula identity error:",
    max(abs(from_censored-from_truncated)),"\n")

# Null-copula likelihood reduction check on replication 1:
# at rho=0 SDM3 = SEM3 and FDM2 = FEM2 (up to numerical tolerance).
test_dat <- as_model_data(sim_dat.pol.list[[1]], sim_dat.clm.list[[1]],
                          sim_dat.clmpos.list[[1]])
dz <- d_pack(D_table[1:4,1],D_table[5:6,1])
sz <- c(SEM3_table[1:3,1], log(SEM3_table[4,1]-1))
fz <- c(FEM2T_table[1:3,1], log(FEM2T_table[4,1]))
null_check <- data.frame(
  Identity=c("SEM3 = SDM3(rhoY=0)","FEM2 = FDM2(rhoN=rhoY=0)"),
  Difference=c(
    stage_nll("SDM3",c(sz,0),list(D=dz),test_dat) -
      stage_nll("SEM3",sz,list(),test_dat),
    stage_nll("FDM2",c(fz,0),list(D=dz,SDM3=c(sz,0)),test_dat) -
      stage_nll("FEM2",fz,list(SEM3=sz),test_dat)
  ))
print(null_check)
write.csv(null_check,"Case4_Null_Copula_Check.csv",row.names=FALSE)

#save.image(file="Case4_2nd_round.Rdata")

#### Deductible tex ##### 

# ============================================================
# Export MAE and RMSRE directly to Overleaf
# ============================================================

case_id <- 4L  # Case 1, 2, 3, 4에 맞게 변경

# Overleaf 표에 넣을 모형 순서 및 표시 이름
model_order <- PREMIUM_PAIRS$model

# 예: 6개 모형이라면 XR1~ER3
# Case 3처럼 8개 모형이면 XR1~XR4, ER1~ER4
#n_models <- length(model_order)
#n_xr <- ceiling(n_models / 2)

#row_labels <- c(
# paste0("XR", seq_len(n_xr)),
#  paste0("ER", seq_len(n_models - n_xr))
#)

deductibles <- c(500, 1000, 2500, 5000, 10000)

# 출력 자릿수
digits_MAE <- 2
digits_RMSRE <- 2

get_metric <- function(model, metric, d) {
  idx <- with(
    by_deductible_results,
    Model == model &
      Metric == metric &
      Deductible == d
  )
  
  values <- by_deductible_results$Mean[idx]
  
  if (length(values) != 1L) {
    stop("Missing or duplicated result: ", model, ", ", metric, ", ", d)
  }
  
  values
}

fmt <- function(x, digits) {
  ifelse(
    is.finite(x),
    sprintf(paste0("%.", digits, "f"), x),
    "--"
  )
}

latex_rows <- length(model_order)

for (j in seq_along(model_order)) {
  
  model <- model_order[j]
  
  mae <- vapply(
    deductibles,
    function(d) get_metric(model, "MAE_premium", d),
    numeric(1)
  )
  
  rmsre <- vapply(
    deductibles,
    function(d) get_metric(model, "RMSRE_premium", d),
    numeric(1)
  )
  
  vals <- c(
    fmt(mae, digits_MAE),
    fmt(rmsre, digits_RMSRE)
  )
  
  prefix <- if (j == 1L) {
    paste0("\\multirow{", model, "}{*}{Case ", case_id, "}")
  } else {
    ""
  }
  
  latex_rows[j] <- paste0(
    prefix,
    " & ",
    model_order[j],
    " & ",
    paste(vals, collapse = " & "),
    " \\\\"
  )
}

cat(paste(latex_rows, collapse = "\n"), "\n")

### Plotting new paremter

# Case 1: 2 x 5 parameter-estimate boxplots (base R)
# Paste after all Monte Carlo estimation loops, replacing the original
# `plot_parameter_box` through its two calls. No CSV files are needed.
# Expected objects: FEM1_table, FEM2G_table, FEM2T_table, FDM1_table,
# FDM2_table, SEM1_table, SEM2_table, SEM3_table, SDM1_table,
# SDM2_table, SDM3_table, beta_N, beta_Y, alpha2, r, rho_N, rho_Y.

pdf_boxes_case4 <- function(file = "Case4_Styled_Parameter_Boxplots.pdf") {
  severity <- list(SEM1=SEM1_table, SEM2=SEM2_table, SEM3=SEM3_table,
                   SDM1=SDM1_table, SDM2=SDM2_table, SDM3=SDM3_table)
  frequency <- list(FEM1=FEM1_table, FEM2=FEM2T_table, 
                    FDM1=FDM1_table, FDM2=FDM2_table)
  
  draw_box <- function(tables, parameter_index, title, true_value) {
    vals <- lapply(tables, function(M) {
      if (nrow(M) < parameter_index) return(numeric(0))
      x <- as.numeric(M[parameter_index, ])
      x[is.finite(x)]
    })
    present <- lengths(vals) > 0L
    if (!any(present)) {
      plot.new(); title(main=title)
      text(0.5, 0.5, "No finite estimates")
      return(invisible(NULL))
    }
    vals <- vals[present]
    model_names <- names(vals)
    cols <- ifelse(grepl("DM", model_names), "lightsteelblue", "grey90")
    # Include the true value in the y-axis limits even when estimates are far away.
    ylim <- range(c(unlist(vals, use.names=FALSE), true_value), finite=TRUE)
    if (diff(ylim) < 1e-10) ylim <- ylim + c(-0.5, 0.5)
    boxplot(vals, names=model_names, las=2, main=title,
            outline=TRUE, pch=20, cex=0.5, col=cols,
            ylim=ylim, ylab="")
    abline(h=true_value, col="red", lty=2, lwd=1.5)
  }
  
  pdf(file, width=11, height=6.5)
  oldpar <- par(no.readonly=TRUE)
  on.exit({par(oldpar); dev.off()}, add=TRUE)
  par(mfrow=c(2,5), mar=c(4.2,3.2,2.2,0.6),
      mgp=c(2,0.7,0), cex.axis=0.8)
  
  # Top row: severity parameters
  draw_box(severity, 1, expression(beta[Y*","*0]), beta_Y[1])
  draw_box(severity, 2, expression(beta[Y*","*1]), beta_Y[2])
  draw_box(severity, 3, expression(beta[Y*","*2]), beta_Y[3])
  draw_box(severity, 4, expression(alpha), alpha2)
  draw_box(severity[c("SDM1","SDM2","SDM3")], 5,
           expression(rho[Y]), rho_Y)
  
  # Bottom row: frequency parameters
  draw_box(frequency, 1, expression(beta[N*","*0]), beta_N[1])
  draw_box(frequency, 2, expression(beta[N*","*1]), beta_N[2])
  draw_box(frequency, 3, expression(beta[N*","*2]), beta_N[3])
  draw_box(frequency, 4, expression(r), r)
  draw_box(frequency[c("FDM1","FDM2")], 5,
           expression(rho[N]), rho_N)
  invisible(file)
}

pdf_boxes_case4()

######### Estimates of the deductible parameters #################

# ============================================================
# Case 1: Ordinal deductible parameter summary
# ============================================================

true_D <- c(
  kappa1 = 0.1,
  kappa2 = 1.0,
  kappa3 = 1.8,
  kappa4 = 2.6,
  beta_D1 = -0.5,
  beta_D2 = -0.5
)

# D_table: parameters x Monte Carlo replications
D_est <- as.matrix(D_table)

# Ensure correct orientation
if (nrow(D_est) != length(true_D)) {
  stop("D_table must have 6 rows (one per parameter).")
}

rownames(D_est) <- names(true_D)

# Remove replications with missing parameter estimates
valid <- colSums(is.finite(D_est)) == length(true_D)
D_est <- D_est[, valid, drop = FALSE]

R_valid <- ncol(D_est)

if (R_valid < 2L) {
  stop("At least two valid replications are required.")
}

# Summary statistics
mean_est <- rowMeans(D_est)

bias_est <- mean_est - true_D

rmse_est <- sqrt(
  rowMeans(
    sweep(D_est, 1, true_D, "-")^2
  )
)

mc_se <- apply(D_est, 1, sd) / sqrt(R_valid)

# Final summary
D_summary <- data.frame(
  Parameter = names(true_D),
  True = unname(true_D),
  Mean = unname(mean_est),
  RMSE = unname(rmse_est),
  SE = unname(mc_se),
  Bias = unname(bias_est),
  row.names = NULL
)

cat("Valid Monte Carlo replications:", R_valid, "\n")

print(
  D_summary,
  digits = 4,
  row.names = FALSE
)


# ============================================================
# CASE 1: Export parameter estimation tables to LaTeX
# Run AFTER all Monte Carlo replications have finished.
# ============================================================

# 1. Model estimates stored as parameters x replications
estimates <- list(
  SEM1 = SEM1_table,
  SEM2 = SEM2_table,
  SEM3 = SEM3_table,
  FEM1 = FEM1_table,
  FEM2T = FEM2T_table,
  SDM1 = SDM1_table,
  SDM2 = SDM2_table,
  SDM3 = SDM3_table,
  FDM1 = FDM1_table,
  FDM2 = FDM2_table
)

# 2. True parameters
true_values <- list(
  severity = c(beta_Y, alpha2, rho_Y),
  frequency = c(beta_N, r, rho_N)
)

# 3. Summarize estimates for each model
summarize_model <- function(model) {
  
  M <- as.matrix(estimates[[model]])
  
  is_severity <- startsWith(model, "S")
  
  truth <- if (is_severity) {
    true_values$severity
  } else {
    true_values$frequency
  }
  
  # Marginal models have 4 parameters;
  # deductible-dependent models have 5.
  p <- if (grepl("DM", model)) 5L else 4L
  
  M <- M[seq_len(p), , drop = FALSE]
  truth <- truth[seq_len(p)]
  
  # Retain only replications with finite estimates
  valid <- colSums(is.finite(M)) == p
  M <- M[, valid, drop = FALSE]
  
  R_valid <- ncol(M)
  
  if (R_valid < 2L) {
    stop("Insufficient valid replications for ", model)
  }
  
  mean_est <- rowMeans(M)
  bias_est <- mean_est - truth
  
  rmse_est <- sqrt(
    rowMeans(sweep(M, 1, truth, "-")^2)
  )
  
  mcse_est <- apply(M, 1, sd) / sqrt(R_valid)
  
  data.frame(
    Parameter = seq_len(p),
    Mean = mean_est,
    RMSE = rmse_est,
    MCSE = mcse_est,
    Bias = bias_est,
    Valid_reps = R_valid
  )
}

summaries <- lapply(names(estimates), summarize_model)
names(summaries) <- names(estimates)

# 4. Generate a LaTeX tabular block
make_latex_tabular <- function(models, parameter_labels) {
  
  n <- length(models)
  p <- length(parameter_labels)
  
  fmt <- function(x) sprintf("%.4f", x)
  
  out <- character()
  
  add <- function(...) {
    out <<- c(out, paste0(...))
  }
  
  add("\\begin{tabular}{l", paste(rep("cc", n), collapse = ""), "}")
  add("\\toprule")
  
  header <- paste0(
    " & ",
    paste(
      sprintf("\\multicolumn{2}{c}{%s}", models),
      collapse = " & "
    ),
    " \\\\"
  )
  add(header)
  
  cmid <- vapply(seq_len(n), function(j) {
    a <- 2 * j
    b <- a + 1
    sprintf("\\cmidrule(lr){%d-%d}", a, b)
  }, character(1))
  add(paste(cmid, collapse = " "))
  
  add(
    " & ",
    paste(rep("Mean & RMSE", n), collapse = " & "),
    " \\\\"
  )
  add(
    " & ",
    paste(rep("(MCSE) & (Bias)", n), collapse = " & "),
    " \\\\"
  )
  add("\\midrule")
  
  for (i in seq_len(p)) {
    
    main_values <- character()
    sub_values <- character()
    
    for (model in models) {
      
      S <- summaries[[model]]
      
      if (i <= nrow(S)) {
        main_values <- c(
          main_values,
          fmt(S$Mean[i]),
          fmt(S$RMSE[i])
        )
        
        sub_values <- c(
          sub_values,
          paste0("(", fmt(S$MCSE[i]), ")"),
          paste0("(", fmt(S$Bias[i]), ")")
        )
      } else {
        main_values <- c(main_values, "", "")
        sub_values <- c(sub_values, "", "")
      }
    }
    
    add(
      parameter_labels[i], " & ",
      paste(main_values, collapse = " & "),
      " \\\\"
    )
    
    add(
      " & ",
      paste(sub_values, collapse = " & "),
      " \\\\"
    )
    
    if (i < p) add("\\midrule")
  }
  
  add("\\bottomrule")
  add("\\end{tabular}")
  
  paste(out, collapse = "\n")
}

# 5. Parameter labels
sev_labels <- c(
  "$\\beta_{Y,0}$",
  "$\\beta_{Y,1}$",
  "$\\beta_{Y,2}$",
  "$\\alpha$",
  "$\\rho_Y$"
)

freq_labels <- c(
  "$\\beta_{N,0}$",
  "$\\beta_{N,1}$",
  "$\\beta_{N,2}$",
  "$r$",
  "$\\rho_N$"
)

# 6. Generate the four tabular blocks
tex_SEM <- make_latex_tabular(
  c("SEM1", "SEM2", "SEM3"),
  sev_labels[1:4]
)

tex_FEM <- make_latex_tabular(
  c("FEM1", "FEM2T"),
  freq_labels[1:4]
)

tex_SDM <- make_latex_tabular(
  c("SDM1", "SDM2", "SDM3"),
  sev_labels
)

tex_FDM <- make_latex_tabular(
  c("FDM1", "FDM2"),
  freq_labels
)

# 7. Print to console
cat("\n========== SEM ==========\n", tex_SEM)
cat("\n\n========== FEM ==========\n", tex_FEM)
cat("\n\n========== SDM ==========\n", tex_SDM)
cat("\n\n========== FDM ==========\n", tex_FDM, "\n")

# 8. Save LaTeX blocks
writeLines(tex_SEM, "Case1_SEM.tex")
writeLines(tex_FEM, "Case1_FEM.tex")
writeLines(tex_SDM, "Case1_SDM.tex")
writeLines(tex_FDM, "Case1_FDM.tex")

# 9. Check valid replications
valid_rep_summary <- data.frame(
  Model = names(summaries),
  Valid_reps = vapply(
    summaries,
    function(x) x$Valid_reps[1],
    numeric(1)
  )
)

print(valid_rep_summary, row.names = FALSE)

load(file="Case4_2nd_metric.Rdata")
