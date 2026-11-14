setwd('/Users/idongha/Library/CloudStorage/OneDrive-SimonFraserUniversity(1sfu)/SFU/Research/Deductible ratemaking/Submitted/2nd/New Empirical/261113/0609')
################################################################################
# INTEGRATED EMPIRICAL ANALYSIS | Case 1 and 4 + 2010 OOS sensitivity
# Based on endogenous_deductible_pipeline.R supplied by the user.
# Run: source("Empirical_Severity_First_Premium.R")
# Dependencies: MASS, VineCopula, numDeriv; optional knitr/kableExtra.
# Data: poldat.RData, claimdat.RData. No data files bundled with this code.
#
# Workflow
#  1. D + SEM1/SEM2/SEM3; common-axis mu-versus-log(Coverage) figure.
#  2. AIC-select SCM1/2/3 and SDM1/2/3 before fitting any frequency model.
#  3. Fit FEM1, FEM2_SEM3, FEM2G_SEM1, FEM2G_cont, FEM2G_disc.
#  4. AIC-select FCM1/2 and FDM1/2 (without assuming Case 4 a priori).
#  5. Compare only XR1-3 (Case 1), ER1-3 (Cases 4c and 4d).
#  6. Optional Q-Q, frequency GOF, and full-chain Godambe inference.
#
# CAUTION: negative slope of SEM3 mu against log(Coverage) is suggestive,
#          NOT proof of endogeneity. Severity GB2 first moment is infinite
#          if alpha2 <= sigma. Unlimited premium is then unavailable.
#          PRIMARY observed-payment target: per-occurrence ground-up coverage
#          C_i, converted from the policy Coverage field into dollars. The
#          contractual interpretation is an explicit assumption, not verified.
################################################################################

################################################################################
##  Insurance Ratemaking with Endogenous Deductibles -- revised empirical pipeline
##
##  Merged / streamlined replacement for the five original scripts
##     Marginal_models__exogenous_deductible__T1_2_3_4_5_F1_2.R
##     copula_selection_contDedSev.R      copula_selection_discDedSev.R
##     copula_selection_contDedFreq.R     copula_selection_discDedFreq.R
##
##  What is new
##   * Every likelihood is written ONCE as a vector of per-unit log-likelihood
##     contributions, with the earlier-stage parameters (deductible model,
##     severity model) passed as ARGUMENTS instead of being read from globals.
##     This is what allows the Godambe (sandwich) covariance of the sequential
##     estimator to be computed (Section 3 below), and it removes the
##     9 copulas x 3 scenarios x 2 deductible types code duplication.
##   * Standard errors for the final models come in three flavours:
##       naive     : stage-wise Hessian, earlier stages treated as known
##       robust    : stage-wise sandwich  H_kk^{-1} J_kk H_kk^{-1}
##       sandwich  : full Godambe sandwich  H^{-1} J H^{-T}  (accounts for the
##                   uncertainty in the parameters estimated at earlier stages)
##
##  Layout
##   0. settings, data, design matrices
##   1. distributions and copula helpers
##   2. per-unit log-likelihood contributions (deductible / severity / frequency)
##   3. estimation of one stage + Godambe sandwich for a chain of stages
##   4. marginal models (D, SEM1-3, FEM1, FEM2, FEM2G) and severity family table
##   5. copula selection (9 families x {cont, disc} x {complete, censored, truncated})
##   6. final joint models (SCM/SDM 1-3, FCM/FDM 1-2) with sandwich s.e.
##   7. goodness-of-fit tables and figures
################################################################################

suppressPackageStartupMessages({
  library(VineCopula)   # BiCopCDF / BiCopPDF / BiCopHfunc / BiCopTau2Par
  library(numDeriv)     # jacobian / hessian for the sandwich
  library(MASS)         # polr(): starting values of the ordinal deductible model
})
has_kable <- requireNamespace("knitr", quietly = TRUE) && requireNamespace("kableExtra", quietly = TRUE)
if (has_kable) { library(knitr); library(kableExtra); options(knitr.table.format = "latex") }

## ---- 0a. settings ------------------------------------------------------------
cfg <- list(
  years_in      = 2006:2009,
  ded_levels    = c(500, 1000, 2500, 5000, 10000, 15000, 25000, 50000, 75000, 100000),
  kmax          = 16,          # extra terms kept in the k-sum of eq. (10) (as in the original code)
  latent_mu     = "groundup",  # "groundup": F_{N*} in eq. (10) has mean exp(x'beta_N)      [paper, eqs. (10)-(11)]
                               # "thinned" : F_{N*} has mean p(d;x) exp(x'beta_N)           [as coded in the original scripts]
  run_selection = TRUE,        # run the 9-family copula selection (slow); FALSE = load cached results
  sel_families  = NULL,        # NULL = all nine families; e.g. c(1, 2, 5) to restrict (family codes of VineCopula)
  sel_dedtypes  = c("cont", "disc"),   # deductible types to run the selection for (run each in its own R process to parallelise)
  stop_after_selection = FALSE,        # TRUE: save the selection cache and quit (used by the parallel drivers)
  auto_select   = TRUE,        # TRUE: each joint model uses its own AIC-best copula family
                               # FALSE: use the families fixed in paper_families below
  paper_families = list(SCM1 = 13, SCM2 = 13, SCM3 = 13, SDM1 = 2, SDM2 = 2, SDM3 = 2, FCM1 = 2, FCM2 = 2, FDM1 = 1, FDM2 = 1),
  sandwich      = TRUE,        # compute Godambe sandwich s.e. for the final models
  hess_r        = 2,           # Richardson order for numDeriv::hessian (2 = faster, 4 = more accurate)
  hess_d        = 0.01,        # relative step of numDeriv::hessian (0.01 keeps the ordinal thresholds of the deductible model ordered)
  nmax_gof      = 9,           # frequency goodness-of-fit tables: categories 0,1,...,nmax-1, >= nmax
  cache_file    = "copula_selection_cache.RData",   # a vector of files is merged when loading, e.g. c("sel_cont.RData", "sel_disc.RData")
  sw_preload    = NULL,        # optional RData file(s) with a list `sw` of precomputed sandwich objects (chains found there are not recomputed)
  out_file      = "Empirical_Revised_Results.RData",
  data_dir      = ".",                   # directory containing poldat.RData and claimdat.RData
  output_dir    = "empirical_merged_outputs",
  holdout_years = integer(0),           # e.g. 2009 => train 2006-2008, evaluate 2009
  coverage_to_dollars = 1e6,            # Coverage is in MILLIONS of dollars; set 1 if already in USD
  premium_gl_n  = 48L,                  # Gauss-Legendre quadrature for coverage-limited stop-loss
  capped_quad_check = TRUE,             # 48- vs 96-node check on a stratified subsample
  capped_quad_tol = 0.02,
  unlimited_gl_n = 48L,                 # unconstrained mean checked again with twice the nodes
  unlimited_rel_tol = 0.03,
  compute_unlimited = TRUE,
  count_tail_eps = 1e-9,
  count_max_k  = 600L,
  include_independence = TRUE,          # baseline for AIC family selection
  run_severity_family_table = TRUE,    # optional original GB2 family comparison
  run_qq       = TRUE,                  # optional, appendix diagnostics only
  run_gof      = TRUE,                 # optional, appendix count GOF
  # Empirical manuscript reporting (Tables 1-3, Figures 1-2, Appendix)
  reporting_enabled = TRUE,
  report_risk_groups = 10L,
  report_bootstrap_reps = 199L,
  report_bootstrap_seed = 20261110L,
  report_min_bootstrap_clusters = 15L,
  report_op_plot_cap = 3,
  report_common_copula_sensitivity = TRUE,
  report_common_family_codes = c(1L, 2L), # Gaussian and Student-t
  selection_require_convergence = TRUE,
  # Same frozen training estimates for both 2010 evaluations.
  run_oos_full = TRUE,
  run_oos_exclude_max = TRUE,
  censored_init_from_complete = TRUE,
  plot_width   = 11,
  plot_height  = 7
)
if (exists("cfg_override")) cfg[names(cfg_override)] <- cfg_override
if (cfg$latent_mu != "groundup")
  stop("Pure-premium prediction requires latent_mu = 'groundup'; 'thinned' has a different latent-count interpretation")
if (!setequal(cfg$sel_dedtypes, c("cont", "disc")))
  stop("This comparative analysis requires both cont and disc deductible selection")
stopifnot(length(cfg$coverage_to_dollars) == 1L,
          is.numeric(cfg$coverage_to_dollars),
          is.finite(cfg$coverage_to_dollars), cfg$coverage_to_dollars > 0)
stopifnot(all(cfg$holdout_years %in% cfg$years_in))
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- function(name) file.path(cfg$output_dir, name)
cfg$cache_file <- output_file(if (length(cfg$holdout_years))
  paste0("selection_train_excluding_", paste(cfg$holdout_years, collapse = "_"), ".RData")
  else "selection_train_all_years_init_from_complete_v2.RData")
cfg$out_file <- output_file(cfg$out_file)
write.csv(data.frame(option = names(cfg), value = vapply(cfg, function(v)
  paste(as.character(v), collapse = ","), character(1))),
  output_file("run_config.csv"), row.names = FALSE)


## ---- 0b. data ----------------------------------------------------------------
load(file.path(cfg$data_dir, "poldat.RData"))          # policy-year data:   PolicyNum, Year, Type*, Coverage, Deduct, Freq
load(file.path(cfg$data_dir, "claimdat.RData"))        # claim data:         PolicyNum, Year, Type*, Coverage, Deduct, LossBeforeDeductible

pol_all <- poldat[poldat$Year %in% cfg$years_in, , drop = FALSE]
clm_all <- claimdat[claimdat$Year %in% cfg$years_in, , drop = FALSE]
if (length(cfg$holdout_years)) {
  pol <- pol_all[!pol_all$Year %in% cfg$holdout_years, , drop = FALSE]
  clm <- clm_all[!clm_all$Year %in% cfg$holdout_years, , drop = FALSE]
  pol_evaluate_raw <- pol_all[pol_all$Year %in% cfg$holdout_years, , drop = FALSE]
  clm_evaluate_raw <- clm_all[clm_all$Year %in% cfg$holdout_years, , drop = FALSE]
} else {
  pol <- pol_all
  clm <- clm_all
  pol_evaluate_raw <- pol_all
  clm_evaluate_raw <- clm_all
}
if (!nrow(pol) || !nrow(pol_evaluate_raw) || !nrow(clm))
  stop("Empty training or evaluation sample; check years and input data")
if (any(!is.finite(pol$Coverage) | pol$Coverage <= 0)) stop("Invalid training coverage")
if (any(!is.finite(pol_evaluate_raw$Coverage) | pol_evaluate_raw$Coverage <= 0)) stop("Invalid evaluation coverage")
if (any(!is.finite(clm$LossBeforeDeductible) | clm$LossBeforeDeductible <= 0))
  stop("GB2 requires strictly positive ground-up losses in the training sample")
cat("Training policy-years:", nrow(pol), " Evaluation policy-years:",
    nrow(pol_evaluate_raw), "\n")

pol$discDed  <- match(pol$Deduct, cfg$ded_levels)      # ordinal deductible level 1..K
clm$discDed  <- match(clm$Deduct, cfg$ded_levels)
stopifnot(!anyNA(pol$discDed), !anyNA(clm$discDed))
pol$dedratio <- pol$Deduct / pol$Coverage               # continuous deductible (ratio)
clm$dedratio <- clm$Deduct / clm$Coverage

clm$LossAfterDeductible <- pmax(0, clm$LossBeforeDeductible - clm$Deduct)   # censored payment
clm$posFreq             <- as.integer(clm$LossBeforeDeductible > clm$Deduct)

## number of payments per policy-year (truncated frequency)
pay <- aggregate(posFreq ~ PolicyNum + Year, data = clm, sum)
pol <- merge(pol, pay, all.x = TRUE)
pol$posFreq[is.na(pol$posFreq)] <- 0

## the independent unit of the likelihood (7) is the policy-year; claims are mapped to their policy-year
pol$unit <- seq_len(nrow(pol))
clm$unit <- match(paste(clm$PolicyNum, clm$Year), paste(pol$PolicyNum, pol$Year))
if (anyNA(clm$unit)) {                       # claims without a policy record: give them their own unit
  warning(sum(is.na(clm$unit)), " claims have no matching policy-year record; treated as separate units")
  clm$unit[is.na(clm$unit)] <- nrow(pol) + seq_len(sum(is.na(clm$unit)))
}
clmpos  <- clm[clm$LossAfterDeductible > 0, ]           # truncated claim data
m_units <- max(nrow(pol), clm$unit)

## design matrices --------------------------------------------------------------
fml   <- ~ TypeCity + TypeCounty + TypeSchool + TypeTown + TypeVillage + log(Coverage)
X_pol <- model.matrix(fml, pol); X_clm <- model.matrix(fml, clm); X_pos <- model.matrix(fml, clmpos)
XD    <- function(X, dedtype) {                    # covariates of the deductible model
  if (dedtype == "cont") X[, colnames(X) != "log(Coverage)", drop = FALSE]   # ratio absorbs the coverage
  else                   X[, -1, drop = FALSE]                                # ordinal logit: no intercept
}
K_ded <- length(cfg$ded_levels)

dat_pol <- list(X = X_pol, n_cen = pol$Freq, n_tru = pol$posFreq, d = pol$Deduct,
                dratio = pol$dedratio, dlev = pol$discDed, unit = pol$unit)
dat_clm <- list(X = X_clm, y = clm$LossAfterDeductible, ystar = clm$LossBeforeDeductible,
                d = clm$Deduct, dratio = clm$dedratio, dlev = clm$discDed, unit = clm$unit)
dat_pos <- list(X = X_pos, y = clmpos$LossAfterDeductible, ystar = clmpos$LossBeforeDeductible,
                d = clmpos$Deduct, dratio = clmpos$dedratio, dlev = clmpos$discDed, unit = clmpos$unit)

################################################################################
## 1. distributions and copula helpers
################################################################################
log1pexp <- function(x) ifelse(x > 30, x, log1p(exp(x)))

## GB2 in the (mu, sigma, alpha1, alpha2) parametrisation of the paper (Appendix)
dGB2 <- function(y, mu, sigma, a1, a2, log = FALSE) {
  z  <- (log(y) - mu) / sigma
  ld <- a1 * z - log(y) - log(abs(sigma)) - lbeta(a1, a2) - (a1 + a2) * log1pexp(z)
  if (log) ld else exp(ld)
}
pGB2 <- function(y, mu, sigma, a1, a2) pf((a2 / a1) * exp((log(y) - mu) / sigma), 2 * a1, 2 * a2)

## generalized gamma (flexsurv's "orig" parametrisation, used for the severity family comparison)
dgengamma_orig <- function(x, shape, scale, k, log = FALSE) {
  ld <- log(shape) - lgamma(k) - log(scale) + (shape * k - 1) * log(x / scale) - (x / scale)^shape
  if (log) ld else exp(ld)
}
pgengamma_orig <- function(x, shape, scale, k) pgamma((x / scale)^shape, shape = k)

## continuous deductible model: GB2 regression for the ratio D/Coverage
D_cont_cdf    <- function(theta, X, r) { k <- ncol(X); pGB2(r, drop(X %*% theta[1:k]), theta[k + 1], theta[k + 2], theta[k + 3]) }
D_cont_logpdf <- function(theta, X, r) { k <- ncol(X); dGB2(r, drop(X %*% theta[1:k]), theta[k + 1], theta[k + 2], theta[k + 3], log = TRUE) }

## discrete deductible model: proportional-odds logit  P(D <= j | x) = logistic(a_j + x'b), j = 1..K-1
D_disc_cdf <- function(theta, X, j) {
  a <- theta[1:(K_ded - 1)]; b <- theta[K_ded:(K_ded - 1 + ncol(X))]
  eta <- drop(X %*% b); jj <- pmin(pmax(j, 1), K_ded - 1)
  out <- plogis(a[jj] + eta); out[j <= 0] <- 0; out[j >= K_ded] <- 1
  out
}
D_disc_pmf <- function(theta, X, j) D_disc_cdf(theta, X, j) - D_disc_cdf(theta, X, j - 1)

## -- fast t-copula h-function: VineCopula's C routine costs ~5 microseconds per point, which makes the discrete
##    deductible models (2e5 - 4e6 points per likelihood evaluation) prohibitively slow. qt(., nu) and pt(., nu)
##    are replaced by cubic splines on fixed grids (logit / asinh scale), built once per nu; max abs error ~1e-8.
.t_cache <- new.env()
qt_fun <- function(nu) { key <- paste0("q", format(nu, digits = 15))
  if (is.null(f <- .t_cache[[key]])) { g <- seq(-28, 28, by = 0.02); f <- splinefun(g, qt(plogis(g), df = nu), method = "natural"); assign(key, f, envir = .t_cache) }
  function(p) f(qlogis(p)) }
pt_fun <- function(nu) { key <- paste0("p", format(nu, digits = 15))
  if (is.null(f <- .t_cache[[key]])) { a <- seq(-20, 20, by = 0.01); f <- splinefun(a, pt(sinh(a), df = nu), method = "natural"); assign(key, f, envir = .t_cache) }
  function(x) pmin(pmax(f(asinh(x)), 0), 1) }
h1_t <- function(u, v, rho, nu) {          # P(V <= v | U = u) for the t copula (exchangeable: h2(u, v) = h1_t(v, u))
  QT <- qt_fun(nu); PT <- pt_fun(nu + 1); x <- QT(u); y <- QT(v)
  PT((y - rho * x) / sqrt((nu + x^2) * (1 - rho^2) / (nu + 1)))
}

## bivariate copula pieces (VineCopula family codes; 0 = independence)
clamp01 <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)
tau2par <- function(fam, tau) if (fam == 0) 0 else BiCopTau2Par(fam, tau)
cop <- function(what, u, v, fam, par, par2 = 0) {
  ## what = "C": C(u,v)   "c": c(u,v)   "h1": dC/du = P(V<=v|U=u)   "h2": dC/dv = P(U<=u|V=v)
  if (!length(u)) return(numeric(0))
  if (fam == 0) return(switch(what, C = u * v, c = rep(1, length(u)), h1 = v, h2 = u))
  u <- as.vector(clamp01(u)); v <- as.vector(clamp01(v))
  switch(what,
         C  = if (fam == 1) pbvnorm(qnorm(u), qnorm(v), par) else BiCopCDF(u, v, fam, par, par2),
         c  = BiCopPDF(u, v, fam, par, par2),
         h1 = if (fam == 2) h1_t(u, v, par, par2) else BiCopHfunc(u, v, fam, par, par2)$hfunc1,
         h2 = if (fam == 2) h1_t(v, u, par, par2) else BiCopHfunc(u, v, fam, par, par2)$hfunc2)
}
cop_families <- data.frame(
  name    = c("Gaussian", "Student's t", "Clayton", "Gumbel", "Frank", "Joe",
              "Survival Clayton", "Survival Gumbel", "Survival Joe"),
  fam     = c(1, 2, 3, 4, 5, 6, 13, 14, 16),
  tau_lo  = c(-0.999, -0.999, 1e-3, 1e-3, -0.999, 1e-3, 1e-3, 1e-3, 1e-3),
  tau_hi  = rep(0.999, 9),
  two_par = c(FALSE, TRUE, rep(FALSE, 7)), stringsAsFactors = FALSE)
fam_name <- function(fam) cop_families$name[match(fam, cop_families$fam)]

## -- exact, vectorised bivariate normal cdf P(X <= h, Y <= k), Genz (2004) / Drezner-Wesolowsky:
##    6 / 12 / 20 Gauss-Legendre points for |rho| < 0.925 and the Taylor-corrected 20-point rule otherwise
gl_nodes <- function(n) { i <- 1:(n - 1); b <- i / sqrt(4 * i^2 - 1); A <- matrix(0, n, n)
  A[cbind(i, i + 1)] <- b; A[cbind(i + 1, i)] <- b; e <- eigen(A, symmetric = TRUE); list(x = e$values, w = 2 * e$vectors[1, ]^2) }
GLQ <- list(gl_nodes(6), gl_nodes(12), gl_nodes(20))
pbvnorm <- function(h, k, rho) {
  n <- max(length(h), length(k)); h <- rep_len(h, n); k <- rep_len(k, n)
  if (abs(rho) < 0.925) {
    G <- GLQ[[1 + (abs(rho) >= 0.3) + (abs(rho) >= 0.75)]]; ng <- length(G$x)
    a <- asin(rho); th <- a * (G$x + 1) / 2; w <- G$w * a / 2
    S <- outer(h^2 + k^2, rep(1, ng)) - 2 * outer(h * k, sin(th))
    return(pnorm(h) * pnorm(k) + drop(exp(-S / (2 * matrix(cos(th)^2, n, ng, byrow = TRUE))) %*% w) / (2 * pi))
  }
  ## |rho| >= 0.925: BVNU(-h, -k, rho) = P(X > -h, Y > -k) = P(X <= h, Y <= k)
  H <- -h; K <- -k; r <- rho
  if (r < 0) K <- -K
  HK <- H * K; bvn <- rep(0, n)
  if (abs(r) < 1) {
    AS <- (1 - r) * (1 + r); A <- sqrt(AS); BS <- (H - K)^2; C <- (4 - HK) / 8; D <- (12 - HK) / 16
    ASR <- -(BS / AS + HK) / 2
    ok <- ASR > -100
    bvn[ok] <- A * exp(ASR[ok]) * (1 - C[ok] * (BS[ok] - AS) * (1 - D[ok] * BS[ok] / 5) / 3 + C[ok] * D[ok] * AS^2 / 5)
    ok2 <- -HK < 100; B <- sqrt(BS)
    bvn[ok2] <- bvn[ok2] - exp(-HK[ok2] / 2) * sqrt(2 * pi) * pnorm(-B[ok2] / A) * B[ok2] * (1 - C[ok2] * BS[ok2] * (1 - D[ok2] * BS[ok2] / 5) / 3)
    A2 <- A / 2; G <- GLQ[[3]]
    for (i in seq_along(G$x)) {
      XS <- (A2 * (G$x[i] + 1))^2; RS <- sqrt(1 - XS)
      ASR <- -(BS / XS + HK) / 2; ok <- ASR > -100
      bvn[ok] <- bvn[ok] + A2 * G$w[i] * exp(ASR[ok]) * (exp(-HK[ok] * (1 - RS) / (2 * (1 + RS))) / RS - (1 + C[ok] * XS * (1 + D[ok] * XS)))
    }
    bvn <- -bvn / (2 * pi)
  }
  if (r > 0) {
    bvn <- bvn + pnorm(-pmax(H, K))
  } else {
    bvn <- -bvn
    i <- K > H
    bvn[i] <- bvn[i] + ifelse(H[i] < 0, pnorm(K[i]) - pnorm(H[i]), pnorm(-H[i]) - pnorm(-K[i]))
  }
  pmin(pmax(bvn, 0), 1)
}
## -- tanh-sinh (double exponential) nodes on [-1, 1]: used for the t copula, whose cdf in VineCopula
##    is evaluated row by row AND rounds the degrees of freedom to an integer (BiCopCDF documentation)
ts_nodes <- function(h = 0.15, tmax = 3.3) { t <- seq(-tmax, tmax, by = h); x <- tanh(pi / 2 * sinh(t))
  w <- h * pi / 2 * cosh(t) / cosh(pi / 2 * sinh(t))^2; keep <- abs(x) < 1 - 1e-15; list(x = x[keep], w = w[keep]) }
TSQ <- ts_nodes(h = 0.25, tmax = 3.25)      # 25 nodes: absolute error ~1e-8 for the t copula (41 nodes: 1e-12)
## -- P(um < U <= up, V <= v) = C(up, v) - C(um, v) = int_{um}^{up} h1(s, v) ds : the probability of the atom
##    (F_D(d-1), F_D(d)] of a discrete deductible jointly with {V <= v}; exact for every family, no cancellation
cop_rect <- function(um, up, v, fam, par, par2 = 0) {
  if (!length(v)) return(numeric(0))
  if (fam == 0) return((up - um) * v)
  um <- clamp01(um); up <- clamp01(up); v <- clamp01(v)
  if (fam == 1) return(pbvnorm(qnorm(up), qnorm(v), par) - pbvnorm(qnorm(um), qnorm(v), par))
  if (fam != 2) return(BiCopCDF(up, v, fam, par, par2) - BiCopCDF(um, v, fam, par, par2))
  ## t copula: tanh-sinh quadrature on the atom; the grid is shared by all observations with the same atom
  ## (e.g. the 17 terms of the k-sum of a policy), and the h-function uses the spline-based qt / pt
  nn <- length(TSQ$x); key <- paste(um, up); first <- !duplicated(key); uidx <- match(key, key[first])
  ru <- um[first]; w <- up[first] - ru
  S  <- matrix(ru, length(ru), nn) + outer(w, (TSQ$x + 1) / 2)
  QT <- qt_fun(par2); PT <- pt_fun(par2 + 1)
  X  <- matrix(QT(clamp01(S)), nrow(S))[uidx, , drop = FALSE]; y <- QT(v)
  H  <- matrix(PT((y - par * X) / sqrt((par2 + X^2) * (1 - par^2) / (par2 + 1))), nrow(X))
  drop(H %*% TSQ$w) * w[uidx] / 2
}

################################################################################
## 2. per-unit log-likelihood contributions
##    Every function returns a VECTOR of contributions (one per policy-year or
##    per claim) and takes the earlier-stage parameters through `prev`
##    (prev$D = deductible model, prev$S = severity model).
################################################################################

## -- conditional severity pieces:  log f_{Y*|D}(y*|d) at y* and log F_{Y*|D}(d|d) (only on rows `need`)
sev_pieces <- function(theta, prev, dedtype, fam, dat, ystar, need = seq_along(ystar)) {
  X <- dat$X; k <- ncol(X)
  beta <- theta[1:k]; sigma <- theta[k + 1]; a1 <- theta[k + 2]; a2 <- theta[k + 3]
  par <- par2 <- 0
  if (dedtype != "none") { par <- tau2par(fam, theta[k + 4]); if (fam == 2) par2 <- theta[k + 5] }
  mu  <- drop(X %*% beta)
  lf  <- dGB2(ystar, mu, sigma, a1, a2, log = TRUE)             # log f_{Y*}(y*)
  FY  <- pGB2(ystar, mu, sigma, a1, a2)                         # F_{Y*}(y*)
  FYd <- pGB2(dat$d[need], mu[need], sigma, a1, a2)             # F_{Y*}(d)
  lFdc <- rep(NA_real_, length(ystar))
  if (dedtype == "none") {                                      # exogenous: f_{Y*|D} = f_{Y*}
    lfc <- lf; lFdc[need] <- log(FYd)
  } else if (dedtype == "cont") {                               # eq. (4)
    u   <- D_cont_cdf(prev$D, XD(X, "cont"), dat$dratio)
    lfc <- lf + log(cop("c", u, FY, fam, par, par2))
    lFdc[need] <- log(cop("h1", u[need], FYd, fam, par, par2))
  } else {                                                      # eq. (6)
    XDd <- XD(X, "disc")
    up  <- D_disc_cdf(prev$D, XDd, dat$dlev); um <- D_disc_cdf(prev$D, XDd, dat$dlev - 1)
    lfD <- log(up - um)
    lfc <- lf + log(pmax(cop("h2", up, FY, fam, par, par2) - cop("h2", um, FY, fam, par, par2), 1e-300)) - lfD
    lFdc[need] <- log(pmax(cop_rect(um[need], up[need], FYd, fam, par, par2), 1e-8)) - lfD[need]
  }
  list(lfc = lfc, lFdc = lFdc)
}

## -- severity: one function for SEM1-3 / SCM1-3 / SDM1-3   (eqs. (8), (9), (10), severity parts)
sev_loglik_i <- function(theta, prev, scenario, dedtype, fam, dat) {
  ystar <- if (scenario == "complete") dat$ystar else dat$y + dat$d
  need  <- switch(scenario, complete = integer(0), censored = which(dat$y == 0), truncated = seq_along(ystar))
  pc    <- sev_pieces(theta, prev, dedtype, fam, dat, ystar, need)
  switch(scenario,
         complete  = pc$lfc,
         censored  = ifelse(dat$y > 0, pc$lfc, pc$lFdc),
         truncated = pc$lfc - log1p(-pmin(exp(pc$lFdc), 1 - 1e-12)))
}

## -- retention probability p(d;x) = 1 - F_{Y*|D}(d|d;x) at the policy level, from a fitted severity model.
##    Memoised: during the frequency-stage optimisation (and for the theta_N block of the Hessian) the severity and
##    deductible parameters do not change, so p(d;x) is computed once and re-used.
.p_cache <- new.env()
retention_prob <- function(theta_S, theta_D, sev_spec, dat) {
  key <- paste(c(format(c(theta_S, theta_D), digits = 17), sev_spec$dedtype, sev_spec$fam, length(dat$d), sum(dat$d)), collapse = "|")
  if (!is.null(p <- .p_cache[[key]])) return(p)
  need <- seq_along(dat$d)
  pc <- sev_pieces(theta_S, list(D = theta_D), sev_spec$dedtype, sev_spec$fam, dat, ystar = dat$d, need = need)
  p <- pmin(pmax(1 - exp(pc$lFdc), 1e-10), 1)
  if (length(ls(.p_cache)) > 200) rm(list = ls(.p_cache), envir = .p_cache)
  assign(key, p, envir = .p_cache); p
}

## -- log g(n|d;x) under truncated reporting with a non-trivial copula: the k-sum of eq. (10)
##    F_{N*} and the copula are evaluated once on the grid k = n-1, n, ..., n+kmax ((kmax+2) values per policy)
##    and differenced, which halves the number of copula evaluations compared with evaluating F(k) and F(k-1) separately
trunc_freq_logpmf <- function(n, r, mu, p, dedtype, u = NULL, up = NULL, um = NULL, fam, par, par2, kmax) {
  K <- kmax + 1; m <- length(n)
  kmat  <- matrix(0:kmax, K, m) + matrix(n, K, m, byrow = TRUE)          # k = n, n+1, ..., n+kmax
  Kg    <- rbind(kmat[1, ] - 1, kmat)                                     # grid k = n-1, ..., n+kmax   ((K+1) x m)
  FNg   <- clamp01(pnbinom(Kg, size = r, mu = matrix(mu, K + 1, m, byrow = TRUE)), 1e-15)   # F_{N*}(k) on the grid
  if (dedtype == "cont") {                                              # eq. (3)
    Hg  <- matrix(cop("h1", matrix(u, K + 1, m, byrow = TRUE), FNg, fam, par, par2), K + 1, m)
    fNk <- Hg[-1, , drop = FALSE] - Hg[-(K + 1), , drop = FALSE]
  } else {                                                              # eq. (5)
    Rg  <- matrix(cop_rect(rep(um, each = K + 1), rep(up, each = K + 1), as.vector(FNg), fam, par, par2), K + 1, m)
    fNk <- (Rg[-1, , drop = FALSE] - Rg[-(K + 1), , drop = FALSE]) / matrix(up - um, K, m, byrow = TRUE)
  }
  w <- dbinom(matrix(n, K, m, byrow = TRUE), size = kmat, prob = matrix(p, K, m, byrow = TRUE))
  g <- colSums(fNk * w)
  log(pmax(g, 1e-8))
}

## -- frequency: one function for FEM1-2 / FCM1-2 / FDM1-2  (eqs. (8), (10)-(11), frequency parts)
freq_loglik_i <- function(theta, prev, scenario, dedtype, fam, sev_spec, dat) {
  X <- dat$X; k <- ncol(X)
  beta <- theta[1:k]; r <- theta[k + 1]
  par <- par2 <- 0
  if (dedtype != "none") { par <- tau2par(fam, theta[k + 2]); if (fam == 2) par2 <- theta[k + 3] }
  mu <- drop(exp(X %*% beta))
  u <- up <- um <- NULL
  if (dedtype == "cont") u <- D_cont_cdf(prev$D, XD(X, "cont"), dat$dratio)
  if (dedtype == "disc") { XDd <- XD(X, "disc"); up <- D_disc_cdf(prev$D, XDd, dat$dlev); um <- D_disc_cdf(prev$D, XDd, dat$dlev - 1) }

  if (scenario == "censored") {                                      # full / censored reporting: N = N*
    n <- dat$n_cen
    if (dedtype == "none") return(dnbinom(n, size = r, mu = mu, log = TRUE))
    FNp <- pnbinom(n, size = r, mu = mu); FNm <- pnbinom(n - 1, size = r, mu = mu)
    if (dedtype == "cont")
      return(log(pmax(cop("h1", u, FNp, fam, par, par2) - cop("h1", u, FNm, fam, par, par2), 1e-300)))
    return(log(pmax(cop_rect(um, up, FNp, fam, par, par2) - cop_rect(um, up, FNm, fam, par, par2), 1e-300)) - log(up - um))
  }
  ## truncated reporting: N = number of payments; p(d;x) recomputed from the severity parameters in prev$S
  n <- dat$n_tru
  p <- retention_prob(prev$S, prev$D, sev_spec, dat)
  if (dedtype == "none") return(dnbinom(n, size = r, mu = mu * p, log = TRUE))     # eq. (11)
  mu_lat <- if (cfg$latent_mu == "groundup") mu else mu * p
  trunc_freq_logpmf(n, r, mu_lat, p, dedtype, u, up, um, fam, par, par2, cfg$kmax)
}

################################################################################
## 3. stages, estimation and the Godambe sandwich
################################################################################
## A "stage" bundles a contribution function with its data, bounds and names.
make_stage <- function(role, name, contrib, unit, lower, upper, par_names, spec = list())
  list(role = role, name = name, contrib = contrib, unit = unit, lower = lower, upper = upper,
       par_names = par_names, spec = spec)

stage_D <- function(dedtype) {
  X <- XD(X_pol, dedtype)
  if (dedtype == "cont")
    make_stage("D", "D_cont", function(th, prev) D_cont_logpdf(th, X, pol$dedratio), pol$unit,
               lower = c(rep(-Inf, ncol(X)), 1e-3, 1e-3, 1e-3), upper = rep(Inf, ncol(X) + 3),
               par_names = c(colnames(X), "sigma", "alpha1", "alpha2"), spec = list(dedtype = "cont"))
  else
    make_stage("D", "D_disc", function(th, prev) log(D_disc_pmf(th, X, pol$discDed)), pol$unit,
               lower = rep(-Inf, K_ded - 1 + ncol(X)), upper = rep(Inf, K_ded - 1 + ncol(X)),
               par_names = c(paste0(1:(K_ded - 1), "|", 2:K_ded), colnames(X)), spec = list(dedtype = "disc"))
}
cop_bounds <- function(fam) {                      # (lower, upper) for (tau[, nu])
  i <- match(fam, cop_families$fam)
  if (fam == 0) return(list(lo = NULL, hi = NULL, nm = NULL))
  lo <- cop_families$tau_lo[i]; hi <- cop_families$tau_hi[i]; nm <- "tau"
  if (cop_families$two_par[i]) { lo <- c(lo, 2.01); hi <- c(hi, 100); nm <- c(nm, "nu") }
  list(lo = lo, hi = hi, nm = nm)
}
stage_sev <- function(name, scenario, dedtype, fam = 0) {
  dat <- if (scenario == "truncated") dat_pos else dat_clm
  cb  <- cop_bounds(fam); k <- ncol(dat$X)
  make_stage("S", name, function(th, prev) sev_loglik_i(th, prev, scenario, dedtype, fam, dat), dat$unit,
             lower = c(rep(-Inf, k), 1e-3, 1e-3, 1e-3, cb$lo), upper = c(rep(Inf, k), Inf, Inf, Inf, cb$hi),
             par_names = c(colnames(dat$X), "sigma", "alpha1", "alpha2", cb$nm),
             spec = list(scenario = scenario, dedtype = dedtype, fam = fam))
}
stage_freq <- function(name, scenario, dedtype, fam = 0, sev_spec = NULL) {
  cb <- cop_bounds(fam); k <- ncol(dat_pol$X)
  make_stage("N", name, function(th, prev) freq_loglik_i(th, prev, scenario, dedtype, fam, sev_spec, dat_pol), dat_pol$unit,
             lower = c(rep(-Inf, k), 1e-3, cb$lo), upper = c(rep(Inf, k), Inf, cb$hi),
             par_names = c(colnames(dat_pol$X), "r", cb$nm),
             spec = list(scenario = scenario, dedtype = dedtype, fam = fam, sev_spec = sev_spec))
}

## -- minimiser used everywhere: L-BFGS-B followed by an nlminb polish (keeps the better of the two)
minimize <- function(nll, init, lower, upper, control = list(maxit = 3000)) {
  f <- function(th) { v <- nll(th); if (is.finite(v)) v else 1e10 }
  init <- pmin(pmax(init, lower + 1e-6), upper - 1e-6)
  o1 <- optim(init, f, method = "L-BFGS-B", lower = lower, upper = upper, control = control)
  o2 <- tryCatch(nlminb(o1$par, f, lower = lower, upper = upper, control = list(eval.max = 400, iter.max = 100)), error = function(e) NULL)
  ## convergence = 0 if either optimiser reports convergence (nlminb often reports "false convergence" when it
  ## is started at the L-BFGS-B optimum and can only improve the objective by a negligible amount)
  if (!is.null(o2) && is.finite(o2$objective) && o2$objective < o1$value)
    list(par = o2$par, value = o2$objective, convergence = as.integer(o1$convergence != 0 && o2$convergence != 0))
  else list(par = o1$par, value = o1$value, convergence = o1$convergence)
}
## -- maximum likelihood for one stage (earlier-stage estimates fixed in `prev`)
fit_stage <- function(stage, init, prev = list(), naive_hessian = TRUE) {
  nll  <- function(th) -sum(stage$contrib(th, prev))
  best <- minimize(nll, init, stage$lower, stage$upper)
  stage$par  <- setNames(best$par, stage$par_names); stage$prev <- prev
  stage$nll  <- best$value; stage$loglik <- -best$value; stage$convergence <- best$convergence
  stage$n    <- length(stage$contrib(stage$par, prev))
  stage$AIC  <- 2 * best$value + 2 * length(best$par)
  stage$BIC  <- 2 * best$value + log(stage$n) * length(best$par)
  if (naive_hessian) {                      # stage-wise Hessian (earlier stages treated as known)
    stage$hess     <- tryCatch(numDeriv::hessian(nll, stage$par, method.args = list(d = cfg$hess_d, r = cfg$hess_r)), error = function(e) NULL)
    stage$se_naive <- tryCatch(sqrt(pmax(diag(solve(stage$hess)), 0)), error = function(e) rep(NA_real_, length(best$par)))
  }
  stage
}
## `prev` list for stage k of a chain, from the parameters of stages 1..k-1
prev_from <- function(chain, k, pars = lapply(chain, `[[`, "par")) {
  pv <- list(); for (j in seq_len(k - 1)) pv[[chain[[j]]$role]] <- pars[[j]]; pv
}
rowsum_full <- function(G, unit, m) { out <- matrix(0, m, ncol(G)); rs <- rowsum(G, unit); out[as.integer(rownames(rs)), ] <- rs; out }

## -- Godambe (sandwich) covariance of a sequentially estimated chain of stages
##    chain = list(stage_1, ..., stage_K) in estimation order, each fitted with fit_stage()
##    and with stage_k$prev equal to the parameters of stages 1..k-1.
##    Returns V = H^{-1} J H^{-T}, where
##      psi_i = (d l_{1i}/d theta_1, ..., d l_{Ki}/d theta_K)     stacked per-unit stage scores
##      J     = sum_i psi_i psi_i'                                (block-diagonal in expectation)
##      H     = -d(sum_i psi_i)/d theta'                          (block lower-triangular)
sandwich_se <- function(chain, r = cfg$hess_r, verbose = TRUE) {
  roles <- sapply(chain, `[[`, "role"); stopifnot(!anyDuplicated(roles))
  pars  <- lapply(chain, `[[`, "par"); P <- lengths(pars); Ptot <- sum(P)
  idx   <- split(seq_len(Ptot), rep(seq_along(P), P))
  theta <- unlist(pars); names(theta) <- unlist(lapply(chain, function(s) paste(s$name, s$par_names, sep = ":")))
  Psi <- matrix(0, m_units, Ptot); H <- matrix(0, Ptot, Ptot)
  for (k in seq_along(chain)) {
    s <- chain[[k]]; prev_k <- prev_from(chain, k, pars)
    if (verbose) cat(sprintf("  sandwich: stage %d (%s) - scores ...", k, s$name))
    G <- numDeriv::jacobian(function(th) s$contrib(th, prev_k), s$par)         # n_k x P_k per-unit scores
    Psi[, idx[[k]]] <- rowsum_full(G, s$unit, m_units)
    if (verbose) cat(" Hessian block row ...")
    th_all <- unlist(pars[1:k]); loc <- split(seq_along(th_all), rep(1:k, P[1:k]))
    nll_k <- function(th) {
      pv <- list(); for (j in seq_len(k - 1)) pv[[roles[j]]] <- th[loc[[j]]]
      v <- -sum(s$contrib(th[loc[[k]]], pv)); if (is.finite(v)) v else NA
    }
    Hk <- numDeriv::hessian(nll_k, th_all, method.args = list(d = cfg$hess_d, r = r))
    H[idx[[k]], unlist(idx[1:k])] <- Hk[loc[[k]], , drop = FALSE]
    if (verbose) cat(" done\n")
  }
  J    <- crossprod(Psi)
  Hinv <- solve(H)
  V    <- Hinv %*% J %*% t(Hinv)
  se_naive <- se_robust <- rep(NA_real_, Ptot)
  for (k in seq_along(chain)) {
    ii  <- idx[[k]]; Hkk_inv <- solve(H[ii, ii])
    se_naive[ii]  <- sqrt(pmax(diag(Hkk_inv), 0))
    se_robust[ii] <- sqrt(pmax(diag(Hkk_inv %*% J[ii, ii] %*% Hkk_inv), 0))
  }
  tab <- data.frame(model = rep(sapply(chain, `[[`, "name"), P), parameter = unlist(lapply(chain, `[[`, "par_names")),
                    estimate = as.numeric(theta), se_naive = se_naive, se_robust = se_robust,
                    se_sandwich = sqrt(pmax(diag(V), 0)), stringsAsFactors = FALSE)
  tab$ratio <- tab$se_sandwich / tab$se_naive
  list(table = tab, V = V, H = H, J = J, idx = idx, theta = theta)
}
## estimates and sandwich s.e. of the LAST stage of a chain (the model of interest)
last_stage_table <- function(sw) { k <- length(sw$idx); sw$table[sw$idx[[k]], ] }

## -- LaTeX / console table helper
tex_table <- function(x, caption, digits = 3, header = NULL, bottom = NULL) {
  if (!has_kable) { print(round(x, digits)); return(invisible(NULL)) }
  if (is.null(bottom)) bottom <- paste0("\\hhline{", strrep("=", ncol(x) + 1), "}")
  kb <- kable(x, digits = digits, caption = caption, booktabs = TRUE, escape = FALSE, bottomrule = bottom, linesep = "") %>%
    kable_styling(latex_options = c("hold_position", "scale_down"))
  if (!is.null(header)) kb <- add_header_above(kb, header)
  print(kb)
}
est_se_matrix <- function(...) {                     # cbind of (estimate, s.e.) columns from last_stage_table()s; rows aligned by parameter name
  L <- list(...); rn <- unique(unlist(lapply(L, `[[`, "parameter")))
  out <- do.call(cbind, lapply(L, function(t) { m <- matrix(NA_real_, length(rn), 2); i <- match(t$parameter, rn); m[i, 1] <- t$estimate; m[i, 2] <- t$se_sandwich; m }))
  rownames(out) <- rn; colnames(out) <- rep(c("Estimate", "s.e."), length(L)); out
}

################################################################################
## 4. D and SEM severity-only baseline (before any FEM)
################################################################################
cat("\n==== 1. Deductible models and marginal SEM1/2/3 ====\n")
fit_Dc <- fit_stage(stage_D("cont"),
  c(coef(lm(log(dedratio) ~ TypeCity + TypeCounty + TypeSchool + TypeTown + TypeVillage,
            data = pol)), 1, 1, 1))
po <- MASS::polr(factor(discDed, levels = 1:K_ded) ~
  TypeCity + TypeCounty + TypeSchool + TypeTown + TypeVillage + log(Coverage),
  data = pol, Hess = FALSE)
fit_Dd <- fit_stage(stage_D("disc"), c(po$zeta, -coef(po)))

init_sev <- function(dat) c(coef(lm(log(dat$ystar) ~ dat$X - 1)), 1, 1, 1)
SEM1 <- fit_stage(stage_sev("SEM1", "complete", "none"), init_sev(dat_clm))
SEM2 <- fit_stage(stage_sev("SEM2", "censored", "none"), SEM1$par)  # same INITIAL parameter vector as SEM1
SEM3 <- fit_stage(stage_sev("SEM3", "truncated", "none"), init_sev(dat_pos))

# Plot fitted GB2 mu = X'beta; identical X/Y ranges across all SIX entity panels.
# Put the legend in its OWN top strip, not on the red SEM1 line.

plot_sem_mu <- function(file = output_file("mu_vs_logcoverage_unified.pdf")) {
  
  d <- data.frame(
    TypeCity = 0,
    TypeCounty = 0,
    TypeSchool = 0,
    TypeTown = 0,
    TypeVillage = 0,
    Coverage = 1
  )
  
  types <- c("Misc", "City", "County", "School", "Town", "Village")
  
  lc_range <- range(log(pol$Coverage), finite = TRUE)
  lc_range <- lc_range + c(-1, 1) * 0.02 * diff(lc_range)
  
  lgrid <- seq(lc_range[1], lc_range[2], length.out = 120L)
  
  fits <- list(
    SEM1 = SEM1,
    SEM2 = SEM2,
    SEM3 = SEM3
  )
  
  cols <- c("firebrick3", "forestgreen", "royalblue3")
  
  lines_by_type <- vector("list", length(types))
  
  for (i in seq_along(types)) {
    
    nd <- d[rep(1L, length(lgrid)), , drop = FALSE]
    nd$Coverage <- exp(lgrid)
    
    if (i > 1L) {
      nd[[paste0("Type", types[i])]] <- 1
    }
    
    Xm <- model.matrix(fml, data = nd)
    
    lines_by_type[[i]] <- vapply(
      fits,
      function(f) {
        as.numeric(Xm %*% f$par[seq_len(ncol(Xm))])
      },
      numeric(length(lgrid))
    )
  }
  
  # Common Y-axis range across all six panels
  common_y <- range(unlist(lines_by_type), finite = TRUE)
  
  yrange <- max(diff(common_y), 0.1)
  
  # Reserve space above the fitted curves for each legend
  common_y <- common_y + c(
    -0.085 * yrange,
    0.25  * yrange
  )
  
  grDevices::pdf(
    file,
    width = cfg$plot_width,
    height = cfg$plot_height
  )
  
  on.exit(grDevices::dev.off(), add = TRUE)
  
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  
  # Six panels: 2 rows x 3 columns
  par(mfrow = c(2, 3))
  
  for (i in seq_along(types)) {
    
    par(
      mar = c(4, 4.3, 3.2, 0.7),
      mgp = c(2.3, 0.75, 0)
    )
    
    plot(
      lgrid,
      lines_by_type[[i]][, 1],
      type = "l",
      lwd = 2,
      col = cols[1],
      xlim = lc_range,
      ylim = common_y,
      main = types[i],
      xlab = "log(Coverage)",
      ylab = expression(mu)
    )
    
    lines(
      lgrid,
      lines_by_type[[i]][, 2],
      col = cols[2],
      lwd = 2
    )
    
    lines(
      lgrid,
      lines_by_type[[i]][, 3],
      col = cols[3],
      lwd = 2
    )
    
    # Individual legend at the top of each panel
    legend(
      "top",
      legend = names(fits),
      col = cols,
      lty = 1,
      lwd = 2,
      horiz = TRUE,
      bty = "n",
      cex = 1,
      inset = c(0, 0.01)
    )
    
    box()
  }
  
  invisible(NULL)
}

# Generate PDF
plot_sem_mu()

# Export log(Coverage) slopes
mu_slopes <- data.frame(
  Model = c("SEM1", "SEM2", "SEM3"),
  Slope_logCoverage = vapply(
    list(SEM1, SEM2, SEM3),
    function(f) unname(f$par["log(Coverage)"]),
    numeric(1)
  )
)

write.csv(
  mu_slopes,
  output_file("SEM_mu_logcoverage_slopes.csv"),
  row.names = FALSE
)

cat("SEM mu slopes (descriptive; NOT a formal endogeneity test):\n")
print(mu_slopes, row.names = FALSE)

# Preserve the full marginal regression table for manuscript Section 4.
sem_parameter_table <- do.call(rbind, lapply(list(SEM1, SEM2, SEM3), function(f) {
  data.frame(Model = f$name, Parameter = names(f$par),
    Estimate = unname(f$par),
    Stagewise_SE = if (!is.null(f$se_naive)) f$se_naive else NA_real_,
    AIC = f$AIC, Convergence = f$convergence)
}))
write.csv(sem_parameter_table, output_file("SEM_Parameter_Estimates.csv"), row.names = FALSE)
d_parameter_table <- do.call(rbind, lapply(list(fit_Dc, fit_Dd), function(f) {
  data.frame(Model = f$name, Parameter = names(f$par),
    Estimate = unname(f$par),
    Stagewise_SE = if (!is.null(f$se_naive)) f$se_naive else NA_real_,
    AIC = f$AIC, Convergence = f$convergence)
}))
write.csv(d_parameter_table, output_file("D_Parameter_Estimates.csv"), row.names = FALSE)
## 4d. severity family comparison (Appendix table: GB2 family and GG family) -----------
sev_dists <- list(   # extra parameters beyond the regression coefficients; densities / cdfs
  GB2     = list(np = 3, f = function(y, mu, t, log) dGB2(y, mu, t[1], t[2], t[3], log), F = function(y, mu, t) pGB2(y, mu, t[1], t[2], t[3])),
  GP      = list(np = 2, f = function(y, mu, t, log) dGB2(y, mu, 1, t[1], t[2], log),    F = function(y, mu, t) pGB2(y, mu, 1, t[1], t[2])),
  Burr    = list(np = 2, f = function(y, mu, t, log) dGB2(y, mu, t[1], 1, t[2], log),    F = function(y, mu, t) pGB2(y, mu, t[1], 1, t[2])),
  IBurr   = list(np = 2, f = function(y, mu, t, log) dGB2(y, mu, t[1], t[2], 1, log),    F = function(y, mu, t) pGB2(y, mu, t[1], t[2], 1)),
  Pareto  = list(np = 1, f = function(y, mu, t, log) dGB2(y, mu, 1, 1, t[1], log),       F = function(y, mu, t) pGB2(y, mu, 1, 1, t[1])),
  Loglogi = list(np = 1, f = function(y, mu, t, log) dGB2(y, mu, t[1], 1, 1, log),       F = function(y, mu, t) pGB2(y, mu, t[1], 1, 1)),
  GG      = list(np = 2, f = function(y, mu, t, log) dgengamma_orig(y, t[1], exp(mu), t[2], log), F = function(y, mu, t) pgengamma_orig(y, t[1], exp(mu), t[2])),
  Weibull = list(np = 1, f = function(y, mu, t, log) dgengamma_orig(y, t[1], exp(mu), 1, log),    F = function(y, mu, t) pgengamma_orig(y, t[1], exp(mu), 1)),
  Gamma   = list(np = 1, f = function(y, mu, t, log) dgamma(y, shape = t[1], scale = exp(mu), log = log), F = function(y, mu, t) pgamma(y, shape = t[1], scale = exp(mu))))
sev_family_table <- function() {
  out <- list()
  for (sc in c("complete", "censored", "truncated")) {
    dat <- if (sc == "truncated") dat_pos else dat_clm; k <- ncol(dat$X)
    ystar <- if (sc == "complete") dat$ystar else dat$y + dat$d
    res <- sapply(names(sev_dists), function(dn) {
      D <- sev_dists[[dn]]
      nll <- function(th) {
        mu <- drop(dat$X %*% th[1:k]); t <- th[-(1:k)]
        v <- switch(sc,
          complete  = D$f(ystar, mu, t, TRUE),
          censored  = ifelse(dat$y > 0, D$f(ystar, mu, t, TRUE), log(D$F(dat$d, mu, t))),
          truncated = D$f(ystar, mu, t, TRUE) - log1p(-pmin(D$F(dat$d, mu, t), 1 - 1e-12)))
        -sum(v)
      }
      o <- minimize(nll, c(coef(lm(log(ystar) ~ dat$X - 1)), rep(1, D$np)), lower = c(rep(-Inf, k), rep(1e-3, D$np)), upper = rep(Inf, k + D$np))
      c(Likelihood = -o$value, AIC = 2 * o$value + 2 * (k + D$np), BIC = 2 * o$value + log(length(ystar)) * (k + D$np))
    })
    out[[sc]] <- res
  }
  out
}
sev_family_tab <- sev_family_table()
for (sc in names(sev_family_tab)) tex_table(sev_family_tab[[sc]], paste("In-sample model selection criteria,", sc, "data"), digits = 1)

## Kendall's tau between F_D and the marginal CDF of the response, used as starting value
marg_cdf_sev <- function(fit, dat, scenario) {
  k <- ncol(dat$X); mu <- drop(dat$X %*% fit$par[1:k]); s <- fit$par[k + 1]; a1 <- fit$par[k + 2]; a2 <- fit$par[k + 3]
  if (scenario == "truncated") { Fd <- pGB2(dat$d, mu, s, a1, a2); (pGB2(dat$y + dat$d, mu, s, a1, a2) - Fd) / (1 - Fd) }
  else pGB2(if (scenario == "complete") dat$ystar else dat$y + dat$d, mu, s, a1, a2)
}
marg_cdf_freq <- function(fit, scenario, theta_S) {
  k <- ncol(X_pol); mu <- drop(exp(X_pol %*% fit$par[1:k]))
  if (scenario == "truncated") pnbinom(dat_pol$n_tru, size = fit$par[k + 1], mu = mu * retention_prob(theta_S, NULL, sev_none, dat_pol))
  else pnbinom(dat_pol$n_cen, size = fit$par[k + 1], mu = mu)
}
F_D_of <- function(dedtype, dat, theta_D)
  if (dedtype == "cont") D_cont_cdf(theta_D, XD(dat$X, "cont"), dat$dratio) else D_disc_cdf(theta_D, XD(dat$X, "disc"), dat$dlev)

## checkpointing: resume only when data and model settings match.
.ckpt <- new.env(); .ckpt$sel <- list(); .ckpt$partial <- list()
cache_signature <- list(
  policy_md5 = unname(tools::md5sum(file.path(cfg$data_dir, "poldat.RData"))),
  claim_md5 = unname(tools::md5sum(file.path(cfg$data_dir, "claimdat.RData"))),
  years = cfg$years_in, holdout = cfg$holdout_years, dedlevels = cfg$ded_levels,
  latent_mu = cfg$latent_mu, kmax = cfg$kmax,
  selection_families = cfg$sel_families,
  include_independence = cfg$include_independence,
  censored_init_from_complete = cfg$censored_init_from_complete,
  init_strategy_version = "censored_from_complete_selected_v2")
ckpt_load <- function() if (file.exists(cfg$cache_file[1])) {
  e <- new.env(parent = emptyenv()); load(cfg$cache_file[1], envir = e)
  if (!identical(e$cache_signature, cache_signature))
    stop("Copula cache is inconsistent with data/settings; delete the cache or use a separate output_dir")
  .ckpt$sel <- if (!is.null(e$sel)) e$sel else list()
  .ckpt$partial <- if (!is.null(e$partial)) e$partial else list()
}
ckpt_save <- function() {
  sel <- .ckpt$sel; partial <- .ckpt$partial
  save(sel, partial, cache_signature, file = cfg$cache_file[1])
}

select_copula <- function(component, scenario, dedtype, init_marg, tau0, prev, sev_spec = NULL, families = cfg$sel_families, key = NULL, nu0 = 5) {
  fams <- if (is.null(families)) cop_families$fam else families
  if (cfg$include_independence) fams <- unique(c(0L, fams))
  done <- if (!is.null(key)) .ckpt$partial[[key]] else NULL
  fits <- lapply(fams, function(fam) {
    if (!is.null(done[[fam_name(fam)]])) { cat("   (checkpoint)", component, scenario, dedtype, fam_name(fam), "\n"); return(done[[fam_name(fam)]]) }
    cb <- cop_bounds(fam)
    init <- if (fam == 0L) init_marg else
      c(init_marg,
        min(max(tau0, cb$lo[1] + 1e-3), cb$hi[1] - 1e-3),
        if (fam == 2L) nu0)
    st <- if (component == "sev") stage_sev(paste(scenario, dedtype, fam), scenario, dedtype, fam)
          else stage_freq(paste(scenario, dedtype, fam), scenario, dedtype, fam, sev_spec)
    t0 <- proc.time()[3]
    fit <- tryCatch(fit_stage(st, init, prev, naive_hessian = FALSE), error = function(e) { message("   ", fam_name(fam), ": ", conditionMessage(e)); NULL })
    cat(sprintf("   %-4s %-10s %-5s %-16s %6.0fs  loglik %s\n", component, scenario, dedtype, fam_name(fam), proc.time()[3] - t0,
                if (is.null(fit)) "NA" else sprintf("%.2f", fit$loglik))); flush.console()
    if (!is.null(key) && !is.null(fit)) { p <- .ckpt$partial[[key]]; if (is.null(p)) p <- list(); p[[fam_name(fam)]] <- fit; .ckpt$partial[[key]] <- p; ckpt_save() }
    fit
  })
  names(fits) <- fam_name(fams)
  tab <- t(sapply(fits, function(f) if (is.null(f)) c(Likelihood = NA, AIC = NA, BIC = NA, tau = NA)
                                    else c(Likelihood = f$loglik, AIC = f$AIC, BIC = f$BIC, tau = unname(f$par["tau"]))))
  list(table = tab, fits = fits)
}
best_fit <- function(selres, name) {
  # Compare AIC only among numerical fits that have converged.
  okay <- vapply(selres$fits, function(f)
    !is.null(f) && is.finite(f$AIC) &&
      !is.null(f$convergence) && f$convergence == 0L &&
      all(is.finite(f$par)), logical(1))
  ii <- which(okay)
  if (!length(ii)) {
    if (isTRUE(cfg$selection_require_convergence))
      stop("No converged copula candidate for ", name,
           "; inspect cached selection diagnostics or refit candidates.")
    ii <- which(is.finite(selres$table[, "AIC"]))
    if (!length(ii)) stop("No finite copula candidate for ", name)
    warning("Using a nonconverged finite-AIC fit for ", name)
  }
  i <- ii[which.min(selres$table[ii, "AIC"])]
  f <- selres$fits[[i]]
  f$name <- name
  f
}

################################################################################
## 5. Severity copula selection FIRST: six SCM/SDM models (AIC best)
################################################################################
cat("\n==== 2. Severity copula selection before FEM ====\n")

# Preserve the original nine named families, and optionally add independence.
fam_name_old <- fam_name
fam_name <- function(fam) {
  nm <- fam_name_old(fam)
  nm[fam == 0L] <- "Independence"
  nm
}

# Bind the checkpoint's provenance to the selected training period.
# The cache file should be cleared when changing model definitions, data, or settings.
fit_selection_key <- function(key, ...) {
  if (is.null(.ckpt$sel[[key]])) {
    .ckpt$sel[[key]] <- select_copula(..., key = key)
    ckpt_save()
  }
  .ckpt$sel[[key]]
}
if (cfg$run_selection) ckpt_load() else {
  if (!file.exists(cfg$cache_file)) stop("No selection cache available: ", cfg$cache_file)
  ckpt_load()
}

if (cfg$run_selection) {
  for (dt in cfg$sel_dedtypes) {
    fitD <- if (dt == "cont") fit_Dc else fit_Dd
    for (sc in c("complete", "censored", "truncated")) {
      marg <- switch(sc, complete = SEM1, censored = SEM2, truncated = SEM3)
      dat <- if (sc == "truncated") dat_pos else dat_clm
      tau0 <- suppressWarnings(cor(F_D_of(dt, dat, fitD$par),
        marg_cdf_sev(marg, dat, sc), method = "kendall"))
      if (!is.finite(tau0)) tau0 <- 0.1
      init_marginal <- marg$par
      init_nu <- 5
      if (identical(sc, "censored") && isTRUE(cfg$censored_init_from_complete)) {
        # SCM2/SDM2 start from SCM1/SDM1 (complete-data fitted values).
        # Use the selected complete family; candidate copula dimensions may
        # differ, so transfer *marginal* parameters and Kendall's tau,
        # plus t degrees of freedom only when the candidate is Student-t.
        complete_key <- paste("sev", dt, "complete", sep = "_")
        base_name <- if (dt == "cont") "SCM1" else "SDM1"
        if (isTRUE(cfg$auto_select)) {
          base_fit <- best_fit(.ckpt$sel[[complete_key]], base_name)
        } else {
          fname <- fam_name(cfg$paper_families[[base_name]])
          base_fit <- .ckpt$sel[[complete_key]]$fits[[fname]]
          if (is.null(base_fit)) stop("Missing complete-model starting fit: ", base_name)
        }
        init_marginal <- unname(base_fit$par[seq_along(marg$par)])
        if ("tau" %in% names(base_fit$par))
          tau0 <- unname(base_fit$par["tau"])
        if ("nu" %in% names(base_fit$par))
          init_nu <- unname(base_fit$par["nu"])
        if (any(!is.finite(init_marginal)) || !is.finite(tau0))
          stop("Invalid initial values inherited from ", base_name)
        cat("  ", if (dt == "cont") "SCM2" else "SDM2",
            " initialized from ", base_name, " (complete fit)\n", sep = "")
      }
      fit_selection_key(paste("sev", dt, sc, sep = "_"),
        "sev", sc, dt, init_marginal, tau0,
        prev = list(D = fitD$par), nu0 = init_nu)
    }
  }
}
sel <- .ckpt$sel
required_sev <- as.vector(outer(c("cont", "disc"),
                              c("complete", "censored", "truncated"),
                              function(dt, sc) paste("sev", dt, sc, sep = "_")))
if (!all(required_sev %in% names(sel)))
  stop("Severity copula selection is incomplete. Missing: ",
       paste(setdiff(required_sev, names(sel)), collapse = ", "))

add_naive_hessian <- function(f) {
  if (!cfg$sandwich) return(f)
  nll <- function(th) -sum(f$contrib(th, f$prev))
  f$hess <- tryCatch(numDeriv::hessian(nll, f$par,
    method.args = list(d = cfg$hess_d, r = cfg$hess_r)), error = function(e) NULL)
  f$se_naive <- tryCatch(sqrt(pmax(diag(solve(f$hess)), 0)),
    error = function(e) rep(NA_real_, length(f$par)))
  f
}
final_model <- function(name, key) {
  if (cfg$auto_select) f <- best_fit(sel[[key]], name)
  else {
    wanted <- fam_name(cfg$paper_families[[name]])
    f <- sel[[key]]$fits[[wanted]]
    if (is.null(f)) stop("Missing requested family ", wanted, " for ", name)
    f$name <- name
  }
  cat(sprintf("   %-5s <- %-17s  AIC %.2f  tau %s\n", name,
    fam_name(f$spec$fam), f$AIC,
    if ("tau" %in% names(f$par)) sprintf("%.3f", f$par["tau"]) else "NA"))
  add_naive_hessian(f)
}
SCM1 <- final_model("SCM1", "sev_cont_complete")
SCM2 <- final_model("SCM2", "sev_cont_censored")
SCM3 <- final_model("SCM3", "sev_cont_truncated")
SDM1 <- final_model("SDM1", "sev_disc_complete")
SDM2 <- final_model("SDM2", "sev_disc_censored")
SDM3 <- final_model("SDM3", "sev_disc_truncated")

# All candidates are reported; no family is hard-coded to the manuscript's old winner.
selection_table <- function(keys, type) {
  do.call(rbind, lapply(keys, function(key) {
    z <- sel[[key]]$table
    fitted <- sel[[key]]$fits
    par1 <- vapply(fitted, function(f) {
      if (is.null(f)) return(NA_real_)
      if (f$spec$fam == 0L) return(0)
      unname(tau2par(f$spec$fam, f$par["tau"]))
    }, numeric(1))
    nu <- vapply(fitted, function(f) {
      if (is.null(f) || !"nu" %in% names(f$par)) return(NA_real_)
      unname(f$par["nu"])
    }, numeric(1))
    data.frame(Component = type, Specification = key, Family = rownames(z),
      logLik = z[, "Likelihood"], AIC = z[, "AIC"], BIC = z[, "BIC"],
      Kendall_tau = z[, "tau"], Copula_par = unname(par1),
      Degrees_of_freedom = unname(nu), row.names = NULL)
  }))
}
severity_candidates <- selection_table(required_sev, "severity")
write.csv(severity_candidates, output_file("Severity_Copula_AIC_AllFamilies.csv"), row.names = FALSE)
selected_severity <- do.call(rbind, lapply(list(
  SCM1 = SCM1, SCM2 = SCM2, SCM3 = SCM3,
  SDM1 = SDM1, SDM2 = SDM2, SDM3 = SDM3), function(f) {
    data.frame(Model = f$name, Deductible = f$spec$dedtype,
      Scenario = f$spec$scenario, Family = fam_name(f$spec$fam),
      Family_code = f$spec$fam, Kendall_tau = if ("tau" %in% names(f$par))
        unname(f$par["tau"]) else 0,
      Copula_par = if (f$spec$fam == 0L) 0 else
        unname(tau2par(f$spec$fam, f$par["tau"])),
      nu = if ("nu" %in% names(f$par)) unname(f$par["nu"]) else NA_real_,
      AIC = f$AIC, BIC = f$BIC, convergence = f$convergence)
  }))
write.csv(selected_severity, output_file("Severity_Copula_Selected.csv"), row.names = FALSE)
################################################################################
## 6. Frequency models AFTER severity copula estimation
################################################################################
cat("\n==== 3. FEM thinning source comparisons ====\n")
init_freq <- c(coef(glm(pol$Freq ~ X_pol - 1, family = poisson)), 1)
sev_none <- list(dedtype = "none", fam = 0L)
FEM1 <- fit_stage(stage_freq("FEM1", "censored", "none"), init_freq)

# Original SEM3 thinning under marginal exogeneity (XR3)
FEM2_SEM3 <- fit_stage(stage_freq("FEM2_SEM3", "truncated", "none",
  sev_spec = sev_none), FEM1$par, prev = list(S = SEM3$par))
# Additional benchmark: SEM1 thinning (requested FEM2G_SEM1)
FEM2G_SEM1 <- fit_stage(stage_freq("FEM2G_SEM1", "truncated", "none",
  sev_spec = sev_none), FEM1$par, prev = list(S = SEM1$par))
# Conditional, automatically selected SCM3 and SDM3 thinning
FEM2G_cont <- fit_stage(stage_freq("FEM2G_cont", "truncated", "none",
  sev_spec = list(dedtype = "cont", fam = SCM3$spec$fam)),
  FEM1$par, prev = list(D = fit_Dc$par, S = SCM3$par))
FEM2G_disc <- fit_stage(stage_freq("FEM2G_disc", "truncated", "none",
  sev_spec = list(dedtype = "disc", fam = SDM3$spec$fam)),
  FEM1$par, prev = list(D = fit_Dd$par, S = SDM3$par))

fem_table <- do.call(rbind, lapply(list(FEM1, FEM2_SEM3, FEM2G_SEM1,
                                      FEM2G_cont, FEM2G_disc), function(f) {
  data.frame(Model = f$name, Parameter = names(f$par), Estimate = unname(f$par),
             Stagewise_SE = if (!is.null(f$se_naive)) f$se_naive else NA_real_,
             logLik = f$loglik, AIC = f$AIC, convergence = f$convergence)
}))
write.csv(fem_table, output_file("FEM_Thinning_Comparisons.csv"), row.names = FALSE)

cat("\n==== 4. Frequency copula selection (Case 4 alternative) ====\n")
if (cfg$run_selection) {
  for (dt in cfg$sel_dedtypes) {
    fitD <- if (dt == "cont") fit_Dc else fit_Dd
    S3 <- if (dt == "cont") SCM3 else SDM3
    f3 <- if (dt == "cont") FEM2G_cont else FEM2G_disc
    ss <- list(dedtype = dt, fam = S3$spec$fam)
    tau0 <- suppressWarnings(cor(F_D_of(dt, dat_pol, fitD$par),
      pnbinom(dat_pol$n_cen, size = FEM1$par[length(FEM1$par)],
        mu = drop(exp(X_pol %*% FEM1$par[seq_len(ncol(X_pol))]))),
      method = "kendall"))
    if (!is.finite(tau0)) tau0 <- -0.1
    fit_selection_key(paste("freq", dt, "censored", sep = "_"),
      "freq", "censored", dt, FEM1$par, tau0,
      prev = list(D = fitD$par))
    # Initial tau0 below is based on the observed deductible/count rank association.
    # The truncation likelihood itself uses p(d;x) from the chosen S3 model.
    tau0 <- suppressWarnings(cor(F_D_of(dt, dat_pol, fitD$par), dat_pol$n_tru,
      method = "kendall"))
    if (!is.finite(tau0)) tau0 <- -0.1
    fit_selection_key(paste("freq", dt, "truncated", sep = "_"),
      "freq", "truncated", dt, f3$par, tau0,
      prev = list(D = fitD$par, S = S3$par), sev_spec = ss)
  }
}
sel <- .ckpt$sel
required_freq <- as.vector(outer(c("cont", "disc"), c("censored", "truncated"),
  function(dt, sc) paste("freq", dt, sc, sep = "_")))
if (!all(required_freq %in% names(sel))) stop("Frequency selection incomplete: ",
  paste(setdiff(required_freq, names(sel)), collapse = ", "))
FCM1 <- final_model("FCM1", "freq_cont_censored")
FCM2 <- final_model("FCM2", "freq_cont_truncated")
FDM1 <- final_model("FDM1", "freq_disc_censored")
FDM2 <- final_model("FDM2", "freq_disc_truncated")

frequency_candidates <- selection_table(required_freq, "frequency")
write.csv(frequency_candidates, output_file("Frequency_Copula_AIC_AllFamilies.csv"), row.names = FALSE)
selected_frequency <- do.call(rbind, lapply(list(FCM1 = FCM1, FCM2 = FCM2,
                                                FDM1 = FDM1, FDM2 = FDM2), function(f) {
  data.frame(Model = f$name, Deductible = f$spec$dedtype,
    Scenario = f$spec$scenario, Family = fam_name(f$spec$fam),
    Family_code = f$spec$fam, Kendall_tau = if ("tau" %in% names(f$par))
      unname(f$par["tau"]) else 0,
    nu = if ("nu" %in% names(f$par)) unname(f$par["nu"]) else NA_real_,
    AIC = f$AIC, BIC = f$BIC, convergence = f$convergence)
}))
write.csv(selected_frequency, output_file("Frequency_Copula_Selected.csv"), row.names = FALSE)

# Full fitted parameter vectors (regression and dependence parameters) for all
# selected SCM/SDM/FCM/FDM stages; tau is on Kendall's-tau scale.
selected_fit_list <- list(SCM1 = SCM1, SCM2 = SCM2, SCM3 = SCM3,
  SDM1 = SDM1, SDM2 = SDM2, SDM3 = SDM3,
  FCM1 = FCM1, FCM2 = FCM2, FDM1 = FDM1, FDM2 = FDM2)
selected_parameters <- do.call(rbind, lapply(selected_fit_list, function(f)
  data.frame(Model = f$name, Family = fam_name(f$spec$fam),
             Parameter = names(f$par), Estimate = unname(f$par),
             AIC = f$AIC, Convergence = f$convergence)))
write.csv(selected_parameters,
  output_file("Selected_Copula_AllParameters.csv"), row.names = FALSE)

# Preserve original moment/sandwich analysis, but do not run expensive numeric
# derivatives by default. For FCM2/FDM2, the chain is D -> SCM3/SDM3 -> N.
sw <- list()
if (cfg$sandwich) {
  cat("\n==== Full Godambe standard errors ====\n")
  chains <- list(
    SEM1 = list(SEM1), SEM2 = list(SEM2), SEM3 = list(SEM3),
    SCM1 = list(fit_Dc, SCM1), SCM2 = list(fit_Dc, SCM2), SCM3 = list(fit_Dc, SCM3),
    SDM1 = list(fit_Dd, SDM1), SDM2 = list(fit_Dd, SDM2), SDM3 = list(fit_Dd, SDM3),
    FEM1 = list(FEM1), FEM2_SEM3 = list(SEM3, FEM2_SEM3),
    FCM1 = list(fit_Dc, FCM1), FCM2 = list(fit_Dc, SCM3, FCM2),
    FDM1 = list(fit_Dd, FDM1), FDM2 = list(fit_Dd, SDM3, FDM2)
  )
  for (nm in names(chains)) {
    cat("   Godambe chain:", nm, "\n")
    sw[[nm]] <- tryCatch(sandwich_se(chains[[nm]], verbose = FALSE),
      error = function(e) { warning("Godambe ", nm, ": ", conditionMessage(e)); NULL })
    saveRDS(sw, output_file("Godambe_Checkpoint.rds"))
  }
  se_tables <- do.call(rbind, lapply(names(sw), function(nm) {
    if (is.null(sw[[nm]])) return(NULL)
    cbind(Chain = nm, last_stage_table(sw[[nm]]))
  }))
  if (!is.null(se_tables))
    write.csv(se_tables, output_file("Godambe_LastStage_Summary.csv"), row.names = FALSE)
}
if (cfg$run_severity_family_table) {
  sev_family_tab <- sev_family_table()
  saveRDS(sev_family_tab, output_file("Severity_Family_Selection.rds"))
}
################################################################################
## 7. In-sample aggregate-payment evaluation (2010 OOS at bottom)
##
## IMPORTANT DISTINCTIONS
##  * The LGPIF data have observed payments, NOT a known true conditional premium.
##    Thus we report MAE_obs / RMSRE_obs / calibration on held-out policy-years
##    where possible; these are not simulation MAE_true / RMSRE_true.
##  * RMSRE_obs uses only observations with S_obs > 0 (because dividing by
##    zero is undefined). N_RMSRE gives the corresponding sample size.
##    MAE_obs and calibration continue to use ALL finite observations.
##  * A first GB2 moment exists iff alpha2 > sigma. When it does not,
##    we refuse to treat a finite numerical integral as an unlimited premium.
##  * PRIMARY: ground-up per-occurrence coverage C_i = Coverage_i *
##    cfg$coverage_to_dollars. Observed and predicted limited payments BOTH use
##    (min(Y_ij, C_i) - D_i)_+; expected claim payment integrates survival over
##    [D_i, C_i], NOT [D_i, D_i + C_i]. This contract limit is assumed.
##  * SECONDARY: unlimited (Y_ij-D_i)_+; premium is reported only when the fitted
##    severity first moment and quadrature permit a reliable calculation.
################################################################################
cat("\n==== 5. Premium comparison, Case 1 vs Cases 4c/4d ====\n")

# ================================================================
# Evaluation data
# ================================================================

# Coverage is a per-occurrence GROUND-UP loss limit.
#
# C_i = Coverage_i * cfg$coverage_to_dollars
#
# Unlimited payment:
#   (Y_ij - D_i)_+
#
# Coverage-limited payment:
#   (min(Y_ij, C_i) - D_i)_+
#
# If C_i <= D_i, limited payments are exactly zero.

prepare_eval <- function(epol, eclm) {
  
  epol <- epol[
    order(epol$Year, epol$PolicyNum),
    , drop = FALSE
  ]
  
  epol$discDed <- match(
    epol$Deduct,
    cfg$ded_levels
  )
  
  if (anyNA(epol$discDed))
    stop("Evaluation has an unrecognized deductible level")
  
  epol$dedratio <- epol$Deduct / epol$Coverage
  
  if (any(!is.finite(epol$dedratio) |
          epol$dedratio <= 0))
    stop("Evaluation has invalid deductible/coverage ratio")
  
  # Do not rescale Coverage in the regression design matrix.
  # Convert only for monetary payment calculations.
  epol$CoverageLimit_Dollars <-
    epol$Coverage * cfg$coverage_to_dollars
  
  if (any(!is.finite(epol$CoverageLimit_Dollars) |
          epol$CoverageLimit_Dollars <= 0))
    stop("Invalid policy-specific coverage limits")
  
  # IMPORTANT CHANGE:
  # Keep policies with C_i <= D_i.
  # They have no positive coverage-limited payment.
  bad_limits <- which(
    epol$CoverageLimit_Dollars <= epol$Deduct
  )
  
  epol$CoverageBelowOrEqualDeductible <-
    epol$CoverageLimit_Dollars <= epol$Deduct
  
  if (length(bad_limits)) {
    warning(
      length(bad_limits),
      " policy-years have C_i <= D_i; ",
      "their coverage-limited payments and premiums ",
      "will be zero. Rows: ",
      paste(head(bad_limits, 10L), collapse = ", ")
    )
  }
  
  # Match claims with policy-years
  keypol <- paste(
    epol$PolicyNum,
    epol$Year,
    sep = "|"
  )
  
  keyclm <- paste(
    eclm$PolicyNum,
    eclm$Year,
    sep = "|"
  )
  
  if (anyDuplicated(keypol))
    stop("PolicyNum+Year not unique in evaluation data")
  
  idx <- match(keyclm, keypol)
  
  if (anyNA(idx)) {
    warning(
      sum(is.na(idx)),
      " evaluation claims lack matching policy-year records"
    )
  }
  
  keep <- which(!is.na(idx))
  idx <- idx[keep]
  
  if (any(eclm$Deduct[keep] != epol$Deduct[idx]))
    stop("Inconsistent deductibles in matched claims")
  
  groundup_y <- eclm$LossBeforeDeductible[keep]
  
  # Unlimited insured payment
  payments <- pmax(
    groundup_y - epol$Deduct[idx],
    0
  )
  
  # Policy-specific coverage-limited insured payment
  payments_limited <- pmax(
    pmin(
      groundup_y,
      epol$CoverageLimit_Dollars[idx]
    ) - epol$Deduct[idx],
    0
  )
  
  if (any(!is.finite(payments)) ||
      any(!is.finite(payments_limited)))
    stop("Nonfinite observed claim payments")
  
  # Aggregate claim-level amounts by policy-year
  by_unit <- function(z) {
    
    out <- numeric(nrow(epol))
    
    if (length(z)) {
      zz <- rowsum(
        z,
        idx,
        reorder = TRUE
      )
      
      out[as.integer(rownames(zz))] <-
        as.numeric(zz[, 1L])
    }
    
    out
  }
  
  epol$S_raw <- by_unit(payments)
  
  epol$S_limited <- by_unit(payments_limited)
  
  epol$n_ground <- by_unit(
    rep(1, length(payments))
  )
  
  epol$n_paid_unlimited <- by_unit(
    as.numeric(payments > 0)
  )
  
  # IMPORTANT CHANGE:
  # Limited payment count is based on the LIMITED payment,
  # not simply on Y > D.
  epol$n_paid <- by_unit(
    as.numeric(payments_limited > 0)
  )
  
  epol
}

eval_pol <- prepare_eval(
  pol_evaluate_raw,
  clm_evaluate_raw
)

# Verify the C_i <= D_i condition
zero_limit <- eval_pol$CoverageBelowOrEqualDeductible

stopifnot(
  all(eval_pol$S_limited[zero_limit] == 0),
  all(eval_pol$n_paid[zero_limit] == 0)
)

if (any(
  eval_pol$S_limited < 0 |
  eval_pol$S_limited > eval_pol$S_raw + 1e-7
)) {
  stop(
    "Limited aggregate losses must lie between ",
    "zero and unlimited aggregate losses"
  )
}

# Coverage limit audit
write.csv(
  data.frame(
    Coverage_original_min = min(eval_pol$Coverage),
    Coverage_original_median = median(eval_pol$Coverage),
    Coverage_original_max = max(eval_pol$Coverage),
    
    Conversion_to_dollars = cfg$coverage_to_dollars,
    
    Coverage_limit_dollars_min =
      min(eval_pol$CoverageLimit_Dollars),
    
    Coverage_limit_dollars_median =
      median(eval_pol$CoverageLimit_Dollars),
    
    Coverage_limit_dollars_max =
      max(eval_pol$CoverageLimit_Dollars),
    
    N_policies = nrow(eval_pol),
    
    N_coverage_below_or_equal_deductible =
      sum(eval_pol$CoverageBelowOrEqualDeductible),
    
    Limit_interpretation =
      "Per-occurrence ground-up loss limit (assumed)"
  ),
  
  output_file("Coverage_Limit_Audit.csv"),
  row.names = FALSE
)

# Evaluation design matrix
X_eval <- model.matrix(
  fml,
  data = eval_pol
)

stopifnot(
  identical(colnames(X_eval), colnames(X_pol))
)

eval_dat <- list(
  X = X_eval,
  d = eval_pol$Deduct,
  coverage_limit = eval_pol$CoverageLimit_Dollars,
  dratio = eval_pol$dedratio,
  dlev = eval_pol$discDed,
  n_cen = eval_pol$n_ground,
  n_tru = eval_pol$n_paid,
  unit = seq_len(nrow(eval_pol))
)

# Quadrature nodes
GL_pay <- gl_nodes(cfg$premium_gl_n)
GL_unlim <- gl_nodes(cfg$unlimited_gl_n)
GL_unlim_hi <- gl_nodes(
  2L * cfg$unlimited_gl_n
)

# ================================================================
# Conditional severity CDF
# SEM, SCM, SDM
# ================================================================

conditional_Y_cdf <- function(fit, theta_D, dat, y) {
  
  th <- fit$par
  sp <- fit$spec
  k <- ncol(dat$X)
  
  mu <- drop(
    dat$X %*% th[seq_len(k)]
  )
  
  Fy <- pGB2(
    y,
    mu,
    th[k + 1L],
    th[k + 2L],
    th[k + 3L]
  )
  
  # Marginal severity
  if (sp$dedtype == "none")
    return(pmin(1, pmax(0, Fy)))
  
  fam <- sp$fam
  
  par <- if (fam == 0L) {
    0
  } else {
    tau2par(fam, th[k + 4L])
  }
  
  par2 <- if (fam == 2L) th[k + 5L] else 0
  
  # Continuous deductible
  if (sp$dedtype == "cont") {
    
    u <- D_cont_cdf(
      theta_D,
      XD(dat$X, "cont"),
      dat$dratio
    )
    
    ans <- cop(
      "h1",
      u,
      Fy,
      fam,
      par,
      par2
    )
    
  } else {
    
    # Discrete deductible: integrate over the D atom
    up <- D_disc_cdf(
      theta_D,
      XD(dat$X, "disc"),
      dat$dlev
    )
    
    um <- D_disc_cdf(
      theta_D,
      XD(dat$X, "disc"),
      dat$dlev - 1L
    )
    
    w <- up - um
    
    if (any(w <= 0))
      stop("Nonpositive deductible atom")
    
    ans <- cop_rect(
      um,
      up,
      Fy,
      fam,
      par,
      par2
    ) / w
  }
  
  pmin(1, pmax(0, ans))
}

# ================================================================
# Conditional ground-up frequency CDF
# FEM, FCM, FDM
# ================================================================

conditional_N_cdf <- function(fit, theta_D, dat, n) {
  
  th <- fit$par
  sp <- fit$spec
  k <- ncol(dat$X)
  
  mu <- drop(
    exp(dat$X %*% th[seq_len(k)])
  )
  
  r <- th[k + 1L]
  
  Fn <- pnbinom(
    n,
    size = r,
    mu = mu
  )
  
  if (sp$dedtype == "none")
    return(Fn)
  
  fam <- sp$fam
  
  par <- if (fam == 0L) {
    0
  } else {
    tau2par(fam, th[k + 2L])
  }
  
  par2 <- if (fam == 2L) th[k + 3L] else 0
  
  if (sp$dedtype == "cont") {
    
    u <- D_cont_cdf(
      theta_D,
      XD(dat$X, "cont"),
      dat$dratio
    )
    
    ans <- cop(
      "h1",
      u,
      Fn,
      fam,
      par,
      par2
    )
    
  } else {
    
    up <- D_disc_cdf(
      theta_D,
      XD(dat$X, "disc"),
      dat$dlev
    )
    
    um <- D_disc_cdf(
      theta_D,
      XD(dat$X, "disc"),
      dat$dlev - 1L
    )
    
    if (any(up <= um))
      stop("Nonpositive deductible atom")
    
    ans <- cop_rect(
      um,
      up,
      Fn,
      fam,
      par,
      par2
    ) / (up - um)
  }
  
  pmin(1, pmax(0, ans))
}

# ================================================================
# Conditional frequency mean E[N | D, X]
# ================================================================

conditional_N_mean <- function(fit, theta_D, dat) {
  
  th <- fit$par
  sp <- fit$spec
  k <- ncol(dat$X)
  
  mu <- drop(
    exp(dat$X %*% th[seq_len(k)])
  )
  
  # Marginal frequency model
  if (sp$dedtype == "none" || sp$fam == 0L)
    return(mu)
  
  if (any(!is.finite(mu)))
    stop("Nonfinite negative-binomial mean")
  
  if (sp$dedtype == "disc") {
    
    up <- D_disc_cdf(
      theta_D,
      XD(dat$X, "disc"),
      dat$dlev
    )
    
    um <- D_disc_cdf(
      theta_D,
      XD(dat$X, "disc"),
      dat$dlev - 1L
    )
    
    minw <- min(up - um)
    
  } else {
    minw <- 1
  }
  
  eps <- pmax(
    1e-12,
    min(1e-6, cfg$count_tail_eps * minw)
  )
  
  pp <- 1 - eps
  
  if (pp >= 1)
    pp <- 1 - 1e-12
  
  K <- qnbinom(
    pp,
    size = th[k + 1L],
    mu = max(mu)
  )
  
  if (!is.finite(K) ||
      K > cfg$count_max_k) {
    
    warning(
      "Frequency cutoff exceeds count_max_k for ",
      fit$name
    )
    
    return(rep(
      NA_real_,
      length(mu)
    ))
  }
  
  res <- numeric(length(mu))
  
  for (j in 0:as.integer(K)) {
    
    H <- conditional_N_cdf(
      fit,
      theta_D,
      dat,
      j
    )
    
    if (any(!is.finite(H)))
      return(rep(
        NA_real_,
        length(mu)
      ))
    
    res <- res + pmax(0, 1 - H)
  }
  
  res
}

# ================================================================
# Coverage-limited expected payment severity
#
# m_i^C =
#   integral_{D_i}^{C_i}
#   [1-F_{Y|D}(y|D_i,X_i)] dy
#
# m_i^C = 0 when C_i <= D_i
# ================================================================

conditional_paid_limited <- function(
    fit, theta_D, dat, gl = GL_pay) {
  
  width <- dat$coverage_limit - dat$d
  
  if (length(width) != nrow(dat$X) ||
      any(!is.finite(width))) {
    stop(
      "Invalid policy-specific deductible or coverage limit"
    )
  }
  
  # Start with zero for every policy.
  ans <- numeric(nrow(dat$X))
  
  # Integrate only when C_i > D_i.
  active <- which(width > 0)
  
  if (!length(active))
    return(ans)
  
  # Subset the entire evaluation data consistently.
  active_dat <- lapply(dat, function(v) {
    
    if (is.matrix(v)) {
      v[active, , drop = FALSE]
    } else {
      v[active]
    }
  })
  
  active_width <- width[active]
  z <- numeric(length(active))
  
  # Change of variables:
  # y = D_i + (C_i-D_i)*t^3
  for (j in seq_along(gl$x)) {
    
    tj <- (1 + gl$x[j]) / 2
    
    y <- active_dat$d +
      active_width * tj^3
    
    sf <- 1 - conditional_Y_cdf(
      fit,
      theta_D,
      active_dat,
      y
    )
    
    if (any(!is.finite(sf))) {
      stop(
        "Nonfinite conditional survival probability ",
        "in coverage-limited premium"
      )
    }
    
    z <- z + (gl$w[j] / 2) *
      sf * 3 * active_width * tj^2
  }
  
  ans[active] <- pmin(
    active_width,
    pmax(z, 0)
  )
  
  ans
}


# ================================================================
# GB2 first-moment diagnostics
# ================================================================

severity_moment <- function(fit) {
  k <- ncol(X_eval)
  sigma <- unname(fit$par[k + 1L])
  a2 <- unname(fit$par[k + 3L])
  
  data.frame(
    Model = fit$name,
    sigma = sigma,
    alpha2 = a2,
    Moment_index = a2 / sigma,
    Marginal_mean_finite =
      is.finite(a2 / sigma) && a2 > sigma
  )
}

# Marginal GB2 stop-loss expectation, when finite
GB2_stoploss <- function(mu, d, sigma, a1, a2) {
  
  if (a2 <= sigma)
    return(rep(NA_real_, length(d)))
  
  v <- plogis((log(d) - mu) / sigma)
  
  EYtail <- exp(
    mu + lbeta(a1 + sigma, a2 - sigma) -
      lbeta(a1, a2)
  ) * pbeta(
    v, a1 + sigma, a2 - sigma,
    lower.tail = FALSE
  )
  
  p_tail <- pbeta(
    v, a1, a2,
    lower.tail = FALSE
  )
  
  pmax(0, EYtail - d * p_tail)
}

# ================================================================
# Unlimited expected payment severity
# Supplementary theoretical analysis
# ================================================================

conditional_paid_unlimited_once <- function(
    fit, theta_D, dat, gl) {
  
  k <- ncol(dat$X)
  th <- fit$par
  
  if (th[k + 3L] <= th[k + 1L])
    return(rep(NA_real_, nrow(dat$X)))
  
  mu <- drop(
    dat$X %*% th[seq_len(k)]
  )
  
  if (fit$spec$dedtype == "none" ||
      fit$spec$fam == 0L) {
    
    return(GB2_stoploss(
      mu,
      dat$d,
      th[k + 1L],
      th[k + 2L],
      th[k + 3L]
    ))
  }
  
  scale <- exp(mu)
  ans <- numeric(nrow(dat$X))
  
  # Unbounded transformation:
  # y = d + scale*t/(1-t)
  for (j in seq_along(gl$x)) {
    
    tj <- (1 + gl$x[j]) / 2
    
    y <- dat$d + scale * tj / (1 - tj)
    
    jac <- scale / (1 - tj)^2
    
    sf <- pmax(
      0,
      1 - conditional_Y_cdf(
        fit, theta_D, dat, y
      )
    )
    
    ans <- ans +
      (gl$w[j] / 2) * sf * jac
  }
  
  ans
}

conditional_paid_unlimited <- function(
    fit, theta_D, dat) {
  
  m <- severity_moment(fit)
  
  if (!m$Marginal_mean_finite)
    return(rep(NA_real_, nrow(dat$X)))
  
  a <- conditional_paid_unlimited_once(
    fit, theta_D, dat, GL_unlim
  )
  
  if (fit$spec$dedtype == "none" ||
      fit$spec$fam == 0L)
    return(a)
  
  b <- conditional_paid_unlimited_once(
    fit, theta_D, dat, GL_unlim_hi
  )
  
  bad <- !is.finite(a) |
    !is.finite(b) |
    abs(a - b) >
    cfg$unlimited_rel_tol * pmax(1, abs(b))
  
  if (any(bad)) {
    warning(
      "Unbounded stop-loss quadrature unstable for ",
      fit$name, " in ", sum(bad),
      " policies; predictions excluded"
    )
  }
  
  b[bad] <- NA_real_
  b
}

# ================================================================
# Collect fitted models
# ================================================================

frequency_fits <- list(
  FEM1 = FEM1,
  FEM2_SEM3 = FEM2_SEM3,
  FEM2G_SEM1 = FEM2G_SEM1,
  FEM2G_cont = FEM2G_cont,
  FEM2G_disc = FEM2G_disc,
  FCM1 = FCM1,
  FCM2 = FCM2,
  FDM1 = FDM1,
  FDM2 = FDM2
)

severity_fits <- list(
  SEM1 = SEM1,
  SEM2 = SEM2,
  SEM3 = SEM3,
  SCM1 = SCM1,
  SCM2 = SCM2,
  SCM3 = SCM3,
  SDM1 = SDM1,
  SDM2 = SDM2,
  SDM3 = SDM3
)

moments <- do.call(
  rbind,
  lapply(severity_fits, severity_moment)
)

write.csv(
  moments,
  output_file("GB2_Moment_Existence.csv"),
  row.names = FALSE
)

cat(
  "\nGB2 first-moment diagnostics ",
  "(FALSE = unlimited mean unavailable):\n"
)

print(moments, row.names = FALSE)

Dtheta_for <- function(fit) {
  switch(
    fit$spec$dedtype,
    cont = fit_Dc$par,
    disc = fit_Dd$par,
    none = NULL
  )
}

# ================================================================
# Predict conditional frequency and severity components
# ================================================================

cat("\nPredicting frequency conditional means...\n")

# Only six frequency components are used in the nine requested premium pairs.
# Extra FEM2G fits remain available for diagnostic/selection calculations.
premium_frequency_models <- c("FEM1", "FEM2_SEM3", "FCM1", "FCM2", "FDM1", "FDM2")
freq_mean_pred <- lapply(
  frequency_fits[premium_frequency_models],
  function(f) {
    conditional_N_mean(
      f, Dtheta_for(f), eval_dat
    )
  }
)

cat(
  "Predicting policy-coverage-limited ",
  "severity payments...\n"
)

sev_limited_pred <- lapply(
  severity_fits,
  function(f) {
    conditional_paid_limited(
      f, Dtheta_for(f), eval_dat
    )
  }
)

# IMPORTANT:
# For C_i <= D_i, all severity expectations must be zero.
stopifnot(
  all(vapply(
    sev_limited_pred,
    function(z) all(z[zero_limit] == 0),
    logical(1)
  ))
)

# ================================================================
# Limited-premium quadrature convergence checks
# ================================================================

quad_checks <- NULL

if (cfg$capped_quad_check) {
  
  idx_by_d <- split(
    seq_len(nrow(eval_pol)),
    eval_pol$Deduct
  )
  
  ix <- sort(unique(unlist(
    lapply(idx_by_d, function(i) {
      
      j <- unique(round(seq(
        1,
        length(i),
        length.out = min(5L, length(i))
      )))
      
      i[j]
    }),
    use.names = FALSE
  )))
  
  small_dat <- lapply(eval_dat, function(v) {
    
    if (is.matrix(v))
      return(v[ix, , drop = FALSE])
    
    v[ix]
  })
  
  gl_high <- gl_nodes(
    2L * cfg$premium_gl_n
  )
  
  quad_checks <- do.call(
    rbind,
    lapply(severity_fits, function(f) {
      
      low <- conditional_paid_limited(
        f, Dtheta_for(f),
        small_dat, GL_pay
      )
      
      high <- conditional_paid_limited(
        f, Dtheta_for(f),
        small_dat, gl_high
      )
      
      rel <- abs(low - high) /
        pmax(1, abs(high))
      
      data.frame(
        Model = f$name,
        N_checked = length(ix),
        Max_relative_change = max(
          rel, na.rm = TRUE
        ),
        Stable = all(
          is.finite(rel) &
            rel <= cfg$capped_quad_tol
        )
      )
    })
  )
  
  write.csv(
    quad_checks,
    output_file(
      "CoverageLimited_StopLoss_Quadrature_Check.csv"
    ),
    row.names = FALSE
  )
  
  if (any(!quad_checks$Stable)) {
    warning(
      "Coverage-limited quadrature unstable for: ",
      paste(
        quad_checks$Model[!quad_checks$Stable],
        collapse = ", "
      ),
      ". Increase cfg$premium_gl_n"
    )
  }
}

# ================================================================
# Unlimited predictions (supplementary)
# ================================================================

if (cfg$compute_unlimited) {
  
  cat(
    "Predicting unlimited severity payments ",
    "with moment and convergence checks...\n"
  )
  
  sev_unlimited_pred <- lapply(
    severity_fits,
    function(f) {
      conditional_paid_unlimited(
        f, Dtheta_for(f), eval_dat
      )
    }
  )
  
} else {
  
  sev_unlimited_pred <- lapply(
    severity_fits,
    function(f) {
      rep(NA_real_, nrow(eval_pol))
    }
  )
}

# ================================================================
# Pure-premium model combinations
#
# XR: marginal / exogenous
# Case 1: XR1--3 (exogenous frequency and severity)
# Case 4c/4d: ER1--3 (both components dependent on deductible)
# ================================================================

model_pairs <- do.call(rbind, list(
  data.frame(Case = "1", Framework = "XR_marginal",
             Scheme = paste0("XR", 1:3),
             Frequency = c("FEM1", "FEM1", "FEM2_SEM3"),
             Severity = c("SEM1", "SEM2", "SEM3")),
  data.frame(Case = "4c", Framework = "Case4_both_cont",
             Scheme = paste0("ER", 1:3),
             Frequency = c("FCM1", "FCM1", "FCM2"),
             Severity = c("SCM1", "SCM2", "SCM3")),
  data.frame(Case = "4d", Framework = "Case4_both_disc",
             Scheme = paste0("ER", 1:3),
             Frequency = c("FDM1", "FDM1", "FDM2"),
             Severity = c("SDM1", "SDM2", "SDM3"))
))
stopifnot(nrow(model_pairs) == 9L,
          identical(as.character(unique(model_pairs$Case)), c("1", "4c", "4d")))

model_pairs$Model <- paste(
  model_pairs$Framework,
  model_pairs$Scheme,
  sep = "__"
)

write.csv(
  model_pairs,
  output_file("Premium_Model_Pairs.csv"),
  row.names = FALSE
)

# ================================================================
# Observed aggregate-loss metrics
# MAE_obs and portfolio bias use all finite predictions and observations.
# RMSRE_obs = sqrt(mean(((pred - obs) / obs)^2)) is evaluated ONLY
# when obs > 0, since relative error is undefined when obs == 0.
# N_RMSRE records how many policy-years enter this metric.
# NOTE: RMSRE_obs is conditional on positive observed payments, and is
# not a proper accuracy score for the mean premium over the full portfolio.
# ================================================================

observed_metrics <- function(pred, obs, n_total = length(obs)) {

  ok <- is.finite(pred) & is.finite(obs)
  n_ok <- sum(ok)

  if (!n_ok) {
    return(c(
      N = 0,
      N_RMSRE = 0,
      MAE_obs = NA_real_,
      RMSRE_obs = NA_real_,
      RB_obs_pct = NA_real_,
      Predicted_total = NA_real_,
      Observed_total = NA_real_
    ))
  }

  pr <- pred[ok]
  ob <- obs[ok]

  pos <- ob > 0
  n_pos <- sum(pos)

  rmsre <- if (n_pos > 0L) {
    sqrt(mean(((pr[pos] - ob[pos]) / ob[pos])^2))
  } else {
    NA_real_
  }

  c(
    N = n_ok,
    N_RMSRE = n_pos,
    MAE_obs = mean(abs(pr - ob)),
    RMSRE_obs = rmsre,

    RB_obs_pct = if (sum(ob) > 0) {
      100 * (sum(pr) - sum(ob)) / sum(ob)
    } else {
      NA_real_
    },

    Predicted_total = sum(pr),
    Observed_total = sum(ob)
  )
}

# ================================================================
# Compare predicted and observed aggregate insured loss
# ================================================================

raw_predictions <- list()
premium_rows <- list()
premium_d_rows <- list()

for (i in seq_len(nrow(model_pairs))) {
  
  pp <- model_pairs[i, ]
  
  fn <- pp$Frequency
  sn <- pp$Severity
  
  pred_lim <- as.numeric(
    freq_mean_pred[[fn]] *
      sev_limited_pred[[sn]]
  )
  
  pred_raw <- as.numeric(
    freq_mean_pred[[fn]] *
      sev_unlimited_pred[[sn]]
  )
  
  # Force exact zero even if an intermediate frequency value
  # becomes nonfinite: a contract with C_i <= D_i pays zero.
  pred_lim[zero_limit] <- 0
  
  raw_predictions[[pp$Model]] <- list(
    Limited = pred_lim,
    Unlimited = pred_raw
  )
  
  for (metric_source in c("Limited", "Unlimited")) {
    
    pred <- if (metric_source == "Limited") {
      pred_lim
    } else {
      pred_raw
    }
    
    obs <- if (metric_source == "Limited") {
      eval_pol$S_limited
    } else {
      eval_pol$S_raw
    }
    
    v <- observed_metrics(pred, obs)
    
    premium_rows[[
      length(premium_rows) + 1L
    ]] <- data.frame(
      Case = pp$Case,
      Framework = pp$Framework,
      Scheme = pp$Scheme,
      Model = pp$Model,
      Frequency = fn,
      Severity = sn,
      Payment_definition = metric_source,
      N_all = nrow(eval_pol),
      as.list(v),
      Sample = if (length(cfg$holdout_years)) {
        "Held_out_years"
      } else {
        "In_sample"
      }
    )
    
    for (d in cfg$ded_levels) {
      
      ii <- which(eval_pol$Deduct == d)
      
      if (!length(ii))
        next
      
      vd <- observed_metrics(
        pred[ii],
        obs[ii]
      )
      
      premium_d_rows[[
        length(premium_d_rows) + 1L
      ]] <- data.frame(
        Case = pp$Case,
        Framework = pp$Framework,
        Scheme = pp$Scheme,
        Model = pp$Model,
        Frequency = fn,
        Severity = sn,
        Payment_definition = metric_source,
        Deductible = d,
        N_all = length(ii),
        as.list(vd),
        Sample = if (length(cfg$holdout_years)) {
          "Held_out_years"
        } else {
          "In_sample"
        }
      )
    }
  }
}

premium_results <- do.call(
  rbind,
  premium_rows
)

premium_by_deductible <- do.call(
  rbind,
  premium_d_rows
)

write.csv(
  premium_results,
  output_file("Premium_Performance_Overall.csv"),
  row.names = FALSE
)

write.csv(
  premium_by_deductible,
  output_file("Premium_Performance_ByDeductible.csv"),
  row.names = FALSE
)

saveRDS(
  raw_predictions,
  output_file("Premium_Predictions_ByPolicy.rds")
)

# One row per policy-year x scheme, explicitly labeled by Case.
train_predictions_long <- do.call(rbind, lapply(seq_len(nrow(model_pairs)), function(i) {
  pp <- model_pairs[i, , drop = FALSE]
  data.frame(
    Case = pp$Case, Framework = pp$Framework, Scheme = pp$Scheme,
    Model = pp$Model, PolicyNum = eval_pol$PolicyNum, Year = eval_pol$Year,
    Deductible = eval_pol$Deduct,
    Observed_Aggregate_Limited = eval_pol$S_limited,
    Premium_Limited = raw_predictions[[pp$Model]]$Limited,
    stringsAsFactors = FALSE
  )
}))
write.csv(train_predictions_long,
          output_file("Premium_Predictions_ByPolicy_Long.csv"), row.names = FALSE)

# Policy-level comparison and coverage-limit audit
write.csv(
  data.frame(
    PolicyNum = eval_pol$PolicyNum,
    Year = eval_pol$Year,
    Deductible = eval_pol$Deduct,
    Coverage = eval_pol$Coverage,
    
    CoverageLimit_Dollars =
      eval_pol$CoverageLimit_Dollars,
    
    MaxInsurerPayment_PerOccurrence_Dollars =
      pmax(
        0,
        eval_pol$CoverageLimit_Dollars -
          eval_pol$Deduct
      ),
    
    CoverageBelowOrEqualDeductible =
      eval_pol$CoverageBelowOrEqualDeductible,
    
    Observed_Aggregate_Raw =
      eval_pol$S_raw,
    
    Observed_Aggregate_Limited =
      eval_pol$S_limited
  ),
  
  output_file("Premium_Evaluation_PolicyYears.csv"),
  row.names = FALSE
)

# Explicit check: all LIMITED premiums must be zero
# for policies with C_i <= D_i.
stopifnot(all(vapply(
  raw_predictions,
  function(z) all(z$Limited[zero_limit] == 0),
  logical(1)
)))

cat("\nCoverage-limited premium results:\n")

print(
  premium_results[
    premium_results$Payment_definition == "Limited",
    c(
      "Case", "Framework", "Scheme", "Frequency",
      "Severity", "N", "N_RMSRE", "MAE_obs",
      "RMSRE_obs", "RB_obs_pct"
    )
  ],
  row.names = FALSE,
  digits = 4
)


# ================================================================
# Premium performance plots by deductible level
# ================================================================

plot_premium_metric <- function(
    metric = "MAE_obs",
    payment_definition = "Limited") {
  
  plot_groups <- c("XR_marginal", "Case4_both_cont", "Case4_both_disc")
  
  subsetd <- premium_by_deductible[
    premium_by_deductible$Payment_definition ==
      payment_definition &
      premium_by_deductible$Framework %in% plot_groups,
  ]
  
  ran <- range(
    subsetd[[metric]],
    finite = TRUE
  )
  
  if (!all(is.finite(ran)))
    return(invisible(NULL))
  
  is_bias <- identical(metric, "RB_obs_pct")
  
  if (is_bias) {
    
    bound <- max(abs(ran), 1)
    ylim <- c(-bound, bound) * 1.05
    
  } else {
    
    positive <- subsetd[[metric]][
      is.finite(subsetd[[metric]]) &
        subsetd[[metric]] > 0
    ]
    
    if (!length(positive))
      return(invisible(NULL))
    
    ylim <- range(positive) * c(0.8, 1.25)
  }
  
  fname <- paste0(
    "Premium_ByDeductible_",
    metric, "_",
    payment_definition,
    ".pdf"
  )
  
  grDevices::pdf(
    output_file(fname),
    width = 11,
    height = 8
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  
  par(
    mfrow = c(3, 2),
    mar = c(4, 4.5, 2.5, 1)
  )
  
  cols <- c(
    "firebrick3",
    "steelblue4",
    "darkgreen"
  )
  
  for (group in plot_groups) {
    
    plot(
      NA,
      xlim = range(cfg$ded_levels),
      ylim = ylim,
      log = if (is_bias) "x" else "xy",
      xaxt = "n",
      xlab = "Deductible",
      ylab = if (is_bias) {
        metric
      } else {
        paste0(metric, " (log scale)")
      },
      main = group
    )
    
    axis(
      1,
      at = cfg$ded_levels,
      labels = format(
        cfg$ded_levels,
        big.mark = ","
      ),
      las = 2,
      cex.axis = 0.65
    )
    
    if (is_bias) {
      abline(
        h = 0,
        lty = 3,
        col = "grey50"
      )
    }
    
    for (i in 1:3) {
      
      scheme <- paste0(
        if (group == "XR_marginal") "XR" else "ER",
        i
      )
      
      sub <- subsetd[
        subsetd$Framework == group &
          subsetd$Scheme == scheme,
      ]
      
      if (nrow(sub)) {
        lines(
          sub$Deductible,
          sub[[metric]],
          col = cols[i],
          lty = i,
          pch = 15 + i,
          type = "b"
        )
      }
    }
    
    legend(
      "topleft",
      legend = if (group == "XR_marginal") {
        paste0("XR", 1:3)
      } else {
        paste0("ER", 1:3)
      },
      col = cols,
      lty = 1:3,
      pch = 16:18,
      bty = "n",
      cex = 0.72,
      bg = "white"
    )
  }
  
  plot.new()
  
  text(
    0.5,
    0.6,
    labels = paste0(
      "Same evaluation sample; ",
      payment_definition,
      " payments\n",
      "Policy ground-up occurrence limit: ",
      "C_i = Coverage_i x ",
      format(
        cfg$coverage_to_dollars,
        scientific = TRUE
      ),
      " dollars"
    )
  )
  
  invisible(NULL)
}

for (mm in c(
  "MAE_obs",
  "RMSRE_obs",
  "RB_obs_pct"
)) {
  plot_premium_metric(
    mm,
    "Limited"
  )
}

# ================================================================
# Separate primary and supplementary premium results
# ================================================================

write.csv(
  premium_results[
    premium_results$Payment_definition == "Unlimited",
  ],
  output_file(
    "Premium_Performance_Unlimited_Supplementary.csv"
  ),
  row.names = FALSE
)

write.csv(
  premium_results[
    premium_results$Payment_definition == "Limited",
  ],
  output_file(
    "Premium_Performance_CoverageLimited_Primary.csv"
  ),
  row.names = FALSE
)

# ================================================================
# Optional severity QQ/PIT diagnostics
# ================================================================

if (cfg$run_qq) {
  
  cat(
    "\n==== Appendix: ground-up QQ/PIT diagnostics ====\n"
  )
  
  cdf_train <- function(fit) {
    
    theta_D <- Dtheta_for(fit)
    
    conditional_Y_cdf(
      fit,
      theta_D,
      dat_clm,
      dat_clm$ystar
    )
  }
  
  all_fit <- c(
    list(
      SEM1 = SEM1,
      SEM2 = SEM2,
      SEM3 = SEM3
    ),
    list(
      SCM1 = SCM1,
      SCM2 = SCM2,
      SCM3 = SCM3
    ),
    list(
      SDM1 = SDM1,
      SDM2 = SDM2,
      SDM3 = SDM3
    )
  )
  
  qq_list <- lapply(
    all_fit,
    cdf_train
  )
  
  qq_lim <- c(-5.5, 5.5)
  
  grDevices::pdf(
    output_file(
      "QQ_Severity_AllModels_UnifiedAxes.pdf"
    ),
    width = 11,
    height = 11
  )
  
  par(
    mfrow = c(3, 3),
    pty = "s",
    mar = c(4, 4, 2.3, 0.8)
  )
  
  for (nm in names(qq_list)) {
    
    u <- qq_list[[nm]]
    z <- qnorm(clamp01(u))
    
    qqnorm(
      z,
      main = nm,
      pch = 20,
      cex = 0.3,
      col = "steelblue4",
      xlim = qq_lim,
      ylim = qq_lim,
      xlab = "Theoretical standard-normal quantile",
      ylab = "Fitted PIT normal score"
    )
    
    abline(
      0, 1,
      col = "firebrick3",
      lwd = 1.5
    )
  }
  
  grDevices::dev.off()
}

# ================================================================
# Optional frequency goodness-of-fit diagnostics
# ================================================================

if (cfg$run_gof) {
  
  cat(
    "\n==== Appendix: observed versus fitted frequency counts ====\n"
  )
  
  nmax <- cfg$nmax_gof
  
  empirical_count <- function(n) {
    
    tab <- tabulate(
      pmin(as.integer(n), nmax) + 1L,
      nbins = nmax + 1L
    )
    
    names(tab) <- c(
      as.character(0:(nmax - 1L)),
      paste0(">=", nmax)
    )
    
    tab
  }
  
  freq_pmf_for_fit <- function(
    fit, prev, scenario,
    severity_spec = NULL) {
    
    ss <- if (is.null(severity_spec)) {
      fit$spec$sev_spec
    } else {
      severity_spec
    }
    
    res <- matrix(
      NA_real_,
      nrow = nrow(X_pol),
      ncol = nmax
    )
    
    for (j in 0:(nmax - 1L)) {
      
      d <- dat_pol
      
      d$n_cen <- rep(
        j, nrow(X_pol)
      )
      
      d$n_tru <- rep(
        j, nrow(X_pol)
      )
      
      res[, j + 1L] <- exp(
        freq_loglik_i(
          fit$par,
          prev,
          scenario,
          fit$spec$dedtype,
          fit$spec$fam,
          ss,
          d
        )
      )
    }
    
    E <- c(
      colSums(res),
      sum(pmax(
        0,
        1 - rowSums(res)
      ))
    )
    
    names(E) <- c(
      as.character(0:(nmax - 1L)),
      paste0(">=", nmax)
    )
    
    E
  }
  
  gof_specs <- list(
    
    FEM1_full = list(
      fit = FEM1,
      prev = list(),
      sc = "censored"
    ),
    
    FEM2_SEM3 = list(
      fit = FEM2_SEM3,
      prev = list(S = SEM3$par),
      sc = "truncated"
    ),
    
    FEM2G_SEM1 = list(
      fit = FEM2G_SEM1,
      prev = list(S = SEM1$par),
      sc = "truncated"
    ),
    
    FEM2G_cont = list(
      fit = FEM2G_cont,
      prev = list(
        D = fit_Dc$par,
        S = SCM3$par
      ),
      sc = "truncated"
    ),
    
    FEM2G_disc = list(
      fit = FEM2G_disc,
      prev = list(
        D = fit_Dd$par,
        S = SDM3$par
      ),
      sc = "truncated"
    ),
    
    FCM1_full = list(
      fit = FCM1,
      prev = list(D = fit_Dc$par),
      sc = "censored"
    ),
    
    FDM1_full = list(
      fit = FDM1,
      prev = list(D = fit_Dd$par),
      sc = "censored"
    ),
    
    FCM2_truncated = list(
      fit = FCM2,
      prev = list(
        D = fit_Dc$par,
        S = SCM3$par
      ),
      sc = "truncated"
    ),
    
    FDM2_truncated = list(
      fit = FDM2,
      prev = list(
        D = fit_Dd$par,
        S = SDM3$par
      ),
      sc = "truncated"
    )
  )
  
  gof_rows <- list()
  
  for (nm in names(gof_specs)) {
    
    z <- gof_specs[[nm]]
    
    predicted <- tryCatch(
      freq_pmf_for_fit(
        z$fit,
        z$prev,
        z$sc
      ),
      
      error = function(e) {
        
        warning(
          "Count GOF error in ",
          nm, ": ",
          conditionMessage(e)
        )
        
        rep(
          NA_real_,
          nmax + 1L
        )
      }
    )
    
    actual <- empirical_count(
      if (z$sc == "censored") {
        pol$Freq
      } else {
        pol$posFreq
      }
    )
    
    pearson <- if (
      any(predicted <= 0) ||
      any(!is.finite(predicted))
    ) {
      NA_real_
    } else {
      sum(
        (actual - predicted)^2 /
          predicted
      )
    }
    
    gof_rows[[nm]] <- data.frame(
      Model = nm,
      Count = names(actual),
      Observed = as.numeric(actual),
      Expected = as.numeric(predicted),
      Pearson_Chi2_descriptive = pearson
    )
  }
  
  gof_freq <- do.call(
    rbind,
    gof_rows
  )
  
  write.csv(
    gof_freq,
    output_file("Frequency_GOF_Appendix.csv"),
    row.names = FALSE
  )
  
  warning(
    "Pearson chi-square values are descriptive; ",
    "use grouped cells and adjusted degrees ",
    "of freedom for formal GOF inference"
  )
}


################################################################################
################################################################################
## APPEND AFTER Section 7, "Optional frequency goodness-of-fit diagnostics"
## Compatible with Empirical_2006_2009_Training_2010_OOS_Complete.R
##
## Comparison: SAME observed ground-up count N = pol$Freq for ALL models.
## Truncated-frequency fits are evaluated at their implied LATENT N distribution,
## NOT at their N* likelihood. This relies on cfg$latent_mu == "groundup".
## No re-estimation; descriptive in-sample GOF, not formal chi-square tests.
################################################################################

cfg$run_common_N_gof <- TRUE
if (isTRUE(cfg$run_common_N_gof)) {
  if (!identical(cfg$latent_mu, "groundup")) {
    stop("Common-N GOF requires cfg$latent_mu == 'groundup'. ",
         "Under the 'thinned' convention the meaning of the NB latent mean differs; ",
         "do not use this plug-in without deriving its latent-N distribution.")
  }
  needed <- c("pol", "X_pol", "dat_pol", "freq_loglik_i", "cfg",
              "FEM1", "FEM2_SEM3", "FCM1", "FCM2", "FDM1", "FDM2",
              "fit_Dc", "fit_Dd")
  absent <- needed[!vapply(needed, exists, logical(1), inherits = TRUE)]
  if (length(absent)) stop("Run the full empirical fitting code first: ",
                           paste(absent, collapse = ", "))
  
  # FEM2 uses the original SEM3 thinning-based fitted beta/r, but evaluated
  # as the ground-up NB model rather than the thinned N* distribution.
  # FCM2 and FDM2 likewise use the NB--copula joint model for latent N.
  common_fits <- list(
    FEM1 = list(fit = FEM1, prev = list()),
    FEM2 = list(fit = FEM2_SEM3, prev = list()),
    FCM1 = list(fit = FCM1, prev = list(D = fit_Dc$par)),
    FCM2 = list(fit = FCM2, prev = list(D = fit_Dc$par)),
    FDM1 = list(fit = FDM1, prev = list(D = fit_Dd$par)),
    FDM2 = list(fit = FDM2, prev = list(D = fit_Dd$par))
  )
  
  nmax <- as.integer(cfg$nmax_gof)
  if (length(nmax) != 1L || is.na(nmax) || nmax < 1L)
    stop("cfg$nmax_gof must be a positive integer")
  n <- nrow(X_pol)
  if (!identical(as.integer(length(pol$Freq)), as.integer(n)))
    stop("pol$Freq and X_pol have different numbers of rows")
  if (anyNA(pol$Freq) || any(pol$Freq < 0) || any(pol$Freq != floor(pol$Freq)))
    stop("Ground-up frequency contains missing/negative/noninteger values")
  
  groups <- c(as.character(0:(nmax - 1L)), paste0(">=", nmax))
  observed <- tabulate(pmin(as.integer(pol$Freq), nmax) + 1L,
                       nbins = nmax + 1L)
  names(observed) <- groups
  
  latent_count_probs <- function(z) {
    fit <- z$fit
    if (is.null(fit$spec$dedtype) || is.null(fit$spec$fam))
      stop("Unexpected fit specification for ", fit$name)
    result <- matrix(NA_real_, nrow = n, ncol = nmax + 1L)
    for (k in 0:(nmax - 1L)) {
      d <- dat_pol
      d$n_cen <- rep.int(k, n)
      # IMPORTANT: "censored" selects the latent NB/copula frequency PMF;
      # it does NOT imply that this fit was estimated using complete counts.
      lp <- freq_loglik_i(fit$par, z$prev,
                          scenario = "censored",
                          dedtype = fit$spec$dedtype,
                          fam = fit$spec$fam,
                          sev_spec = NULL,
                          dat = d)
      if (length(lp) != n || any(!is.finite(lp)))
        stop("Nonfinite or malformed latent-N log probabilities for ", fit$name)
      result[, k + 1L] <- exp(lp)
    }
    sum_head <- rowSums(result[, seq_len(nmax), drop = FALSE])
    if (any(!is.finite(sum_head)) || any(sum_head > 1 + 1e-7))
      stop("Invalid latent-N probabilities for ", fit$name,
           ": sum of probabilities for 0,...,", nmax-1L, " exceeds 1")
    result[, nmax + 1L] <- pmax(0, 1 - sum_head)
    colnames(result) <- groups
    result
  }
  
  detail <- lapply(names(common_fits), function(model) {
    prob <- latent_count_probs(common_fits[[model]])
    expected <- colSums(prob)
    if (any(!is.finite(expected)) || any(expected <= 0)) {
      chi <- NA_real_
    } else {
      chi <- sum((as.numeric(observed) - expected)^2 / expected)
    }
    data.frame(
      Model = model,
      Count = groups,
      Observed = as.numeric(observed),
      Expected = as.numeric(expected),
      Pearson_Chi2_descriptive = chi,
      stringsAsFactors = FALSE
    )
  })
  gof_common_N <- do.call(rbind, detail)
  rownames(gof_common_N) <- NULL
  gof_common_summary <- unique(gof_common_N[
    c("Model", "Pearson_Chi2_descriptive")])
  gof_common_summary$Total_Observed <- sum(observed)
  gof_common_summary$Total_Expected <- vapply(
    gof_common_summary$Model,
    function(m) sum(gof_common_N$Expected[gof_common_N$Model == m]),
    numeric(1)
  )
  
  # Existing code's output_file() is used to respect cfg$output_dir.
  write.csv(gof_common_N, output_file("Frequency_GOF_CommonGroundupN.csv"),
            row.names = FALSE)
  write.csv(gof_common_summary,
            output_file("Frequency_GOF_CommonGroundupN_Summary.csv"),
            row.names = FALSE)
  
  # Overleaf-ready table. Same observed N for all six fitted models.
  latex_num <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "--")
  model_order <- names(common_fits)
  row_tex <- vapply(groups, function(g) {
    values <- vapply(model_order, function(m) {
      latex_num(gof_common_N$Expected[
        gof_common_N$Model == m & gof_common_N$Count == g])
    }, character(1))
    paste0(paste(c(if (grepl("^>=", g)) paste0("$\\geq", nmax, "$") else g,
                   as.character(observed[[g]]), values), collapse = " & "),
           " \\\\")
  }, character(1))
  chis <- vapply(model_order, function(m) {
    latex_num(gof_common_summary$Pearson_Chi2_descriptive[
      gof_common_summary$Model == m])
  }, character(1))
  table_tex <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\color{blue}",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrrr}",
    "\\toprule",
    "Count & Empirical & FEM1 & FEM2 & FCM1 & FCM2 & FDM1 & FDM2 \\\\",
    "\\midrule",
    row_tex,
    "\\midrule",
    paste0(paste(c("$\\chi^2$--stat", "--", chis), collapse = " & "),
           " \\\\"),
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    paste0("\\caption{Descriptive goodness-of-fit for the ground-up claim count $N$ ",
           "for all frequency models. Expected counts for FEM2, FCM2, and FDM2 ",
           "are computed from their fitted latent ground-up frequency distributions ",
           "rather than their thinned observed-count distributions. The final ",
           "category combines counts of at least ", nmax, ".}"),
    "\\label{tab:FreqcopulaCommonN}",
    "\\end{table}"
  )
  writeLines(table_tex, output_file("Table_Frequency_GOF_CommonGroundupN.tex"))
  
  cat("\n==== Common ground-up-N GOF: same empirical N for all fits ====\n")
  print(gof_common_summary, row.names = FALSE)
  cat("Wrote Frequency_GOF_CommonGroundupN.csv and ",
      "Table_Frequency_GOF_CommonGroundupN.tex\n")
}

################################################################################
## 8. MANUSCRIPT REPORTING: COPULA TABLES, PREMIUM SCORES AND ROBUSTNESS
##
## Uses fitted objects from Sections 0--7 above. No changes to the likelihood
## formulas or contractual payment definitions. Convergence-aware AIC selection
## may change a selected family if the previous winner did not converge.
## Table 1: scenario-specific SCM/SDM copula selection + candidate diagnostics
## Table 2: FCM/FDM copula selection + candidate diagnostics
## Table 3: coverage-limited premium scoring on a COMMON policy-year sample
## Figure 1: deductible-level O/P ratios and cluster bootstrap intervals
## Figure 2: observed/predicted severity across baseline premium risk groups
## Appendix: positive-payment RMSRE, unlimited payments, common-family
##           severity sensitivity, model contrasts and numerical diagnostics.
##
## These are IN-SAMPLE DESCRIPTIVE results unless cfg$holdout_years is set.
## Bootstrap intervals condition on the fitted premiums; refitting the entire
## model for each bootstrap replication would be needed for full model uncertainty.
################################################################################

manuscript_report <- NULL

if (isTRUE(cfg$reporting_enabled)) {
  cat("\n==== 6. Manuscript empirical reporting ===\n")

  report_dir <- output_file("manuscript_reporting")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  report_file <- function(name) file.path(report_dir, name)

  main_pairs <- model_pairs  # exactly nine Case 1, 4c, 4d premium schemes
  main_names <- main_pairs$Model
  observed_limited <- as.numeric(eval_pol$S_limited)
  nn <- nrow(eval_pol)

  # A single evaluation sample for ALL headline models, so rankings do not
  # change simply because numerical failures remove different policy-years.
  valid_common <- is.finite(observed_limited) & observed_limited >= 0
  invalid_by_model <- data.frame(
    Model = main_names,
    Nonfinite_or_negative_predictions = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (jj in seq_along(main_names)) {
    p <- raw_predictions[[main_names[jj]]]$Limited
    if (length(p) != nn) stop("Prediction length mismatch: ", main_names[jj])
    ok <- is.finite(p) & p >= 0
    invalid_by_model$Nonfinite_or_negative_predictions[jj] <- sum(!ok)
    valid_common <- valid_common & ok
  }
  write.csv(invalid_by_model,
            report_file("Appendix_Prediction_Completeness.csv"),
            row.names = FALSE)
  n_common <- sum(valid_common)
  if (n_common == 0L) {
    stop("No common valid evaluation policy-years across main model pairs. ",
         "Check Appendix_Prediction_Completeness.csv and conditional_N_mean().")
  }
  if (n_common < nn) {
    warning("Main comparison uses ", n_common, " / ", nn,
            " policy-years with valid predictions for EVERY main model.")
  }
  common_ids <- which(valid_common)
  sample_label <- if (length(cfg$holdout_years)) {
    "Held-out evaluation (models trained on other years)"
  } else {
    "In-sample descriptive evaluation"
  }

  # -------------------- core scoring rules ----------------------------------
  op_ratio <- function(y, p) {
    sy <- sum(y); sp <- sum(p)
    if (!is.finite(sy) || !is.finite(sp)) return(NA_real_)
    if (sp > 0) return(sy / sp)
    if (sy > 0) return(Inf)
    NA_real_  # zero / zero provides no information about calibration
  }

  poisson_deviance <- function(y, p) {
    if (length(y) != length(p) ||
        any(!is.finite(y)) || any(!is.finite(p)) ||
        any(y < 0) || any(p < 0))
      stop("Invalid observed or predicted losses in Poisson deviance")
    if (!length(y)) return(NA_real_)
    ans <- numeric(length(y))
    zero_obs <- y == 0
    ans[zero_obs] <- 2 * p[zero_obs]
    contradiction <- y > 0 & p == 0
    ans[contradiction] <- Inf
    good <- y > 0 & p > 0
    ans[good] <- 2 * (y[good] * log(y[good] / p[good]) -
                        (y[good] - p[good]))
    # finite precision can create tiny negative contributions; fix rounding only
    ans[is.finite(ans)] <- pmax(0, ans[is.finite(ans)])
    mean(ans)
  }

  # Tiny deterministic tests for zero-observation and zero-prediction cases.
  stopifnot(
    isTRUE(all.equal(poisson_deviance(c(0, 12), c(0, 12)), 0)),
    isTRUE(all.equal(poisson_deviance(0, 3), 6)),
    is.infinite(poisson_deviance(2, 0))
  )

  report_scores <- function(pred, obs) {
    stopifnot(length(pred) == length(obs))
    ok <- is.finite(pred) & is.finite(obs) & pred >= 0 & obs >= 0
    p <- as.numeric(pred[ok]); y <- as.numeric(obs[ok])
    n <- length(y)
    if (!n) {
      return(data.frame(
        N = 0L, N_positive = 0L,
        Observed_total = NA_real_, Predicted_total = NA_real_,
        Observed_mean = NA_real_, Predicted_mean = NA_real_,
        OP_ratio = NA_real_, RelBias_pct = NA_real_,
        Mean_Poisson_Deviance = NA_real_,
        MAE_obs = NA_real_, RMSE_obs = NA_real_,
        RMSRE_positive = NA_real_
      ))
    }
    positive <- y > 0
    data.frame(
      N = n,
      N_positive = sum(positive),
      Observed_total = sum(y),
      Predicted_total = sum(p),
      Observed_mean = mean(y),
      Predicted_mean = mean(p),
      OP_ratio = op_ratio(y, p),
      RelBias_pct = if (sum(y) > 0) {
        100 * (sum(p) - sum(y)) / sum(y)
      } else {
        NA_real_
      },
      Mean_Poisson_Deviance = poisson_deviance(y, p),
      MAE_obs = mean(abs(p - y)),
      RMSE_obs = sqrt(mean((p - y)^2)),
      RMSRE_positive = if (any(positive)) {
        sqrt(mean(((p[positive] - y[positive]) / y[positive])^2))
      } else {
        NA_real_
      }
    )
  }

  # -------------------- Tables 1 and 2: copula selection --------------------
  severity_models <- list(
    SCM1 = SCM1, SCM2 = SCM2, SCM3 = SCM3,
    SDM1 = SDM1, SDM2 = SDM2, SDM3 = SDM3
  )
  frequency_models <- list(
    FCM1 = FCM1, FCM2 = FCM2, FDM1 = FDM1, FDM2 = FDM2
  )
  key_to_model <- c(
    sev_cont_complete = "SCM1", sev_cont_censored = "SCM2",
    sev_cont_truncated = "SCM3", sev_disc_complete = "SDM1",
    sev_disc_censored = "SDM2", sev_disc_truncated = "SDM3",
    freq_cont_censored = "FCM1", freq_cont_truncated = "FCM2",
    freq_disc_censored = "FDM1", freq_disc_truncated = "FDM2"
  )

  copula_report <- function(keys, selected_fits) {
    candidate_parts <- vector("list", length(keys))
    selected_parts <- vector("list", length(keys))
    for (ii in seq_along(keys)) {
      key <- keys[ii]
      if (is.null(sel[[key]])) stop("Missing cached selection: ", key)
      spl <- strsplit(key, "_", fixed = TRUE)[[1L]]
      model_name <- unname(key_to_model[key])
      fit_sel <- selected_fits[[model_name]]
      ff <- sel[[key]]$fits
      nr <- length(ff)
      if (!nr) stop("Empty candidate list for ", key)
      aic <- vapply(ff, function(f)
        if (is.null(f)) NA_real_ else as.numeric(f$AIC), numeric(1))
      loglik <- vapply(ff, function(f)
        if (is.null(f)) NA_real_ else as.numeric(f$loglik), numeric(1))
      bic <- vapply(ff, function(f)
        if (is.null(f)) NA_real_ else as.numeric(f$BIC), numeric(1))
      conv <- vapply(ff, function(f)
        if (is.null(f) || is.null(f$convergence)) NA_integer_
        else as.integer(f$convergence), integer(1))
      tau <- vapply(ff, function(f) {
        if (is.null(f)) return(NA_real_)
        if (f$spec$fam == 0L) return(0)
        if (!"tau" %in% names(f$par)) return(NA_real_)
        as.numeric(f$par["tau"])
      }, numeric(1))
      nu <- vapply(ff, function(f)
        if (is.null(f) || !"nu" %in% names(f$par)) NA_real_
        else as.numeric(f$par["nu"]), numeric(1))
      n_used <- vapply(ff, function(f)
        if (is.null(f)) NA_integer_ else as.integer(f$n), integer(1))
      converged <- !is.na(conv) & conv == 0L & is.finite(aic)
      delta <- rep(NA_real_, nr)
      wt <- rep(NA_real_, nr)
      if (any(converged)) {
        delta[converged] <- aic[converged] - min(aic[converged])
        ww <- exp(-0.5 * delta[converged])
        wt[converged] <- ww / sum(ww)
      }
      chosen <- fam_name(fit_sel$spec$fam)
      cc <- data.frame(
        Model = model_name,
        Deductible = spl[2L],
        Scenario = spl[3L],
        Specification = key,
        Candidate_family = names(ff),
        LogLik = loglik, AIC = aic, BIC = bic,
        Delta_AIC = delta, Akaike_weight = wt,
        Kendall_tau = tau, Student_t_df = nu,
        Convergence_code = conv, Converged = converged,
        N_likelihood_contributions = n_used,
        Selected = names(ff) == chosen,
        stringsAsFactors = FALSE
      )
      cc <- cc[order(cc$Delta_AIC, cc$AIC, na.last = TRUE),
               , drop = FALSE]
      candidate_parts[[ii]] <- cc
      best <- cc[cc$Selected, , drop = FALSE]
      if (nrow(best) != 1L)
        stop("Selected copula not represented exactly once: ", key)
      if (!best$Converged)
        warning("Selected family has a nonzero convergence code: ", key)
      selected_parts[[ii]] <- data.frame(
        Model = model_name,
        Deductible = spl[2L], Scenario = spl[3L],
        Selected_family = chosen,
        Kendall_tau = best$Kendall_tau,
        Student_t_df = best$Student_t_df,
        LogLik = best$LogLik, AIC = best$AIC,
        Delta_AIC = best$Delta_AIC,
        Akaike_weight = best$Akaike_weight,
        Second_best_Delta_AIC = if (sum(converged) > 1L) {
          sort(cc$Delta_AIC[is.finite(cc$Delta_AIC)])[2L]
        } else NA_real_,
        N_likelihood_contributions = best$N_likelihood_contributions,
        Convergence_code = best$Convergence_code,
        stringsAsFactors = FALSE
      )
    }
    list(selected = do.call(rbind, selected_parts),
         candidates = do.call(rbind, candidate_parts))
  }

  t1 <- copula_report(required_sev, severity_models)
  t2 <- copula_report(required_freq, frequency_models)
  write.csv(t1$selected, report_file("Table1_Severity_Copula_Selected.csv"),
            row.names = FALSE)
  write.csv(t1$candidates,
            report_file("Table1_Severity_Copula_AllCandidates.csv"),
            row.names = FALSE)
  write.csv(t2$selected, report_file("Table2_Frequency_Copula_Selected.csv"),
            row.names = FALSE)
  write.csv(t2$candidates,
            report_file("Table2_Frequency_Copula_AllCandidates.csv"),
            row.names = FALSE)

  # -------------------- Table 3: common-sample main scores ------------------
  main_scores <- do.call(rbind, lapply(seq_len(nrow(main_pairs)), function(i) {
    pp <- main_pairs[i, , drop = FALSE]
    score <- report_scores(
      raw_predictions[[pp$Model]]$Limited[valid_common],
      observed_limited[valid_common]
    )
    cbind(pp[, c("Case", "Framework", "Scheme", "Model", "Frequency", "Severity")],
          score)
  }))
  rownames(main_scores) <- NULL
  table3 <- main_scores[, c(
    "Case", "Framework", "Scheme", "Model", "Frequency", "Severity", "N",
    "OP_ratio", "RelBias_pct", "Mean_Poisson_Deviance",
    "MAE_obs", "RMSE_obs", "Observed_total", "Predicted_total"
  )]
  write.csv(table3, report_file("Table3_Limited_PurePremium_Main.csv"),
            row.names = FALSE)
  write.csv(main_scores,
            report_file("Appendix_Limited_PurePremium_AllMetrics.csv"),
            row.names = FALSE)

  # -------------------- Figure 1: deductible-level calibration -------------
  ded_parts <- list()
  k <- 0L
  for (i in seq_len(nrow(main_pairs))) {
    pp <- main_pairs[i, , drop = FALSE]
    prediction <- raw_predictions[[pp$Model]]$Limited
    for (d in cfg$ded_levels) {
      ix <- which(valid_common & eval_pol$Deduct == d)
      if (!length(ix)) next
      ss <- report_scores(prediction[ix], observed_limited[ix])
      k <- k + 1L
      ded_parts[[k]] <- cbind(
        pp[, c("Case", "Framework", "Scheme", "Model")],
        data.frame(Deductible = d), ss
      )
    }
  }
  ded_scores <- do.call(rbind, ded_parts)
  rownames(ded_scores) <- NULL

  # Bootstrap over policy IDs within each deductible group. The fitted model
  # is held fixed. Repeated policy-years from one policy form one cluster.
  cluster_ci <- function(y, p, ids, B, min_clusters) {
    sums <- rowsum(cbind(y = y, p = p), group = as.character(ids),
                   reorder = FALSE)
    K <- nrow(sums)
    if (K < min_clusters || B < 2L)
      return(c(N_clusters = K, OP_lower = NA_real_, OP_upper = NA_real_))
    v <- rep(NA_real_, B)
    for (bb in seq_len(B)) {
      sampled <- sample.int(K, K, replace = TRUE)
      v[bb] <- op_ratio(sum(sums[sampled, "y"]),
                         sum(sums[sampled, "p"]))
    }
    v <- v[is.finite(v)]
    if (length(v) < max(20L, ceiling(B / 2)))
      return(c(N_clusters = K, OP_lower = NA_real_, OP_upper = NA_real_))
    ci <- as.numeric(quantile(v, c(0.025, 0.975), names = FALSE))
    c(N_clusters = K, OP_lower = ci[1L], OP_upper = ci[2L])
  }

  set.seed(cfg$report_bootstrap_seed)
  ded_scores$N_clusters <- NA_integer_
  ded_scores$OP_lower <- NA_real_
  ded_scores$OP_upper <- NA_real_
  for (i in seq_len(nrow(ded_scores))) {
    d <- ded_scores$Deductible[i]
    nm <- ded_scores$Model[i]
    ix <- which(valid_common & eval_pol$Deduct == d)
    ci <- cluster_ci(observed_limited[ix],
                     raw_predictions[[nm]]$Limited[ix],
                     eval_pol$PolicyNum[ix],
                     cfg$report_bootstrap_reps,
                     cfg$report_min_bootstrap_clusters)
    ded_scores$N_clusters[i] <- as.integer(ci["N_clusters"])
    ded_scores$OP_lower[i] <- ci["OP_lower"]
    ded_scores$OP_upper[i] <- ci["OP_upper"]
  }
  write.csv(ded_scores,
            report_file("Figure1_Deductible_OP_and_95pctClusterCI.csv"),
            row.names = FALSE)

  groups <- unique(main_pairs$Framework)
  scheme_colors <- c("#21618C", "#D35400", "#148F77")
  grDevices::pdf(report_file("Figure1_Deductible_OP_95pctCI.pdf"),
                 width = 11, height = 8.5)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2), mar = c(4.0, 4.2, 2.7, 1),
      oma = c(1.2, 0, 1.2, 0))
  for (g in groups) {
    dtab <- ded_scores[ded_scores$Framework == g, , drop = FALSE]
    if (!nrow(dtab)) { plot.new(); next }
    yhigh <- cfg$report_op_plot_cap
    plot(NA_real_, NA_real_,
         xlim = range(cfg$ded_levels), ylim = c(0, yhigh), log = "x",
         xlab = "Deductible (USD)", ylab = "Observed / predicted",
         main = g, xaxt = "n")
    axis(1, at = cfg$ded_levels,
         labels = format(cfg$ded_levels, scientific = FALSE, big.mark = ","),
         las = 2, cex.axis = 0.58)
    abline(h = 1, lty = 2, col = "grey45")
    schemes <- unique(dtab$Scheme)
    for (s in seq_along(schemes)) {
      a <- dtab[dtab$Scheme == schemes[s], , drop = FALSE]
      a <- a[order(a$Deductible), , drop = FALSE]
      xx <- a$Deductible * (1 + (s - 2) * 0.025)
      yy <- pmin(a$OP_ratio, yhigh)
      yy[!is.finite(yy)] <- NA_real_
      lines(xx, yy, type = "b", pch = 15 + s, col = scheme_colors[s],
            lwd = 1.5)
      goodci <- is.finite(a$OP_lower) & is.finite(a$OP_upper)
      if (any(goodci)) {
        arrows(xx[goodci], pmax(0, pmin(a$OP_lower[goodci], yhigh)),
               xx[goodci], pmin(a$OP_upper[goodci], yhigh),
               angle = 90, code = 3, length = .025,
               col = scheme_colors[s])
      }
    }
    legend("topright", legend = schemes,
           col = scheme_colors[seq_along(schemes)],
           pch = 15 + seq_along(schemes), lty = 1, bty = "n",
           cex = 0.75)
  }
  plot.new()
  text(0.5, 0.68,
       "Common policy-year evaluation sample\n95% descriptive policy-cluster bootstrap intervals",
       cex = 0.95)
  text(0.5, 0.35,
       sprintf("Display truncated above O/P = %s; see CSV for full values",
               cfg$report_op_plot_cap), cex = 0.8)
  mtext(sample_label, outer = TRUE, side = 3, line = 0.05)
  par(oldpar)
  dev.off()

  # -------------------- Figure 2: baseline-risk grouped calibration --------
  # All model predictions use the SAME bins, formed from XR1's predicted
  # pure premium, thus observed payments and group composition are identical.
  baseline_model <- "XR_marginal__XR1"
  if (!baseline_model %in% main_names)
    stop("Baseline XR1 is unavailable for risk grouping")
  base_pred <- raw_predictions[[baseline_model]]$Limited[valid_common]
  G <- as.integer(cfg$report_risk_groups)
  if (G < 2L) stop("report_risk_groups must be >= 2")
  rr <- rank(base_pred, ties.method = "average")
  risk_group <- pmin(G, pmax(1L, ceiling(G * rr / length(rr))))
  group_parts <- list(); k <- 0L
  for (i in seq_len(nrow(main_pairs))) {
    pp <- main_pairs[i, , drop = FALSE]
    prediction <- raw_predictions[[pp$Model]]$Limited[valid_common]
    for (z in sort(unique(risk_group))) {
      ix <- which(risk_group == z)
      ss <- report_scores(prediction[ix], observed_limited[valid_common][ix])
      k <- k + 1L
      group_parts[[k]] <- cbind(
        pp[, c("Case", "Framework", "Scheme", "Model")],
        data.frame(Risk_group = z,
                   Baseline_mean_premium = mean(base_pred[ix])), ss
      )
    }
  }
  risk_scores <- do.call(rbind, group_parts)
  rownames(risk_scores) <- NULL
  write.csv(risk_scores,
            report_file("Figure2_BaselineRiskGroup_Calibration.csv"),
            row.names = FALSE)

  grDevices::pdf(report_file("Figure2_Observed_vs_Predicted_ByRiskGroup.pdf"),
                 width = 11, height = 8.5)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2), mar = c(4, 4.6, 2.5, 0.8),
      oma = c(0.8, 0, 1.2, 0))
  for (g in groups) {
    gg <- risk_scores[risk_scores$Framework == g, , drop = FALSE]
    if (!nrow(gg)) { plot.new(); next }
    yrs <- c(gg$Observed_mean, gg$Predicted_mean)
    yr <- range(yrs, finite = TRUE)
    if (!all(is.finite(yr))) yr <- c(0, 1)
    yr <- c(0, max(1, yr[2L] * 1.1))
    base <- gg[gg$Scheme == unique(gg$Scheme)[1L], , drop = FALSE]
    base <- base[order(base$Risk_group), , drop = FALSE]
    plot(base$Risk_group, base$Observed_mean,
         type = "b", pch = 19, lwd = 2, col = "#2C3E50",
         xlim = c(1, G), ylim = yr,
         xlab = "XR1 baseline premium risk group", ylab = "Mean aggregate payment (USD)",
         main = g)
    schemes <- unique(gg$Scheme)
    for (s in seq_along(schemes)) {
      a <- gg[gg$Scheme == schemes[s], , drop = FALSE]
      a <- a[order(a$Risk_group), , drop = FALSE]
      lines(a$Risk_group, a$Predicted_mean,
            type = "b", col = scheme_colors[s], lwd = 1.5,
            pch = 15 + s, lty = s)
    }
    legend("topleft", legend = c("Observed", schemes),
           col = c("#2C3E50", scheme_colors[seq_along(schemes)]),
           pch = c(19, 15 + seq_along(schemes)),
           lty = c(1, seq_along(schemes)), bty = "n", cex = .7)
  }
  plot.new()
  text(0.5, 0.5,
       "Fixed XR1 baseline risk groups\nIdentical observations and groups for every model\nTies are kept in the same group",
       cex = .94)
  mtext(sample_label, outer = TRUE, side = 3, line = 0.05)
  par(oldpar)
  dev.off()

  # -------------------- Appendix: paired Case 1 vs Case 4 contrasts ----------
  contrasts <- do.call(rbind, lapply(c("4c", "4d"), function(ccase) {
    do.call(rbind, lapply(1:3, function(scheme_idx) {
      base <- model_pairs[model_pairs$Case == "1" &
                            model_pairs$Scheme == paste0("XR", scheme_idx), , drop = FALSE]
      alt <- model_pairs[model_pairs$Case == ccase &
                           model_pairs$Scheme == paste0("ER", scheme_idx), , drop = FALSE]
      data.frame(Case = ccase, Scheme = scheme_idx,
                 Base_Model = base$Model,
                 Alternative_Model = alt$Model)
    }))
  }))
  contrast_parts <- lapply(seq_len(nrow(contrasts)), function(i) {
    z <- contrasts[i, , drop = FALSE]
    pa <- raw_predictions[[z$Base_Model]]$Limited
    pb <- raw_predictions[[z$Alternative_Model]]$Limited
    both <- valid_common & is.finite(pa) & is.finite(pb) & pa >= 0 & pb >= 0
    a <- report_scores(pa[both], observed_limited[both])
    b <- report_scores(pb[both], observed_limited[both])
    cbind(z, data.frame(
      N = a$N, Base_OP = a$OP_ratio, Alternative_OP = b$OP_ratio,
      Base_Poisson_Deviance = a$Mean_Poisson_Deviance,
      Alternative_Poisson_Deviance = b$Mean_Poisson_Deviance,
      Delta_Poisson_Deviance = b$Mean_Poisson_Deviance - a$Mean_Poisson_Deviance,
      Base_MAE = a$MAE_obs, Alternative_MAE = b$MAE_obs,
      Delta_MAE = b$MAE_obs - a$MAE_obs
    ))
  })
  contrast_scores <- do.call(rbind, contrast_parts)
  write.csv(contrast_scores,
            report_file("Appendix_Case1_vs_Case4_Contrasts.csv"), row.names = FALSE)

  # -------------------- Appendix: common-copula sensitivity -----------------
  # Re-use EXISTING AIC-selection candidate fits; no extra optimization.
  # To hold frequency constant, use FEM1 for all six candidate severity models.
  # This is a severity-component sensitivity, NOT an ER3 refit under alternative
  # thinning probabilities. Fits not converged are omitted and flagged.
  common_family_scores <- data.frame()
  if (isTRUE(cfg$report_common_copula_sensitivity)) {
    cat("Common-family severity sensitivity using cached candidates...\n")
    sens_parts <- list(); k <- 0L
    for (key in required_sev) {
      model_name <- unname(key_to_model[key])
      chosen_fit <- severity_models[[model_name]]
      candidate_fits <- sel[[key]]$fits
      chosen_family <- fam_name(chosen_fit$spec$fam)
      family_codes <- unique(c(chosen_fit$spec$fam,
                               as.integer(cfg$report_common_family_codes)))
      for (code in family_codes) {
        fam <- fam_name(code)
        f <- candidate_fits[[fam]]
        available <- !is.null(f) && is.finite(f$AIC) &&
          !is.null(f$convergence) && f$convergence == 0L
        if (!available) {
          k <- k + 1L
          sens_parts[[k]] <- data.frame(
            Model = model_name, Candidate_family = fam,
            AIC_selected_family = chosen_family,
            Used = FALSE, Reason = "Not fitted / nonconverged",
            report_scores(numeric(0), numeric(0))
          )
          next
        }
        pred <- as.numeric(freq_mean_pred$FEM1 *
          conditional_paid_limited(f, Dtheta_for(f), eval_dat))
        pred[zero_limit] <- 0
        use <- valid_common & is.finite(pred) & pred >= 0
        k <- k + 1L
        sens_parts[[k]] <- data.frame(
          Model = model_name, Candidate_family = fam,
          AIC_selected_family = chosen_family,
          Used = TRUE,
          Reason = "FEM1 held fixed (severity-only)",
          report_scores(pred[use], observed_limited[use])
        )
      }
    }
    common_family_scores <- do.call(rbind, sens_parts)
    write.csv(common_family_scores,
              report_file("Appendix_CommonCopula_SeverityOnly.csv"),
              row.names = FALSE)
  }

  # -------------------- Appendix: supplemental pre-existing outputs ----------
  write.csv(moments,
            report_file("Appendix_GB2_FirstMoment_Existence.csv"),
            row.names = FALSE)
  if (!is.null(quad_checks))
    write.csv(quad_checks,
              report_file("Appendix_CoverageLimited_Quadrature.csv"),
              row.names = FALSE)
  write.csv(premium_results[
    premium_results$Payment_definition == "Unlimited", , drop = FALSE],
    report_file("Appendix_Unlimited_Premiums.csv"), row.names = FALSE)
  write.csv(premium_results[
    premium_results$Payment_definition == "Limited", , drop = FALSE],
    report_file("Appendix_OriginalEvaluation_Metrics.csv"),
    row.names = FALSE)

  # Manuscript and implementation notes kept alongside tables and figures.
  notes <- c(
    "EMPIRICAL MANUSCRIPT REPORTING NOTES",
    paste("Sample:", sample_label),
    paste("Main headline sample:", n_common, "/", nn, "policy-years"),
    "",
    "Table 1: severity copula selection, one AIC ranking per deductible type",
    "  and observation scenario (complete, censored, truncated).",
    "Table 2: frequency copula selection, one ranking per deductible type",
    "  and observation scenario (censored, truncated).",
    "Do NOT compare AIC across different observed-data likelihoods or samples.",
    "Only converged candidates enter delta-AIC and Akaike weights.",
    "",
    "Table 3: all primary models are scored on the same valid policy-years.",
    "Limited means (min(Y,C)-D)_+ with the policy's assumed ground-up C.",
    "Mean Poisson deviance is a mean-oriented SCORE and does not require",
    "the aggregate payments themselves to have a Poisson distribution.",
    "Observed/predicted (O/P) near 1 is portfolio calibration, not proof",
    "of individual predictive discrimination.",
    "",
    "Figure 1: policy-level resampling within deductible groups. Fits are",
    "FIXED; intervals are descriptive and do NOT propagate parameter",
    "estimation uncertainty or account for in-sample model selection.",
    paste("Bootstrap replicates:", cfg$report_bootstrap_reps),
    paste("Minimum independent policy clusters:",
          cfg$report_min_bootstrap_clusters),
    "",
    "Figure 2: decile-like groups are based on the same XR1 premium and",
    "are identical for every model; ties remain grouped together.",
    "",
    "Appendix RMSRE_positive excludes policy-years with observed S=0.",
    "RMSE remains available and is not the same as RMSRE.",
    "Common-copula sensitivity holds FEM1 fixed, even for Scenario 3:",
    "it isolates the severity copula and does NOT refit thinning models.",
    "Unlimited premiums are supplementary; their finite mean is checked",
    "conservatively in the original code.",
    "",
    "IMPORTANT MODELING ASSUMPTIONS TO VALIDATE:",
    "1. Coverage is a per-occurrence ground-up limit in dollars after scaling.",
    "2. The claim data represent the observation mechanism in each likelihood.",
    "3. The frequency-severity factorization is substantively justified.",
    "4. Results are in-sample descriptive unless holdout_years is set."
  )
  writeLines(notes, report_file("README_Manuscript_Reporting.txt"))

  manuscript_report <- list(
    sample_label = sample_label,
    n_common = n_common,
    comparison_index = common_ids,
    table1 = t1$selected,
    table2 = t2$selected,
    table3 = table3,
    deductible_op = ded_scores,
    risk_group_calibration = risk_scores,
    component_contrasts = contrast_scores,
    common_family_sensitivity = common_family_scores,
    appendix_all_scores = main_scores
  )

  cat("\nManuscript Table 1 (severity copulas):\n")
  print(t1$selected[, c("Model", "Selected_family", "Kendall_tau",
                        "Delta_AIC", "AIC", "N_likelihood_contributions")],
        row.names = FALSE)
  cat("\nManuscript Table 2 (frequency copulas):\n")
  print(t2$selected[, c("Model", "Selected_family", "Kendall_tau",
                        "Delta_AIC", "AIC", "N_likelihood_contributions")],
        row.names = FALSE)
  cat("\nManuscript Table 3 (limited pure premiums):\n")
  print(table3[, c("Model", "N", "OP_ratio", "RelBias_pct",
                   "Mean_Poisson_Deviance", "MAE_obs", "RMSE_obs")],
        row.names = FALSE, digits = 4)
  cat("Reporting saved to: ", normalizePath(report_dir), "\n", sep = "")
}

# ================================================================
# Save final results
# ================================================================

results <- list(
  
  config = cfg,
  
  training_policy_years = nrow(pol),
  evaluation_policy_years = nrow(eval_pol),
  
  deductible = list(
    continuous = fit_Dc,
    discrete = fit_Dd
  ),
  
  marginal_severity = list(
    SEM1 = SEM1,
    SEM2 = SEM2,
    SEM3 = SEM3
  ),
  
  selected_severity = severity_fits,
  frequency = frequency_fits,
  
  copula_selection = sel,
  
  selected_severity_table =
    selected_severity,
  
  selected_frequency_table =
    selected_frequency,
  
  FEM_table = fem_table,
  mu_slopes = mu_slopes,
  
  moment_existence = moments,
  
  capped_quadrature_check = quad_checks,
  
  premium_model_pairs = model_pairs,
  
  premium_results = premium_results,
  
  premium_by_deductible =
    premium_by_deductible,
  
  godambe = sw,
  manuscript_reporting = manuscript_report
)

saveRDS(
  results,
  output_file(
    "Empirical_Revised_Results.rds"
  )
)

# Save complete R workspace
save.image(cfg$out_file)

cat(
  "\nSUCCESS: revised empirical pipeline completed.\n",
  "Output directory: ",
  normalizePath(cfg$output_dir),
  "\n",
  sep = ""
)

if (!length(cfg$holdout_years)) {
  cat(
    "WARNING: evaluation is IN SAMPLE. ",
    "The appended 2010 validations provide out-of-time performance.\n"
  )
}

cat(
  "NOTE: Limited payment = (min(Y,C_i)-D_i)_+\n",
  "where C_i = Coverage_i * ",
  cfg$coverage_to_dollars,
  " USD.\n",
  sep = ""
)

cat(
  "NOTE: C_i <= D_i policy-years have ",
  "zero coverage-limited pure premiums.\n"
)

cat(
  "NOTE: Unlimited premiums are supplementary ",
  "and require valid finite expectations.\n"
)

# ================================================================
# Continuous deductible-ratio diagnostics (NO quantile groups)
# Place after raw_predictions and model_pairs have been created.
# Uses every policy-year's actual D / Coverage (Coverage in raw data units).
# Required existing objects: eval_pol, model_pairs, raw_predictions,
#                            output_file().
# ================================================================

ratio <- eval_pol$Deduct / eval_pol$Coverage
observed <- eval_pol$S_limited

stopifnot(length(ratio) == nrow(eval_pol),
          length(observed) == nrow(eval_pol))
stopifnot(all(is.finite(ratio) & ratio > 0),
          all(is.finite(observed) & observed >= 0))

# 1) Policy-level diagnostics: no bins and no aggregation.
ratio_policy_rows <- lapply(seq_len(nrow(model_pairs)), function(k) {
  spec <- model_pairs[k, , drop = FALSE]
  pred <- as.numeric(raw_predictions[[spec$Model]]$Limited)
  stopifnot(length(pred) == length(ratio))
  ok <- is.finite(pred) & pred >= 0
  
  # Individual Poisson deviance contribution, for descriptive use.
  dev <- rep(NA_real_, length(pred))
  positive_pred <- ok & pred > 0
  zero_obs <- positive_pred & observed == 0
  positive_obs <- positive_pred & observed > 0
  dev[zero_obs] <- 2 * pred[zero_obs]
  dev[positive_obs] <- 2 * (
    observed[positive_obs] * log(observed[positive_obs] / pred[positive_obs]) -
      (observed[positive_obs] - pred[positive_obs])
  )
  dev[ok & pred == 0 & observed == 0] <- 0
  dev[ok & pred == 0 & observed > 0] <- Inf
  data.frame(
    Case = spec$Case,
    Framework = spec$Framework,
    Scheme = spec$Scheme,
    Model = spec$Model,
    PolicyNum = eval_pol$PolicyNum,
    Year = eval_pol$Year,
    Deductible = eval_pol$Deduct,
    Coverage = eval_pol$Coverage,
    DeductibleRatio = ratio,
    Observed = observed,
    Predicted = pred,
    AbsoluteError = abs(observed - pred),
    SquaredError = (observed - pred)^2,
    PoissonDeviance = dev,
    stringsAsFactors = FALSE
  )
})
ratio_policy <- do.call(rbind, ratio_policy_rows)
write.csv(ratio_policy,
          output_file("Premium_ByDeductibleRatio_PolicyLevel.csv"),
          row.names = FALSE)

# 2) Exact-ratio aggregates (only genuinely identical D/C values pooled).
# No deciles, cut(), quantiles, or arbitrary groups.
make_exact <- function(d) {
  split_idx <- split(seq_len(nrow(d)), d$DeductibleRatio)
  do.call(rbind, lapply(split_idx, function(idx) {
    a <- d[idx, , drop = FALSE]
    p <- sum(a$Predicted)
    o <- sum(a$Observed)
    data.frame(
      Model = a$Model[1],
      Framework = a$Framework[1],
      Scheme = a$Scheme[1],
      DeductibleRatio = a$DeductibleRatio[1],
      N = nrow(a),
      MAE = mean(a$AbsoluteError),
      RMSE = sqrt(mean(a$SquaredError)),
      MPD = mean(a$PoissonDeviance),
      OP_ratio = if (p > 0) o / p else NA_real_,
      ObservedTotal = o,
      PredictedTotal = p
    )
  }))
}
ratio_exact <- do.call(rbind,
                       lapply(split(ratio_policy, ratio_policy$Model), make_exact))
ratio_exact <- ratio_exact[order(ratio_exact$Framework,
                                 ratio_exact$Scheme,
                                 ratio_exact$DeductibleRatio), ]
row.names(ratio_exact) <- NULL
write.csv(ratio_exact,
          output_file("Premium_ByDeductibleRatio_ExactValues.csv"),
          row.names = FALSE)

# 3) Smooth conditional metrics across the continuous D/C axis.
# Gaussian-kernel weights in log10(D/C), not discrete risk groups.
# Bandwidth controls local smoothing and should be reported if published.
smooth_ratio <- function(d, grid, bandwidth) {
  x <- log10(d$DeductibleRatio)
  p <- d$Predicted
  o <- d$Observed
  out <- lapply(grid, function(r) {
    w <- exp(-0.5 * ((x - log10(r)) / bandwidth)^2)
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) return(NULL)
    w <- w / sw
    effective_n <- 1 / sum(w^2)
    pred_sum <- sum(w * p)
    obs_sum <- sum(w * o)
    good_dev <- is.finite(d$PoissonDeviance)
    
    data.frame(
      DeductibleRatio = r,
      EffectiveN = effective_n,
      MAE = sum(w * d$AbsoluteError),
      RMSE = sqrt(sum(w * d$SquaredError)),
      MPD = if (all(good_dev)) sum(w * d$PoissonDeviance)
      else if (any(good_dev))
        sum(w[good_dev] * d$PoissonDeviance[good_dev]) /
        sum(w[good_dev])
      else NA_real_,
      OP_ratio = if (pred_sum > 0) obs_sum / pred_sum else NA_real_
    )
  })
  do.call(rbind, out)
}

log_ratio <- log10(ratio)
spread <- IQR(log_ratio)
# Configurable bandwidth in log10-ratio units. No hard-coded group count.
ratio_bandwidth <- if (is.finite(spread) && spread > 0) 0.20 * spread else 0.20
ratio_grid <- 10^seq(min(log_ratio), max(log_ratio), length.out = 150L)

cont_schemes <- ratio_policy[
  ratio_policy$Framework %in% c("XR_marginal", "Case4_both_cont"), ]
cont_models <- split(cont_schemes, cont_schemes$Model)
ratio_curves <- do.call(rbind, lapply(cont_models, function(d) {
  z <- smooth_ratio(d, ratio_grid, ratio_bandwidth)
  z$Model <- d$Model[1]
  z$Framework <- d$Framework[1]
  z$Scheme <- d$Scheme[1]
  z
}))
row.names(ratio_curves) <- NULL
write.csv(ratio_curves,
          output_file("Premium_ByDeductibleRatio_Smoothed.csv"),
          row.names = FALSE)

# 4) Comparable XRk vs continuous ERk figures: one panel per scheme.
# O/P = 1 indicates local aggregate calibration.
plot_ratio_comparison <- function(metric) {
  stopifnot(metric %in% c("MAE", "RMSE", "MPD", "OP_ratio"))
  outfile <- output_file(paste0("Premium_ByDeductibleRatio_", metric,
                                "_XR_vs_ContinuousER.pdf"))
  grDevices::pdf(outfile, width = 9, height = 10)
  op <- par(no.readonly = TRUE)
  on.exit({par(op); grDevices::dev.off()}, add = TRUE)
  par(mfrow = c(3, 1), mar = c(4.0, 5.1, 2.7, 1.4))
  
  for (k in 1:3) {
    a <- ratio_curves[ratio_curves$Framework == "XR_marginal" &
                        ratio_curves$Scheme == paste0("XR", k), ]
    b <- ratio_curves[ratio_curves$Framework == "Case4_both_cont" &
                        ratio_curves$Scheme == paste0("ER", k), ]
    # Avoid extrapolating curves in extremely sparse tails.
    a <- a[is.finite(a[[metric]]) & a$EffectiveN >= 10, ]
    b <- b[is.finite(b[[metric]]) & b$EffectiveN >= 10, ]
    if (!nrow(a) || !nrow(b)) {
      plot.new()
      title(main = paste0("Scenario ", k, ": insufficient local data"))
      next
    }
    ylimits <- range(c(a[[metric]], b[[metric]]), finite = TRUE)
    if (diff(ylimits) == 0) ylimits <- ylimits + c(-0.1, 0.1)
    plot(a$DeductibleRatio, a[[metric]],
         type = "l", lwd = 2, col = "steelblue4",
         log = "x", ylim = ylimits,
         xlab = "Deductible / Coverage (log scale)",
         ylab = metric, main = paste0("XR", k, " vs continuous ER", k))
    lines(b$DeductibleRatio, b[[metric]],
          col = "firebrick3", lwd = 2, lty = 2)
    if (metric == "OP_ratio") abline(h = 1, lty = 3, col = "gray50")
    legend("topright",
           legend = c(paste0("XR", k), paste0("ER", k, " continuous")),
           col = c("steelblue4", "firebrick3"),
           lty = c(1, 2), lwd = 2, bty = "n")
  }
  invisible(outfile)
}

invisible(lapply(c("MAE", "RMSE", "MPD", "OP_ratio"), plot_ratio_comparison))
cat("Continuous deductible-ratio outputs saved under: ",
    cfg$output_dir, "\n", sep = "")
cat("Smoothing bandwidth on log10(D/C): ", signif(ratio_bandwidth, 4), "\n")


# ================================================================
# 2010 OOS: Premium performance over deductible ratio
#
# Correct ratio:
#   Deductible (USD) / Coverage (million USD) / 1e6
#
# Full vs ExcludeMax:
#   - Same x-axis range
#   - Same smoothing bandwidth
#   - Same evaluation methodology
#
# No deductible-ratio grouping
# ================================================================

oos_deductible_ratio_diagnostics <- function(
    oos_pol,
    oos_predictions,
    model_pairs,
    oos_valid,
    oos_file,
    exclude_max = FALSE,
    coverage_to_dollars = 1e6,
    ratio_xlim_pct = NULL,
    bandwidth = NULL,
    min_ess = 15) {
  
  idx <- which(oos_valid)
  
  if (!length(idx)) {
    stop("No valid 2010 observations.")
  }
  
  # --------------------------------------------------------------
  # 0. Correct deductible ratio
  # --------------------------------------------------------------
  
  coverage_dollars <-
    oos_pol$Coverage[idx] * coverage_to_dollars
  
  ratio <- oos_pol$Deduct[idx] / coverage_dollars
  
  obs <- oos_pol$S_limited[idx]
  
  stopifnot(
    all(is.finite(ratio) & ratio > 0),
    all(is.finite(obs) & obs >= 0)
  )
  
  tag <- if (exclude_max) "ExcludeMax" else "Full"
  
  # Plot percentage, store dimensionless ratio
  ratio_pct <- 100 * ratio
  
  # Same plotting range for both OOS scenarios
  if (is.null(ratio_xlim_pct)) {
    ratio_xlim_pct <- range(ratio_pct)
  }
  
  stopifnot(
    length(ratio_xlim_pct) == 2L,
    all(is.finite(ratio_xlim_pct)),
    all(ratio_xlim_pct > 0),
    diff(ratio_xlim_pct) > 0
  )
  
  # --------------------------------------------------------------
  # 1. Policy-level prediction errors
  # --------------------------------------------------------------
  
  rows <- lapply(seq_len(nrow(model_pairs)), function(j) {
    
    pp <- model_pairs[j, , drop = FALSE]
    
    pred <- as.numeric(
      oos_predictions[[pp$Model]]$Limited[idx]
    )
    
    stopifnot(
      length(pred) == length(obs),
      all(is.finite(pred) & pred >= 0)
    )
    
    # Poisson deviance contribution
    dev <- rep(Inf, length(pred))
    
    zero <- obs == 0
    pos <- obs > 0 & pred > 0
    
    dev[zero] <- 2 * pred[zero]
    
    dev[pos] <- 2 * (
      obs[pos] * log(obs[pos] / pred[pos]) -
        (obs[pos] - pred[pos])
    )
    
    dev[is.finite(dev)] <-
      pmax(0, dev[is.finite(dev)])
    
    data.frame(
      Framework = pp$Framework,
      Scheme = pp$Scheme,
      Model = pp$Model,
      PolicyNum = oos_pol$PolicyNum[idx],
      Year = oos_pol$Year[idx],
      Deductible = oos_pol$Deduct[idx],
      Coverage = oos_pol$Coverage[idx],
      Coverage_Dollars = coverage_dollars,
      DeductibleRatio = ratio,
      DeductibleRatio_pct = ratio_pct,
      Observed = obs,
      Predicted = pred,
      AbsoluteError = abs(obs - pred),
      SquaredError = (obs - pred)^2,
      PoissonDeviance = dev
    )
  })
  
  policy <- do.call(rbind, rows)
  
  write.csv(
    policy,
    oos_file(paste0(
      "DeductibleRatio_PolicyLevel_2010_",
      tag, ".csv"
    )),
    row.names = FALSE
  )
  
  # --------------------------------------------------------------
  # 2. Continuous kernel smoothing
  # --------------------------------------------------------------
  
  x <- log10(ratio)
  
  if (is.null(bandwidth)) {
    bandwidth <- max(0.08, 0.20 * IQR(x))
  }
  
  stopifnot(
    is.finite(bandwidth),
    bandwidth > 0
  )
  
  # Common grid based on the same full-sample axis limits
  grid <- 10^seq(
    log10(ratio_xlim_pct[1] / 100),
    log10(ratio_xlim_pct[2] / 100),
    length.out = 160L
  )
  
  smooth_one <- function(d) {
    
    x_d <- log10(d$DeductibleRatio)
    
    out <- lapply(grid, function(r) {
      
      # Stable Gaussian kernel weights
      logw <- -0.5 * (
        (x_d - log10(r)) / bandwidth
      )^2
      
      w <- exp(logw - max(logw))
      w <- w / sum(w)
      
      ess <- 1 / sum(w^2)
      
      weighted_obs <- sum(w * d$Observed)
      weighted_pred <- sum(w * d$Predicted)
      
      # Preserve infinite Poisson deviance
      # when observed > 0 and predicted = 0
      mpd <- if (any(
        is.infinite(d$PoissonDeviance) & w > 0
      )) {
        Inf
      } else {
        sum(w * d$PoissonDeviance)
      }
      
      data.frame(
        DeductibleRatio = r,
        DeductibleRatio_pct = 100 * r,
        EffectiveN = ess,
        MAE = sum(w * d$AbsoluteError),
        RMSE = sqrt(sum(w * d$SquaredError)),
        MPD = mpd,
        OP_ratio = if (weighted_pred > 0) {
          weighted_obs / weighted_pred
        } else {
          NA_real_
        }
      )
    })
    
    result <- do.call(rbind, out)
    
    result$Framework <- d$Framework[1]
    result$Scheme <- d$Scheme[1]
    result$Model <- d$Model[1]
    
    result
  }
  
  smooth <- do.call(
    rbind,
    lapply(
      split(policy, policy$Model),
      smooth_one
    )
  )
  
  rownames(smooth) <- NULL
  
  write.csv(
    smooth,
    oos_file(paste0(
      "DeductibleRatio_Smoothed_2010_",
      tag, ".csv"
    )),
    row.names = FALSE
  )
  
  # --------------------------------------------------------------
  # 3. Plot XR versus continuous ER
  # --------------------------------------------------------------
  
  plot_metric <- function(metric) {
    
    filename <- oos_file(paste0(
      "DeductibleRatio_", metric,
      "_XR_vs_ContinuousER_2010_",
      tag, ".pdf"
    ))
    
    grDevices::pdf(
      filename,
      width = 9,
      height = 10
    )
    
    old_par <- par(no.readonly = TRUE)
    
    on.exit({
      par(old_par)
      grDevices::dev.off()
    }, add = TRUE)
    
    par(
      mfrow = c(3, 1),
      mar = c(4.5, 5, 3, 1)
    )
    
    for (k in 1:3) {
      
      xr <- smooth[
        smooth$Framework == "XR_marginal" &
          smooth$Scheme == paste0("XR", k),
        ,
        drop = FALSE
      ]
      
      er <- smooth[
        smooth$Framework == "Case4_both_cont" &
          smooth$Scheme == paste0("ER", k),
        ,
        drop = FALSE
      ]
      
      # Suppress low-information tails
      xr <- xr[
        is.finite(xr[[metric]]) &
          xr$EffectiveN >= min_ess,
        ,
        drop = FALSE
      ]
      
      er <- er[
        is.finite(er[[metric]]) &
          er$EffectiveN >= min_ess,
        ,
        drop = FALSE
      ]
      
      if (!nrow(xr) || !nrow(er)) {
        
        plot.new()
        
        title(
          main = paste(
            "Scenario", k,
            ": insufficient data"
          )
        )
        
        next
      }
      
      ylim <- range(
        c(xr[[metric]], er[[metric]]),
        finite = TRUE
      )
      
      if (diff(ylim) == 0) {
        delta <- max(0.1, abs(ylim[1]) * 0.1)
        ylim <- ylim + c(-delta, delta)
      }
      
      # Deductible ratio displayed as percentage
      plot(
        xr$DeductibleRatio_pct,
        xr[[metric]],
        type = "l",
        col = "steelblue4",
        lwd = 2,
        log = "x",
        xlim = ratio_xlim_pct,
        ylim = ylim,
        xlab = "Deductible / Coverage (%)",
        ylab = metric,
        main = paste0(
          "2010 ", tag,
          ": XR", k,
          " vs ER", k
        )
      )
      
      lines(
        er$DeductibleRatio_pct,
        er[[metric]],
        col = "firebrick3",
        lwd = 2,
        lty = 2
      )
      
      if (metric == "OP_ratio") {
        
        abline(
          h = 1,
          lty = 3,
          col = "gray50"
        )
      }
      
      legend(
        "topright",
        legend = c(
          paste0("XR", k),
          paste0("Continuous ER", k)
        ),
        col = c(
          "steelblue4",
          "firebrick3"
        ),
        lty = c(1, 2),
        lwd = 2,
        bty = "n"
      )
    }
    
    invisible(filename)
  }
  
  invisible(lapply(
    c("MAE", "RMSE", "MPD", "OP_ratio"),
    plot_metric
  ))
  
  cat(
    "\n2010 deductible-ratio diagnostics completed:",
    tag,
    "\nCoverage multiplier:", coverage_to_dollars,
    "\nRatio range (%):",
    paste(signif(ratio_xlim_pct, 4), collapse = " to "),
    "\nBandwidth:", signif(bandwidth, 4),
    "\n"
  )
  
  invisible(smooth)
}



################################################################################
# OUT-OF-SAMPLE VALIDATION: POLICY YEAR 2010
# MERGED INTO THE EMPIRICAL 2006-2009 ESTIMATION SCRIPT.
# Alternatively, source this file after restoring the completed training workspace.
#
# NO refitting, NO copula reselection, NO Godambe recomputation.
# 2010 outcomes are accessed ONLY for observed-loss evaluation and scoring.
# Nine headline model pairs are evaluated: Cases 1, 4c and 4d.
# The primary loss: sum_j (min(Y_ij, C_i) - D_i)_+, C_i = Coverage_i * scale.
################################################################################

run_oos_2010 <- function(exclude_max = FALSE, reuse_full = NULL) {
cat("\n================ OUT-OF-SAMPLE: YEAR 2010 | ",
    if (exclude_max) "EXCLUDE MAX" else "FULL SAMPLE",
    " ================\n", sep = "")

oos_year <- 2010L
needed <- c("cfg", "poldat", "claimdat", "pol", "clm", "X_pol", "fml",
            "model_pairs", "frequency_fits", "severity_fits", "fit_Dc", "fit_Dd",
            "prepare_eval", "conditional_N_mean", "conditional_paid_limited",
            "conditional_paid_unlimited", "gl_nodes", "raw_predictions")
missing_items <- needed[!vapply(needed, exists, logical(1), inherits = TRUE)]
if (length(missing_items)) stop("Run the existing empirical script first. Missing: ",
                                paste(missing_items, collapse = ", "))
if (!identical(as.integer(cfg$years_in), 2006:2009)) {
  stop("Expected fitted training years cfg$years_in = 2006:2009; current setting is: ",
       paste(cfg$years_in, collapse = ", "))
}
if (any(pol$Year == oos_year) || any(clm$Year == oos_year))
  stop("2010 data appear in training objects. Out-of-sample assessment is invalid.")
if (length(cfg$holdout_years))
  stop("This appended analysis expects the original 2006-2009 in-sample fit, with holdout_years = integer(0).")

# Keep original cfg, eval_pol, eval_dat, and in-sample predictions unchanged.
oos_dir <- file.path(cfg$output_dir,
                     if (exclude_max) "out_of_sample_2010_exclude_max_loss"
                     else "out_of_sample_2010")
oos_report_dir <- file.path(oos_dir, "manuscript_reporting")
dir.create(oos_report_dir, recursive = TRUE, showWarnings = FALSE)
oos_file <- function(name) file.path(oos_dir, name)
oos_report_file <- function(name) file.path(oos_report_dir, name)

oos2010_policies <- poldat[poldat$Year == oos_year, , drop = FALSE]
oos2010_claims <- claimdat[claimdat$Year == oos_year, , drop = FALSE]
if (!exclude_max) {
  # Descriptive 2010-only statistics, unmodified by the sensitivity exclusion.
  pos_2010 <- pmax(0, oos2010_claims$LossBeforeDeductible - oos2010_claims$Deduct)
  nn2010 <- as.integer(oos2010_claims$LossBeforeDeductible > oos2010_claims$Deduct)
  count_2010 <- if (nrow(oos2010_claims))
    stats::aggregate(nn2010,
      by = list(PolicyNum = oos2010_claims$PolicyNum, Year = oos2010_claims$Year), FUN = sum)
    else data.frame(PolicyNum = character(0), Year = integer(0), x = integer(0))
  keys_pol <- paste(oos2010_policies$PolicyNum, oos2010_policies$Year)
  keys_claim <- paste(count_2010$PolicyNum, count_2010$Year)
  nstar <- count_2010$x[match(keys_pol, keys_claim)]
  nstar[is.na(nstar)] <- 0
  stat2010 <- function(x) c(Mean = if (length(x)) mean(x) else NA_real_,
                            SD = if (length(x) > 1L) sd(x) else NA_real_)
  summary_2010 <- t(rbind(
    N = stat2010(oos2010_policies$Freq), Nstar = stat2010(nstar),
    Y = stat2010(oos2010_claims$LossBeforeDeductible),
    Y_censored = stat2010(pos_2010),
    Y_truncated = stat2010(pos_2010[pos_2010 > 0])))
  dir.create(file.path(cfg$output_dir, "out_of_sample_2010"), showWarnings = FALSE)
  write.csv(summary_2010, file.path(cfg$output_dir, "out_of_sample_2010", "Summary_Statistics_2010.csv"))
}
if (!nrow(oos2010_policies)) stop("There are no 2010 policies in poldat.")
if (anyDuplicated(paste(oos2010_policies$PolicyNum, oos2010_policies$Year)))
  stop("Duplicate PolicyNum-Year keys in 2010 policy records.")
if (any(!is.finite(oos2010_policies$Coverage) | oos2010_policies$Coverage <= 0))
  stop("Invalid 2010 policy coverage; cannot construct design/premium inputs.")
if (any(!is.finite(oos2010_policies$Deduct) | oos2010_policies$Deduct <= 0))
  stop("Invalid 2010 policy deductible.")
if (anyNA(match(oos2010_policies$Deduct, cfg$ded_levels)))
  stop("2010 contains deductible levels outside cfg$ded_levels; do not silently drop these policies.")
if (nrow(oos2010_claims)) {
  bad_claims <- !paste(oos2010_claims$PolicyNum, oos2010_claims$Year) %in%
    paste(oos2010_policies$PolicyNum, oos2010_policies$Year)
  if (any(bad_claims)) stop(sum(bad_claims),
    " 2010 claims have no matching 2010 policy record. Resolve linkage before scoring.")
  if (any(!is.finite(oos2010_claims$LossBeforeDeductible) |
          oos2010_claims$LossBeforeDeductible < 0))
    stop("Nonfinite or negative 2010 ground-up claim amounts.")
} else {
  warning("There are no 2010 claim records: all observed 2010 aggregate losses will be zero.")
}

# Exact same preprocessing/payment definition used in original evaluation.
oos_pol <- prepare_eval(oos2010_policies, oos2010_claims)
# POST-HOC sensitivity exclusion, same fixed 2006-2009 parameter estimates.
# Exclude all claims for the policy-year with largest limited aggregate S.
full_oos_pol <- oos_pol
oos_sel <- seq_len(nrow(oos_pol))
excluded_policy <- NULL
if (exclude_max) {
  if (nrow(oos_pol) < 2L) stop("Need at least two 2010 policy-years")
  if (any(!is.finite(oos_pol$S_limited)))
    stop("Invalid 2010 observed aggregate losses")
  imax <- which(oos_pol$S_limited == max(oos_pol$S_limited))
  if (length(imax) != 1L) stop("Maximum 2010 loss is tied; specify policy manually")
  excluded_policy <- oos_pol[imax, , drop = FALSE]
  oos_sel <- seq_len(nrow(oos_pol))[-imax]
  excluded_policy$Exclusion_reason <-
    "Largest 2010 observed aggregate coverage-limited payment (post hoc)"
  write.csv(excluded_policy,
            oos_file("Excluded_2010_Extreme_PolicyYear.csv"), row.names = FALSE)
  oos_pol <- oos_pol[-imax, , drop = FALSE]
  rownames(oos_pol) <- NULL
  cat("  Excluded one 2010 policy-year: ",
      as.character(excluded_policy$PolicyNum),
      "; observed limited S = ", excluded_policy$S_limited, "\n", sep = "")
}
oos_X <- model.matrix(fml, data = oos_pol)
if (!identical(colnames(oos_X), colnames(X_pol)))
  stop("2010 design-matrix columns differ from training design-matrix columns.")

oos_dat <- list(
  X = oos_X,
  d = oos_pol$Deduct,
  coverage_limit = oos_pol$CoverageLimit_Dollars,
  dratio = oos_pol$dedratio,
  dlev = oos_pol$discDed,
  n_cen = oos_pol$n_ground,
  n_tru = oos_pol$n_paid,
  unit = seq_len(nrow(oos_pol))
)
oos_zero_limit <- oos_pol$CoverageBelowOrEqualDeductible
stopifnot(all(oos_pol$S_limited[oos_zero_limit] == 0),
          all(oos_pol$S_limited >= 0),
          all(oos_pol$S_limited <= oos_pol$S_raw + 1e-7))

# In this step the fitted theta_D, theta_Y, theta_N, and copula families
# are fixed at their 2006-2009 estimates. Outcomes are not input to predictors.
oos_Dtheta <- function(fit) switch(fit$spec$dedtype,
                                   cont = fit_Dc$par, disc = fit_Dd$par, none = NULL)
if (is.null(reuse_full)) {
  cat("Predicting 2010 conditional frequency means (fixed 2006-2009 fits)...\n")
  oos_freq_mean <- lapply(frequency_fits[premium_frequency_models], function(f)
    as.numeric(conditional_N_mean(f, oos_Dtheta(f), oos_dat)))
  cat("Predicting 2010 coverage-limited severity payments (fixed fits)...\n")
  oos_sev_limited <- lapply(severity_fits, function(f)
    as.numeric(conditional_paid_limited(f, oos_Dtheta(f), oos_dat)))
} else {
  # Predictions come from FULL 2010 inputs; only evaluation rows are removed.
  # No re-estimation AND no repeated prediction integration for the excluded case.
  stopifnot(exclude_max,
            length(oos_sel) == nrow(oos_pol),
            identical(as.character(full_oos_pol$PolicyNum),
                      as.character(reuse_full$full_policy_data$PolicyNum)))
  oos_freq_mean <- lapply(reuse_full$freq_mean, function(z) z[oos_sel])
  oos_sev_limited <- lapply(reuse_full$sev_limited, function(z) z[oos_sel])
}

if (!all(vapply(oos_sev_limited, function(v)
    length(v) == nrow(oos_pol) && all(v[oos_zero_limit] == 0), logical(1))))
  stop("2010 limited severity payments violate the contractual C <= D zero-payment condition.")

# Supplementary: retain the same unlimited calculation and GB2 moment guard.
if (!is.null(reuse_full)) {
  oos_sev_unlimited <- lapply(reuse_full$sev_unlimited, function(z) z[oos_sel])
} else if (isTRUE(cfg$compute_unlimited)) {
  cat("Predicting supplementary unlimited severity payments for 2010...\n")
  oos_sev_unlimited <- lapply(severity_fits, function(f)
    as.numeric(conditional_paid_unlimited(f, oos_Dtheta(f), oos_dat)))
} else {
  oos_sev_unlimited <- lapply(severity_fits, function(f)
    rep(NA_real_, nrow(oos_pol)))
}

# Preserve EXACTLY the nine Case 1 / 4c / 4d model pairings.
oos_main_pairs <- model_pairs  # exactly nine Case 1, 4c, 4d schemes
stopifnot(nrow(oos_main_pairs) == 9L,
          all(oos_main_pairs$Frequency %in% names(oos_freq_mean)),
          all(oos_main_pairs$Severity %in% names(oos_sev_limited)))

oos_predictions <- setNames(vector("list", nrow(model_pairs)), model_pairs$Model)
for (ii in seq_len(nrow(model_pairs))) {
  pp <- model_pairs[ii, ]
  pred_limited <- as.numeric(oos_freq_mean[[pp$Frequency]] *
                              oos_sev_limited[[pp$Severity]])
  pred_unlimited <- as.numeric(oos_freq_mean[[pp$Frequency]] *
                                oos_sev_unlimited[[pp$Severity]])
  pred_limited[oos_zero_limit] <- 0
  if (length(pred_limited) != nrow(oos_pol)) stop("Incorrect prediction length: ", pp$Model)
  oos_predictions[[pp$Model]] <- list(Limited = pred_limited,
                                       Unlimited = pred_unlimited)
}
saveRDS(oos_predictions, oos_file("Premium_Predictions_ByPolicy_2010.rds"))

# Export individual policy-year predictions for reproducible paired comparisons.
oos_prem_wide <- data.frame(
  PolicyNum = oos_pol$PolicyNum, Year = oos_pol$Year,
  Deductible = oos_pol$Deduct, Coverage = oos_pol$Coverage,
  CoverageLimit_Dollars = oos_pol$CoverageLimit_Dollars,
  CoverageBelowOrEqualDeductible = oos_zero_limit,
  Observed_Aggregate_Limited = oos_pol$S_limited,
  Observed_Aggregate_Unlimited = oos_pol$S_raw,
  stringsAsFactors = FALSE
)
for (nm in names(oos_predictions))
  oos_prem_wide[[paste0("Premium_Limited__", nm)]] <- oos_predictions[[nm]]$Limited
write.csv(oos_prem_wide, oos_file("Premium_Predictions_ByPolicy_2010.csv"), row.names = FALSE)
oos_prem_long <- do.call(rbind, lapply(seq_len(nrow(model_pairs)), function(i) {
  pp <- model_pairs[i, , drop = FALSE]
  data.frame(Case = pp$Case, Framework = pp$Framework,
             Scheme = pp$Scheme, Model = pp$Model,
             PolicyNum = oos_pol$PolicyNum, Year = oos_pol$Year,
             Deductible = oos_pol$Deduct,
             Observed_Aggregate_Limited = oos_pol$S_limited,
             Premium_Limited = oos_predictions[[pp$Model]]$Limited,
             stringsAsFactors = FALSE)
}))
write.csv(oos_prem_long,
          oos_file("Premium_Predictions_ByPolicy_Long_2010.csv"), row.names = FALSE)

# -------------------- IDENTICAL METRICS / ZERO CONVENTIONS --------------------
oos_op <- function(y, p) {
  sy <- sum(y); sp <- sum(p)
  if (!is.finite(sy) || !is.finite(sp)) return(NA_real_)
  if (sp > 0) return(sy / sp)
  if (sy > 0) return(Inf)
  NA_real_  # O/P undefined when both observed and predicted totals are zero.
}

oos_poisson_deviance <- function(y, p) {
  stopifnot(length(y) == length(p),
            all(is.finite(y)), all(is.finite(p)), all(y >= 0), all(p >= 0))
  if (!length(y)) return(NA_real_)
  a <- numeric(length(y))
  a[y == 0] <- 2 * p[y == 0]
  a[y > 0 & p == 0] <- Inf
  ok <- y > 0 & p > 0
  a[ok] <- 2 * (y[ok] * log(y[ok] / p[ok]) - (y[ok] - p[ok]))
  a[is.finite(a)] <- pmax(a[is.finite(a)], 0)
  mean(a)
}
stopifnot(isTRUE(all.equal(oos_poisson_deviance(c(0, 12), c(0, 12)), 0)),
          isTRUE(all.equal(oos_poisson_deviance(0, 3), 6)),
          is.infinite(oos_poisson_deviance(2, 0)))

oos_scores <- function(pred, obs) {
  stopifnot(length(pred) == length(obs))
  good <- is.finite(pred) & is.finite(obs) & pred >= 0 & obs >= 0
  p <- as.numeric(pred[good]); y <- as.numeric(obs[good]); n <- length(y)
  if (!n) return(data.frame(
    N = 0L, N_positive = 0L,
    Observed_total = NA_real_, Predicted_total = NA_real_,
    Observed_mean = NA_real_, Predicted_mean = NA_real_,
    OP_ratio = NA_real_, RelBias_pct = NA_real_,
    Mean_Poisson_Deviance = NA_real_, MAE_obs = NA_real_,
    RMSE_obs = NA_real_, RMSRE_positive = NA_real_))
  pos <- y > 0
  data.frame(
    N = n, N_positive = sum(pos),
    Observed_total = sum(y), Predicted_total = sum(p),
    Observed_mean = mean(y), Predicted_mean = mean(p),
    OP_ratio = oos_op(y, p),
    RelBias_pct = if (sum(y) > 0) 100 * (sum(p) - sum(y)) / sum(y) else NA_real_,
    Mean_Poisson_Deviance = oos_poisson_deviance(y, p),
    MAE_obs = mean(abs(p - y)),
    RMSE_obs = sqrt(mean((p - y)^2)),
    RMSRE_positive = if (any(pos))
      sqrt(mean(((p[pos] - y[pos]) / y[pos])^2)) else NA_real_
  )
}

# Main comparison: EXACTLY THE SAME 2010 POLICY-YEARS for all nine models.
oos_valid <- is.finite(oos_pol$S_limited) & oos_pol$S_limited >= 0
# ================================================================
# Deductible-ratio diagnostics: 2010 OOS
# ================================================================

# Common x-axis range based on the FULL 2010 dataset.
# full_oos_pol is already defined earlier in run_oos_2010().

full_ratio <- full_oos_pol$Deduct /
  (full_oos_pol$Coverage * cfg$coverage_to_dollars)

ratio_xlim_pct <- range(
  100 * full_ratio,
  finite = TRUE
)

# Use one bandwidth for both sensitivity scenarios
full_x <- log10(full_ratio)

ratio_bandwidth <- max(
  0.08,
  0.20 * IQR(full_x)
)

oos_ratio_results <- oos_deductible_ratio_diagnostics(
  oos_pol = oos_pol,
  oos_predictions = oos_predictions,
  model_pairs = oos_main_pairs,
  oos_valid = oos_valid,
  oos_file = oos_file,
  exclude_max = exclude_max,
  coverage_to_dollars = cfg$coverage_to_dollars,
  ratio_xlim_pct = ratio_xlim_pct,
  bandwidth = ratio_bandwidth,
  min_ess = 15
)

oos_completeness <- data.frame(Model = oos_main_pairs$Model,
                               Invalid_predictor_rows = NA_integer_)
for (ii in seq_len(nrow(oos_main_pairs))) {
  nm <- oos_main_pairs$Model[ii]
  p <- oos_predictions[[nm]]$Limited
  good <- is.finite(p) & p >= 0
  oos_completeness$Invalid_predictor_rows[ii] <- sum(!good)
  oos_valid <- oos_valid & good
}
write.csv(oos_completeness,
          oos_report_file("Appendix_Prediction_Completeness_2010.csv"), row.names = FALSE)
n_oos <- sum(oos_valid)
oos_ratio_results <- oos_deductible_ratio_diagnostics(
  oos_pol = oos_pol,
  oos_predictions = oos_predictions,
  model_pairs = oos_main_pairs,
  oos_valid = oos_valid,
  oos_file = oos_file,
  exclude_max = exclude_max
)
if (!n_oos) stop("No common valid 2010 policy-year observations across all nine headline models. Check completeness CSV.")
if (n_oos < nrow(oos_pol)) warning("Headline comparison retains ", n_oos,
                                   " / ", nrow(oos_pol), " valid 2010 policy-years.")

# Table 1 and Table 2 are NOT reselected on the 2010 validation sample.
# Copy the chosen TRAINING specifications solely as reference tables.
if (exists("selected_severity")) write.csv(selected_severity,
    oos_report_file("Table1_TrainingSelectedSeverityCopulas_Reference.csv"), row.names = FALSE)
if (exists("selected_frequency")) write.csv(selected_frequency,
    oos_report_file("Table2_TrainingSelectedFrequencyCopulas_Reference.csv"), row.names = FALSE)

oos_all_scores <- do.call(rbind, lapply(seq_len(nrow(oos_main_pairs)), function(i) {
  pp <- oos_main_pairs[i, , drop = FALSE]
  cbind(pp[, c("Case", "Framework", "Scheme", "Model", "Frequency", "Severity"), drop = FALSE],
        oos_scores(oos_predictions[[pp$Model]]$Limited[oos_valid],
                   oos_pol$S_limited[oos_valid]))
}))
rownames(oos_all_scores) <- NULL
oos_table3 <- oos_all_scores[, c(
  "Case", "Framework", "Scheme", "Model", "Frequency", "Severity", "N",
  "OP_ratio", "RelBias_pct", "Mean_Poisson_Deviance",
  "MAE_obs", "RMSE_obs", "Observed_total", "Predicted_total"
)]
write.csv(oos_table3,
          oos_report_file("Table3_2010_Limited_PurePremium_Main.csv"), row.names = FALSE)
write.csv(oos_all_scores,
          oos_report_file("Appendix_2010_Limited_PurePremium_AllMetrics.csv"), row.names = FALSE)

# No Case 3 / thinning-sensitivity ratemaking combinations are reported.

# Additional output matching original raw (limited + unlimited) metrics.
oos_metric_parts <- list(); oos_metric_k <- 0L
oos_d_metric_parts <- list(); oos_d_metric_k <- 0L
for (i in seq_len(nrow(model_pairs))) {
  pp <- model_pairs[i, , drop = FALSE]
  for (kind in c("Limited", "Unlimited")) {
    pred <- oos_predictions[[pp$Model]][[kind]]
    obs <- if (kind == "Limited") oos_pol$S_limited else oos_pol$S_raw
    oos_metric_k <- oos_metric_k + 1L
    oos_metric_parts[[oos_metric_k]] <- cbind(
      pp[, c("Case", "Framework", "Scheme", "Model", "Frequency", "Severity"), drop = FALSE],
      Payment_definition = kind, N_all = nrow(oos_pol),
      oos_scores(pred, obs), Sample = if (exclude_max) "Out_of_sample_2010_exclude_max" else "Out_of_sample_2010"
    )
    for (d in cfg$ded_levels) {
      ix <- which(oos_pol$Deduct == d)
      if (!length(ix)) next
      oos_d_metric_k <- oos_d_metric_k + 1L
      oos_d_metric_parts[[oos_d_metric_k]] <- cbind(
        pp[, c("Case", "Framework", "Scheme", "Model", "Frequency", "Severity"), drop = FALSE],
        Payment_definition = kind, Deductible = d, N_all = length(ix),
        oos_scores(pred[ix], obs[ix]), Sample = if (exclude_max) "Out_of_sample_2010_exclude_max" else "Out_of_sample_2010"
      )
    }
  }
}
oos_overall_metrics <- do.call(rbind, oos_metric_parts)
oos_by_d_metrics <- do.call(rbind, oos_d_metric_parts)
write.csv(oos_overall_metrics, oos_file("Premium_Performance_Overall_2010.csv"), row.names = FALSE)
write.csv(oos_by_d_metrics, oos_file("Premium_Performance_ByDeductible_2010.csv"), row.names = FALSE)
write.csv(oos_overall_metrics[oos_overall_metrics$Payment_definition == "Unlimited", , drop = FALSE],
          oos_report_file("Appendix_2010_Unlimited_Premiums.csv"), row.names = FALSE)

# Keep the original primary/supplementary export naming convention.
write.csv(oos_overall_metrics[oos_overall_metrics$Payment_definition == "Limited", , drop = FALSE],
          oos_file("Premium_Performance_CoverageLimited_Primary.csv"), row.names = FALSE)
write.csv(oos_overall_metrics[oos_overall_metrics$Payment_definition == "Unlimited", , drop = FALSE],
          oos_file("Premium_Performance_Unlimited_Supplementary.csv"), row.names = FALSE)

# Reproduce diagnostics for the THREE frameworks (Cases 1/4c/4d).
oos_plot_deductible_metric <- function(metric, filename) {
  tt <- oos_by_d_metrics[
    oos_by_d_metrics$Payment_definition == "Limited", , drop = FALSE]
  z <- tt[[metric]]
  z <- z[is.finite(z)]
  if (!length(z)) return(invisible(NULL))
  use_ylog <- metric != "RelBias_pct" && all(z > 0)
  ylim <- if (metric == "RelBias_pct") {
    b <- max(1, abs(z)); c(-1.05*b, 1.05*b)
  } else if (use_ylog) {
    r <- range(z); c(r[1]/1.2, r[2]*1.2)
  } else {
    r <- range(z); r + c(-.08, .08) * max(diff(r), 1)
  }
  col <- c("firebrick3", "steelblue4", "darkgreen")
  grDevices::pdf(oos_file(filename), width = 11, height = 8)
  op <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2), mar = c(4.2, 4.4, 2.6, .8))
  for (g in unique(oos_main_pairs$Framework)) {
    plot(NA_real_, NA_real_, xlim = range(cfg$ded_levels), ylim = ylim,
         log = if (use_ylog) "xy" else "x", xaxt = "n",
         xlab = "Deductible (USD)",
         ylab = if (metric == "RMSRE_positive") "RMSRE (observed S > 0)" else metric,
         main = g)
    axis(1, at = cfg$ded_levels,
         labels = format(cfg$ded_levels, scientific = FALSE, big.mark = ","),
         las = 2, cex.axis = .6)
    if (metric == "RelBias_pct") abline(h = 0, lty = 3, col = "grey50")
    for (j in 1:3) {
      sc <- paste0(if (g == "XR_marginal") "XR" else "ER", j)
      aa <- tt[tt$Framework == g & tt$Scheme == sc, , drop = FALSE]
      aa <- aa[order(aa$Deductible), , drop = FALSE]
      if (nrow(aa)) lines(aa$Deductible, aa[[metric]], type = "b",
                           col = col[j], lty = j, pch = 15+j)
    }
    legend("topleft", legend = paste0(if (g == "XR_marginal") "XR" else "ER", 1:3),
           col = col, lty = 1:3, pch = 16:18, bty = "n", cex = .72)
  }
  plot.new()
  text(.5, .65, "2010 OUT-OF-SAMPLE", cex = 1)
  text(.5, .4, "Training 2006-2009; fixed model parameters", cex = .92)
  par(op); dev.off()
  invisible(NULL)
}
oos_plot_deductible_metric("MAE_obs", "Premium_ByDeductible_MAE_obs_Limited.pdf")
oos_plot_deductible_metric("RMSRE_positive", "Premium_ByDeductible_RMSRE_obs_Limited.pdf")
oos_plot_deductible_metric("RelBias_pct", "Premium_ByDeductible_RB_obs_pct_Limited.pdf")

# -------------------- DEDUCTIBLE-LEVEL O/P AND 95% BOOTSTRAP CI -------------
oos_ded_parts <- list(); kk <- 0L
for (i in seq_len(nrow(oos_main_pairs))) {
  pp <- oos_main_pairs[i, , drop = FALSE]
  p <- oos_predictions[[pp$Model]]$Limited
  for (d in cfg$ded_levels) {
    ix <- which(oos_valid & oos_pol$Deduct == d)
    if (!length(ix)) next
    kk <- kk + 1L
    oos_ded_parts[[kk]] <- cbind(
      pp[, c("Case", "Framework", "Scheme", "Model"), drop = FALSE],
      Deductible = d, oos_scores(p[ix], oos_pol$S_limited[ix])
    )
  }
}
if (!length(oos_ded_parts)) stop("No scored deductible groups in 2010.")
oos_ded_scores <- do.call(rbind, oos_ded_parts)
rownames(oos_ded_scores) <- NULL

oos_cluster_ci <- function(y, p, ids, B, min_clusters) {
  sums <- rowsum(cbind(y = y, p = p), as.character(ids), reorder = FALSE)
  K <- nrow(sums)
  if (K < min_clusters || B < 2L)
    return(c(N_clusters = K, OP_lower = NA_real_, OP_upper = NA_real_))
  v <- rep(NA_real_, B)
  for (bb in seq_len(B)) {
    ix <- sample.int(K, K, replace = TRUE)
    v[bb] <- oos_op(sum(sums[ix, "y"]), sum(sums[ix, "p"]))
  }
  v <- v[is.finite(v)]
  if (length(v) < max(20L, ceiling(B / 2)))
    return(c(N_clusters = K, OP_lower = NA_real_, OP_upper = NA_real_))
  limits <- as.numeric(stats::quantile(v, c(.025, .975), names = FALSE))
  c(N_clusters = K, OP_lower = limits[1L], OP_upper = limits[2L])
}

set.seed(cfg$report_bootstrap_seed + 2010L)
oos_ded_scores$N_clusters <- NA_integer_
oos_ded_scores$OP_lower <- NA_real_
oos_ded_scores$OP_upper <- NA_real_
for (i in seq_len(nrow(oos_ded_scores))) {
  d <- oos_ded_scores$Deductible[i]
  nm <- oos_ded_scores$Model[i]
  ix <- which(oos_valid & oos_pol$Deduct == d)
  ci <- oos_cluster_ci(oos_pol$S_limited[ix],
                       oos_predictions[[nm]]$Limited[ix],
                       oos_pol$PolicyNum[ix],
                       B = as.integer(cfg$report_bootstrap_reps),
                       min_clusters = as.integer(cfg$report_min_bootstrap_clusters))
  oos_ded_scores$N_clusters[i] <- as.integer(ci["N_clusters"])
  oos_ded_scores$OP_lower[i] <- ci["OP_lower"]
  oos_ded_scores$OP_upper[i] <- ci["OP_upper"]
}
write.csv(oos_ded_scores,
          oos_report_file("Figure1_Deductible_OP_and_95pctClusterCI.csv"), row.names = FALSE)
write.csv(oos_ded_scores,
          oos_report_file("Figure1B_Deductible_OP_by_ObservationScheme.csv"), row.names = FALSE)

# Figure 1. Each framework gets one panel with the three schemes and CIs.
oos_groups <- unique(oos_main_pairs$Framework)
oos_schemes <- c(1L, 2L, 3L)
oos_scheme_colors <- c("#21618C", "#D35400", "#148F77")
oos_cap <- cfg$report_op_plot_cap

oos_plot_framework <- function() {
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(3, 2), mar = c(4.1, 4.2, 2.5, .8),
      oma = c(.5, 0, 1.4, 0))
  for (g in oos_groups) {
    tt <- oos_ded_scores[oos_ded_scores$Framework == g, , drop = FALSE]
    plot(NA_real_, NA_real_, xlim = range(cfg$ded_levels),
         ylim = c(0, oos_cap), log = "x", xaxt = "n",
         xlab = "Deductible (USD)", ylab = "Observed / predicted",
         main = g)
    axis(1, at = cfg$ded_levels,
         labels = format(cfg$ded_levels, scientific = FALSE, big.mark = ","),
         las = 2, cex.axis = .53)
    abline(h = 1, col = "grey45", lty = 2)
    ss <- unique(tt$Scheme)
    for (j in seq_along(ss)) {
      aa <- tt[tt$Scheme == ss[j], , drop = FALSE]
      aa <- aa[order(aa$Deductible), , drop = FALSE]
      xx <- aa$Deductible * (1 + (j - 2) * .025)
      yy <- pmin(aa$OP_ratio, oos_cap)
      yy[!is.finite(yy)] <- NA_real_
      lines(xx, yy, col = oos_scheme_colors[j], type = "b", lwd = 1.5,
            pch = 15 + j)
      ci <- is.finite(aa$OP_lower) & is.finite(aa$OP_upper)
      if (any(ci)) arrows(xx[ci], pmax(0, pmin(aa$OP_lower[ci], oos_cap)),
                          xx[ci], pmin(aa$OP_upper[ci], oos_cap),
                          angle = 90, code = 3, length = .025,
                          col = oos_scheme_colors[j])
    }
    legend("topright", legend = ss, col = oos_scheme_colors[seq_along(ss)],
           pch = 15 + seq_along(ss), lty = 1, bty = "n", cex = .75)
  }
  plot.new()
  text(.5, .65, "2010 OUT-OF-SAMPLE\n95% descriptive policy-cluster intervals", cex = 1)
  text(.5, .35, sprintf("Display capped at O/P = %s; CSV is uncapped", oos_cap), cex = .88)
  mtext("2010 validation | fitted on 2006-2009", outer = TRUE, side = 3, line = .1)
}
grDevices::pdf(oos_report_file("Figure1_Deductible_OP_95pctCI.pdf"), width = 11, height = 8.5)
oos_plot_framework()
grDevices::dev.off()

# Figure 1B. Each observation scheme gets one plot of THREE frameworks.
oos_framework_cols <- c("#B22222", "#228B22", "#4169E1")[seq_along(oos_groups)]
stopifnot(length(oos_groups) == length(oos_framework_cols))
oos_scheme_num <- function(z) as.integer(sub(".*([123])$", "\\1", z))
oos_ded_scores$Scheme_number <- oos_scheme_num(oos_ded_scores$Scheme)
if (anyNA(oos_ded_scores$Scheme_number) ||
    any(!oos_ded_scores$Scheme_number %in% 1:3))
  stop("Scheme labels must end in 1, 2 or 3.")
oos_scheme_title <- c("Complete", "Censored", "Truncated")
oos_draw_scheme <- function(s, show_legend = TRUE) {
  tt <- oos_ded_scores[oos_ded_scores$Scheme_number == s, , drop = FALSE]
  plot(NA_real_, NA_real_, xlim = range(cfg$ded_levels),
       ylim = c(0, oos_cap), log = "x", xaxt = "n",
       xlab = "Deductible (USD)", ylab = "Observed / predicted",
       main = paste0("Scheme ", s, " - ", oos_scheme_title[s]))
  axis(1, at = cfg$ded_levels,
       labels = format(cfg$ded_levels, scientific = FALSE, big.mark = ","),
       las = 2, cex.axis = .7)
  abline(h = 1, lty = 2, col = "grey45")
  for (j in seq_along(oos_groups)) {
    aa <- tt[tt$Framework == oos_groups[j], , drop = FALSE]
    aa <- aa[order(aa$Deductible), , drop = FALSE]
    yy <- pmin(aa$OP_ratio, oos_cap); yy[!is.finite(yy)] <- NA_real_
    lines(aa$Deductible, yy, type = "b", pch = 14 + j,
          col = oos_framework_cols[j], lty = j, lwd = 1.6)
  }
  if (show_legend) legend("topright", legend = oos_groups,
                          col = oos_framework_cols, lty = seq_along(oos_groups),
                          pch = 14 + seq_along(oos_groups), bty = "n", cex = .68)
}
for (s in 1:3) {
  grDevices::pdf(oos_report_file(sprintf("Figure1B_Scheme%d_%s_Deductible_OP.pdf",
                                       s, oos_scheme_title[s])), width = 10, height = 6.7)
  op <- par(no.readonly = TRUE)
  par(mar = c(6.4, 4.4, 3.2, 1.0))
  oos_draw_scheme(s)
  par(op); dev.off()
}
grDevices::pdf(oos_report_file("Figure1B_AllSchemes_Deductible_OP.pdf"),
               width = 17, height = 6.5)
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3), mar = c(6, 4.3, 3.1, .8), oma = c(.4, .2, 1.2, .2))
for (s in 1:3) oos_draw_scheme(s, show_legend = TRUE)
mtext("2010 OUT-OF-SAMPLE | Cases 1, 4c, 4d",
      outer = TRUE, side = 3, line = .1)
par(op); dev.off()

# -------------------- FIGURE 2: RISK-GROUP CALIBRATION -----------------------
# Define cutpoints using ONLY the 2006-2009 training XR1 predictions.
# Then apply those FIXED cutpoints to 2010; do not refit to 2010 outcomes.
oos_base_name <- "XR_marginal__XR1"
train_xr1 <- raw_predictions[[oos_base_name]]$Limited
holdout_xr1 <- oos_predictions[[oos_base_name]]$Limited
if (is.null(train_xr1) || is.null(holdout_xr1))
  stop("XR1 premiums unavailable for training-defined risk groups.")
train_xr1 <- train_xr1[is.finite(train_xr1) & train_xr1 >= 0]
if (!length(train_xr1)) stop("No finite training XR1 premiums for risk cutpoints.")
G <- as.integer(cfg$report_risk_groups)
if (G < 2L) stop("report_risk_groups must be >= 2")
cut_raw <- as.numeric(stats::quantile(train_xr1,
        probs = seq(0, 1, length.out = G + 1L), names = FALSE, na.rm = TRUE))
cut_raw <- unique(cut_raw)
if (length(cut_raw) < 2L) {
  warning("All training XR1 premiums are equal; using a single 2010 risk group.")
  oos_breaks <- c(-Inf, Inf)
} else {
  oos_breaks <- cut_raw
  oos_breaks[1L] <- -Inf
  oos_breaks[length(oos_breaks)] <- Inf
}
oos_risk_group <- as.integer(cut(holdout_xr1,
    breaks = oos_breaks, include.lowest = TRUE, labels = FALSE))
write.csv(data.frame(Lower = head(oos_breaks, -1), Upper = tail(oos_breaks, -1)),
          oos_report_file("Figure2_TrainingXR1_RiskGroup_Breaks.csv"), row.names = FALSE)

oos_risk_parts <- list(); kk <- 0L
for (i in seq_len(nrow(oos_main_pairs))) {
  pp <- oos_main_pairs[i, , drop = FALSE]
  p <- oos_predictions[[pp$Model]]$Limited
  for (z in sort(unique(oos_risk_group[oos_valid]))) {
    ix <- which(oos_valid & oos_risk_group == z)
    if (!length(ix)) next
    kk <- kk + 1L
    oos_risk_parts[[kk]] <- cbind(
      pp[, c("Case", "Framework", "Scheme", "Model"), drop = FALSE],
      Risk_group = z, Baseline_mean_premium = mean(holdout_xr1[ix]),
      oos_scores(p[ix], oos_pol$S_limited[ix])
    )
  }
}
oos_risk_scores <- do.call(rbind, oos_risk_parts)
rownames(oos_risk_scores) <- NULL
write.csv(oos_risk_scores,
          oos_report_file("Figure2_BaselineRiskGroup_Calibration.csv"), row.names = FALSE)

grDevices::pdf(oos_report_file("Figure2_Observed_vs_Predicted_ByRiskGroup.pdf"),
               width = 11, height = 8.5)
op <- par(no.readonly = TRUE)
par(mfrow = c(3, 2), mar = c(4.1, 4.6, 2.6, .9), oma = c(.9, 0, 1.1, 0))
for (g in oos_groups) {
  gg <- oos_risk_scores[oos_risk_scores$Framework == g, , drop = FALSE]
  yr <- range(c(gg$Observed_mean, gg$Predicted_mean), finite = TRUE)
  ymax <- if (length(yr) && all(is.finite(yr))) max(1, yr[2L] * 1.12) else 1
  base <- gg[gg$Scheme == unique(gg$Scheme)[1L], , drop = FALSE]
  base <- base[order(base$Risk_group), , drop = FALSE]
  plot(base$Risk_group, base$Observed_mean, type = "b", pch = 19, lwd = 2,
       col = "#2C3E50", xlim = c(1, length(oos_breaks) - 1L),
       ylim = c(0, ymax), xlab = "Training-XR1 premium risk group",
       ylab = "Mean aggregate payment (USD)", main = g)
  ss <- unique(gg$Scheme)
  for (j in seq_along(ss)) {
    aa <- gg[gg$Scheme == ss[j], , drop = FALSE]
    aa <- aa[order(aa$Risk_group), , drop = FALSE]
    lines(aa$Risk_group, aa$Predicted_mean, type = "b", lwd = 1.5,
          col = oos_scheme_colors[j], pch = 15 + j, lty = j)
  }
  legend("topleft", legend = c("Observed", ss),
         col = c("#2C3E50", oos_scheme_colors[seq_along(ss)]),
         pch = c(19, 15 + seq_along(ss)),
         lty = c(1, seq_along(ss)), bty = "n", cex = .72)
}
plot.new()
text(.5, .63, "2010 out-of-sample comparison", cex = 1)
text(.5, .40, "Risk-group cutpoints are based on\n2006-2009 XR1 fitted premiums only", cex = .95)
mtext("2010 validation | fixed training-defined risk bins", outer = TRUE, side = 3, line = .05)
par(op); dev.off()

# -------------------- APPENDIX: COMPONENT-WISE CONTRASTS ---------------------
oos_contrasts <- do.call(rbind, lapply(c("4c", "4d"), function(ccase) {
  do.call(rbind, lapply(1:3, function(ss) {
    ba <- model_pairs[model_pairs$Case == "1" &
                        model_pairs$Scheme == paste0("XR", ss), , drop = FALSE]
    al <- model_pairs[model_pairs$Case == ccase &
                        model_pairs$Scheme == paste0("ER", ss), , drop = FALSE]
    data.frame(Case = ccase, Scheme = ss,
               Base_Model = ba$Model, Alternative_Model = al$Model)
  }))
}))
oos_contrast_parts <- lapply(seq_len(nrow(oos_contrasts)), function(i) {
  cc <- oos_contrasts[i, , drop = FALSE]
  pa <- oos_predictions[[cc$Base_Model]]$Limited
  pb <- oos_predictions[[cc$Alternative_Model]]$Limited
  both <- oos_valid & is.finite(pa) & is.finite(pb) & pa >= 0 & pb >= 0
  aa <- oos_scores(pa[both], oos_pol$S_limited[both])
  bb <- oos_scores(pb[both], oos_pol$S_limited[both])
  cbind(cc, N = aa$N, Base_OP = aa$OP_ratio, Alternative_OP = bb$OP_ratio,
        Base_Poisson_Deviance = aa$Mean_Poisson_Deviance,
        Alternative_Poisson_Deviance = bb$Mean_Poisson_Deviance,
        Delta_Poisson_Deviance = bb$Mean_Poisson_Deviance - aa$Mean_Poisson_Deviance,
        Base_MAE = aa$MAE_obs, Alternative_MAE = bb$MAE_obs,
        Delta_MAE = bb$MAE_obs - aa$MAE_obs)
})
oos_contrast_scores <- do.call(rbind, oos_contrast_parts)
write.csv(oos_contrast_scores,
          oos_report_file("Appendix_Case1_vs_Case4_Contrasts.csv"), row.names = FALSE)

# -------------------- APPENDIX: COMMON COPULA SENSITIVITY --------------------
# Reuses selection candidate fits FROM 2006-2009, does not fit candidates in 2010.
# As in the original appendix, FEM1 is fixed to isolate severity dependence.
oos_common_family_scores <- data.frame()
if (isTRUE(cfg$report_common_copula_sensitivity)) {
  oos_key_to_model <- c(
    sev_cont_complete = "SCM1", sev_cont_censored = "SCM2",
    sev_cont_truncated = "SCM3", sev_disc_complete = "SDM1",
    sev_disc_censored = "SDM2", sev_disc_truncated = "SDM3")
  sensitivity_parts <- list(); kk <- 0L
  for (key in names(oos_key_to_model)) {
    mn <- unname(oos_key_to_model[key])
    chosen <- severity_fits[[mn]]
    candidates <- sel[[key]]$fits
    candidate_codes <- unique(c(chosen$spec$fam,
                                as.integer(cfg$report_common_family_codes)))
    for (family_code in candidate_codes) {
      family <- fam_name(family_code)
      f <- candidates[[family]]
      available <- !is.null(f) && is.finite(f$AIC) &&
        !is.null(f$convergence) && f$convergence == 0L
      kk <- kk + 1L
      if (!available) {
        sensitivity_parts[[kk]] <- data.frame(
          Model = mn, Candidate_family = family,
          AIC_selected_family = fam_name(chosen$spec$fam),
          Used = FALSE, Reason = "Not fitted / nonconverged",
          oos_scores(numeric(0), numeric(0)))
      } else {
        p <- oos_freq_mean$FEM1 *
          conditional_paid_limited(f, oos_Dtheta(f), oos_dat)
        p[oos_zero_limit] <- 0
        ix <- oos_valid & is.finite(p) & p >= 0
        sensitivity_parts[[kk]] <- data.frame(
          Model = mn, Candidate_family = family,
          AIC_selected_family = fam_name(chosen$spec$fam),
          Used = TRUE, Reason = "FEM1 fixed; training candidate only",
          oos_scores(p[ix], oos_pol$S_limited[ix]))
      }
    }
  }
  oos_common_family_scores <- do.call(rbind, sensitivity_parts)
  write.csv(oos_common_family_scores,
            oos_report_file("Appendix_CommonCopula_SeverityOnly.csv"), row.names = FALSE)
}

# -------------------- APPENDIX: NUMERICAL AND GB2 DIAGNOSTICS --------------
if (exists("moments"))
  write.csv(moments,
            oos_report_file("Appendix_GB2_FirstMoment_Existence.csv"), row.names = FALSE)
if (isTRUE(cfg$capped_quad_check)) {
  id <- split(seq_len(nrow(oos_pol)), oos_pol$Deduct)
  ix <- sort(unique(unlist(lapply(id, function(z)
    z[unique(round(seq(1, length(z), length.out = min(5L, length(z)))))]),
    use.names = FALSE)))
  oos_dat_small <- lapply(oos_dat, function(v)
    if (is.matrix(v)) v[ix, , drop = FALSE] else v[ix])
  gl_lo <- gl_nodes(cfg$premium_gl_n)
  gl_hi <- gl_nodes(2L * cfg$premium_gl_n)
  oos_quadrature <- do.call(rbind, lapply(severity_fits, function(f) {
    a <- conditional_paid_limited(f, oos_Dtheta(f), oos_dat_small, gl_lo)
    b <- conditional_paid_limited(f, oos_Dtheta(f), oos_dat_small, gl_hi)
    rel <- abs(a - b) / pmax(1, abs(b))
    data.frame(Model = f$name, N_checked = length(ix),
               Max_relative_change = max(rel, na.rm = TRUE),
               Stable = all(is.finite(rel) & rel <= cfg$capped_quad_tol))
  }))
  write.csv(oos_quadrature,
            oos_report_file("Appendix_2010_CoverageLimited_Quadrature.csv"), row.names = FALSE)
  if (any(!oos_quadrature$Stable)) warning("Some 2010 limited premium quadratures are unstable. See appendix CSV.")
}

# -------------------- IN-SAMPLE vs OUT-OF-SAMPLE TABLE -----------------------
if (exists("manuscript_report") && !is.null(manuscript_report$table3)) {
  tr <- manuscript_report$table3
  fields <- c("N", "OP_ratio", "RelBias_pct", "Mean_Poisson_Deviance", "MAE_obs", "RMSE_obs")
  tr <- tr[, c("Model", fields), drop = FALSE]
  names(tr)[-1L] <- paste0("Train_2006_2009_", fields)
  te <- oos_table3[, c("Model", fields), drop = FALSE]
  names(te)[-1L] <- paste0("Test_2010_", fields)
  oos_train_test <- merge(tr, te, by = "Model", sort = FALSE)
  oos_train_test <- oos_train_test[match(oos_main_pairs$Model, oos_train_test$Model), ]
  write.csv(oos_train_test, oos_report_file("Table_Train_vs_2010.csv"), row.names = FALSE)
}

# -------------------- 2010 DATA-AUDIT AND RESULTS ARCHIVE -------------------
oos_audit <- data.frame(
  Training_years = paste(sort(unique(pol$Year)), collapse = ","),
  Test_year = oos_year,
  Train_policy_years = nrow(pol),
  Test_policy_years = nrow(oos_pol),
  Test_claim_records = if (exclude_max)
    sum(oos2010_claims$PolicyNum != excluded_policy$PolicyNum)
    else nrow(oos2010_claims),
  Excluded_policy_years = as.integer(exclude_max),
  Excluded_PolicyNum = if (exclude_max) as.character(excluded_policy$PolicyNum) else NA_character_,
  Full_2010_policy_years = nrow(full_oos_pol),
  Complete_case_test_policy_years = n_oos,
  Test_zero_limited_payment_policy_years = sum(oos_pol$S_limited == 0),
  Test_zero_coverage_above_deductible = sum(oos_zero_limit),
  Coverage_conversion_to_dollars = cfg$coverage_to_dollars,
  Modeling_coverage_interpretation = "Per-occurrence ground-up limit (assumed)",
  Copulas_reselected_on_2010 = FALSE,
  Models_refitted_on_2010 = FALSE,
  stringsAsFactors = FALSE
)
write.csv(oos_audit, oos_file("OutOfSample_2010_Data_Audit.csv"), row.names = FALSE)
write.csv(model_pairs, oos_file("Premium_Model_Pairs_TrainingFixed.csv"), row.names = FALSE)

notes2010 <- c(
  "OUT-OF-SAMPLE EMPIRICAL VALIDATION: YEAR 2010",
  "Training: policy years 2006-2009; validation: policy year 2010.",
  if (exclude_max) "POST-HOC SENSITIVITY: largest 2010 paid-loss policy-year excluded from evaluation only."
    else "PRIMARY OOS EVALUATION: all 2010 policy-years retained.",
  "All model and copula parameters are frozen at 2006-2009 estimates.",
  "No selection, calibration or parameter fitting occurs on the 2010 outcomes.",
  paste0("Common complete-case evaluation n=", n_oos,
         " of ", nrow(oos_pol), " policy-years (all nine headline models)."),
  "Primary loss: per-occurrence (min(Y, C_i)-D_i)_+ summed within policy-year.",
  paste0("C_i = Coverage_i * ", cfg$coverage_to_dollars,
         " dollars; contract interpretation is assumed, not verified."),
  "O/P = sum(observed losses) / sum(predicted pure premiums).",
  "Poisson deviance, MAE, RMSE and relative portfolio bias use zero-loss policies.",
  "RMSRE_positive excludes zero observed limited losses; N_positive is reported.",
  "O/P = NA for observed=predicted=0; O/P = Inf for observed>0 and predicted=0.",
  "Figure 1: policy-cluster bootstrap CIs condition on fixed fitted predictions.",
  "Figure 1B: cases 1, 4c, 4d in each observation scheme.",
  "Figure 2: test groups use only TRAINING-period XR1 risk cutpoints.",
  "Tables 1-2 are fixed training copula specifications: NOT test-year AIC fits.",
  "Appendix common copulas uses cached training candidates and FEM1 fixed.",
  "If some test predictions are invalid, a common subset is used for all nine models.",
  "Out-of-sample does not mean causally identified; frequency-severity factorization remains assumed."
)
writeLines(notes2010, oos_report_file("README_OutOfSample_2010.txt"))

oos_2010_results <- list(Excluded_max_policy = exclude_max,
  Excluded_policy = excluded_policy,
  train_years = 2006:2009, test_year = 2010L,
  audit = oos_audit, valid_common = which(oos_valid),
  model_pairs = model_pairs, table3 = oos_table3,
  appendix_metrics = oos_all_scores,
  deductible_op = oos_ded_scores,
  risk_group_calibration = oos_risk_scores,
  component_contrasts = oos_contrast_scores,
  common_family_sensitivity = oos_common_family_scores,
  predictions = oos_predictions,
  individual_policies = oos_pol
)
saveRDS(oos_2010_results, oos_file("OutOfSample_2010_Results.rds"))
cat("\n2010 out-of-sample Table 3 (all nine headline models):\n")
print(oos_table3[, c("Framework", "Scheme", "N", "OP_ratio", "RelBias_pct",
                       "Mean_Poisson_Deviance", "MAE_obs", "RMSE_obs")],
      row.names = FALSE, digits = 4)
cat("\n2010 evaluation complete. Results saved under: ",
    normalizePath(oos_dir), "\n", sep = "")

# The two OOS branches share the same frozen 2006-2009 model estimates.
  invisible(list(table3 = oos_table3, full_policy_data = full_oos_pol,
                 freq_mean = oos_freq_mean, sev_limited = oos_sev_limited,
                 sev_unlimited = oos_sev_unlimited,
                 excluded_policy = excluded_policy, directory = oos_dir))
}

if (isTRUE(cfg$run_oos_full) || isTRUE(cfg$run_oos_exclude_max)) {
  if (!isTRUE(cfg$run_oos_full) && isTRUE(cfg$run_oos_exclude_max))
    warning("Exclusion requested without full OOS; full analysis is recommended")
  if (isTRUE(cfg$run_oos_full))
    oos_full_result <- run_oos_2010(FALSE)
  if (isTRUE(cfg$run_oos_exclude_max))
    oos_excluded_result <- run_oos_2010(TRUE, reuse_full = if (isTRUE(cfg$run_oos_full)) oos_full_result else NULL)
  if (isTRUE(cfg$run_oos_full) && isTRUE(cfg$run_oos_exclude_max)) {
    aa <- oos_full_result$table3
    bb <- oos_excluded_result$table3
    keys <- c("Case", "Scheme", "Model")
    keep <- c(keys, "N", "OP_ratio", "RelBias_pct",
              "Mean_Poisson_Deviance", "MAE_obs", "RMSE_obs")
    comparison <- merge(aa[, keep], bb[, keep], by = keys,
                        suffixes = c("_Full", "_ExcludeMax"), sort = FALSE)
    comparison <- comparison[match(model_pairs$Model, comparison$Model), , drop = FALSE]
    write.csv(comparison,
      file.path(cfg$output_dir, "Table_2010_Full_vs_ExcludeMax.csv"), row.names = FALSE)
  }
  cat("\nBoth 2010 OOS configurations completed with frozen 2006-2009 fits.\n")
}



######################################################################################################################


################################################################################
# OUT-OF-SAMPLE VALIDATION: POLICY YEAR 2010
# MERGED INTO THE EMPIRICAL 2006-2009 ESTIMATION SCRIPT.
# Alternatively, source this file after restoring the completed training workspace.
#
# NO refitting, NO copula reselection, NO Godambe recomputation.
# 2010 outcomes are accessed ONLY for observed-loss evaluation and scoring.
# Nine headline model pairs are evaluated: Cases 1, 4c and 4d.
# The primary loss: sum_j (min(Y_ij, C_i) - D_i)_+, C_i = Coverage_i * scale.
################################################################################

# ================================================================
# 2010 OOS diagnostics across the continuous deductible ratio D/C
# Uses frozen 2006-2009 predictions, without arbitrary ratio bins.
# Install helper BEFORE run_oos_2010(), and call it just after
# oos_valid has been calculated inside run_oos_2010().
# ================================================================

oos_deductible_ratio_diagnostics <- function(
    oos_pol, oos_predictions, model_pairs, oos_valid, oos_file,
    exclude_max = FALSE, coverage_to_dollars = 1e6,
    ratio_xlim_pct = NULL, bandwidth = NULL, min_ess = 15) {
  
  stopifnot(nrow(oos_pol) == length(oos_valid),
            all(model_pairs$Model %in% names(oos_predictions)),
            length(coverage_to_dollars) == 1L,
            is.finite(coverage_to_dollars), coverage_to_dollars > 0)
  
  sel <- which(oos_valid)
  if (!length(sel))
    stop("No common valid 2010 policies for deductible-ratio diagnostic.")
  
  # Deduct is in USD; Coverage is in millions of USD.
  # The dimensionless ratio must use coverage converted to dollars.
  coverage_dollars <- oos_pol$Coverage[sel] * coverage_to_dollars
  r <- as.numeric(oos_pol$Deduct[sel] / coverage_dollars)
  y <- as.numeric(oos_pol$S_limited[sel])
  if (any(!is.finite(r) | r <= 0) || any(!is.finite(y) | y < 0))
    stop("Invalid 2010 deductible ratios or losses.")
  
  sample_tag <- if (exclude_max) "ExcludeMax" else "Full"
  if (is.null(ratio_xlim_pct)) ratio_xlim_pct <- range(100 * r)
  if (length(ratio_xlim_pct) != 2L ||
      any(!is.finite(ratio_xlim_pct) | ratio_xlim_pct <= 0) ||
      ratio_xlim_pct[1] >= ratio_xlim_pct[2])
    stop("ratio_xlim_pct must contain two increasing positive percentages.")
  
  # Policy-level diagnostics for all nine models; zero-loss policies included.
  row_list <- lapply(seq_len(nrow(model_pairs)), function(j) {
    pp <- model_pairs[j, , drop = FALSE]
    p <- as.numeric(oos_predictions[[pp$Model]]$Limited[sel])
    if (length(p) != length(y) || any(!is.finite(p) | p < 0))
      stop("Invalid 2010 limited premium: ", pp$Model)
    dev <- rep(Inf, length(p))
    dev[y == 0] <- 2 * p[y == 0]
    both_pos <- y > 0 & p > 0
    dev[both_pos] <- 2 * (y[both_pos] * log(y[both_pos] / p[both_pos]) -
                            y[both_pos] + p[both_pos])
    dev[is.finite(dev)] <- pmax(dev[is.finite(dev)], 0)
    
    data.frame(
      Case = pp$Case, Framework = pp$Framework,
      Scheme = pp$Scheme, Model = pp$Model,
      PolicyNum = oos_pol$PolicyNum[sel], Year = oos_pol$Year[sel],
      Deductible = oos_pol$Deduct[sel],
      Coverage = oos_pol$Coverage[sel],
      Coverage_Dollars = coverage_dollars,
      DeductibleRatio = r, DeductibleRatio_pct = 100 * r,
      Observed = y, Predicted = p,
      AbsoluteError = abs(y - p), SquaredError = (y - p)^2,
      PoissonDeviance = dev, stringsAsFactors = FALSE
    )
  })
  policy <- do.call(rbind, row_list)
  rownames(policy) <- NULL
  write.csv(policy,
            oos_file(paste0("DeductibleRatio_PolicyLevel_2010_", sample_tag, ".csv")),
            row.names = FALSE)
  
  # Aggregate only EXACT duplicate ratios; no arbitrary binning.
  exact <- do.call(rbind, lapply(split(policy, policy$Model), function(d) {
    do.call(rbind, lapply(split(seq_len(nrow(d)), d$DeductibleRatio),
                          function(ix) {
                            z <- d[ix, , drop = FALSE]
                            sp <- sum(z$Predicted); sy <- sum(z$Observed)
                            data.frame(
                              Model = z$Model[1], Framework = z$Framework[1],
                              Scheme = z$Scheme[1], DeductibleRatio = z$DeductibleRatio[1],
                              DeductibleRatio_pct = z$DeductibleRatio_pct[1],
                              N = nrow(z), ObservedTotal = sy, PredictedTotal = sp,
                              MAE = mean(z$AbsoluteError),
                              RMSE = sqrt(mean(z$SquaredError)),
                              MPD = mean(z$PoissonDeviance),
                              OP_ratio = if (sp > 0) sy / sp else NA_real_
                            )
                          })
    )
  }))
  exact <- exact[order(exact$Model, exact$DeductibleRatio), ]
  rownames(exact) <- NULL
  write.csv(exact,
            oos_file(paste0("DeductibleRatio_ExactValues_2010_", sample_tag, ".csv")),
            row.names = FALSE)
  
  # Kernel smoothing on log10(dimensionless deductible ratio).
  # Full and ExcludeMax use the same range and bandwidth.
  x <- log10(r)
  if (is.null(bandwidth))
    bandwidth <- max(0.08, 0.20 * stats::IQR(x))
  if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0)
    stop("bandwidth must be a positive log10-ratio number.")
  
  grid <- 10^seq(log10(ratio_xlim_pct[1] / 100),
                 log10(ratio_xlim_pct[2] / 100),
                 length.out = 160L)
  
  smooth_one <- function(d) {
    xd <- log10(d$DeductibleRatio)
    out <- lapply(grid, function(g) {
      logw <- -0.5 * ((xd - log10(g)) / bandwidth)^2
      w <- exp(logw - max(logw))  # avoid numerical underflow
      w <- w / sum(w)
      ess <- 1 / sum(w^2)
      wp <- sum(w * d$Predicted)
      wy <- sum(w * d$Observed)
      # Any positive-weight infinite Poisson deviance makes the mean infinite.
      dev <- if (any(is.infinite(d$PoissonDeviance) & w > 0)) {
        Inf
      } else {
        sum(w * d$PoissonDeviance)
      }
      data.frame(
        DeductibleRatio = g, DeductibleRatio_pct = 100 * g,
        EffectiveN = ess,
        MAE = sum(w * d$AbsoluteError),
        RMSE = sqrt(sum(w * d$SquaredError)),
        MPD = dev,
        OP_ratio = if (wp > 0) wy / wp else NA_real_,
        WeightedObserved = wy, WeightedPredicted = wp
      )
    })
    z <- do.call(rbind, out)
    z$Model <- d$Model[1]
    z$Framework <- d$Framework[1]
    z$Scheme <- d$Scheme[1]
    z
  }
  
  smooth <- do.call(rbind, lapply(split(policy, policy$Model), smooth_one))
  rownames(smooth) <- NULL
  write.csv(smooth,
            oos_file(paste0("DeductibleRatio_Smoothed_2010_", sample_tag, ".csv")),
            row.names = FALSE)
  
  # Four PDFs: XR1/2/3 versus the corresponding continuous ER1/2/3.
  plot_metric <- function(metric) {
    file <- oos_file(paste0(
      "DeductibleRatio_", metric,
      "_XR_vs_ContinuousER_2010_", sample_tag, ".pdf"))
    grDevices::pdf(file, width = 9, height = 10)
    old_par <- par(no.readonly = TRUE)
    on.exit({par(old_par); grDevices::dev.off()}, add = TRUE)
    par(mfrow = c(3, 1), mar = c(4.5, 5.2, 3.0, 1.2))
    
    for (k in 1:3) {
      a <- smooth[smooth$Framework == "XR_marginal" &
                    smooth$Scheme == paste0("XR", k), , drop = FALSE]
      b <- smooth[smooth$Framework == "Case4_both_cont" &
                    smooth$Scheme == paste0("ER", k), , drop = FALSE]
      a <- a[is.finite(a[[metric]]) & a$EffectiveN >= min_ess, , drop = FALSE]
      b <- b[is.finite(b[[metric]]) & b$EffectiveN >= min_ess, , drop = FALSE]
      
      if (!nrow(a) || !nrow(b)) {
        plot.new()
        title(main = paste0("Scenario ", k, ": insufficient local support"))
        next
      }
      yl <- range(c(a[[metric]], b[[metric]]))
      if (diff(yl) < 1e-12)
        yl <- yl + c(-0.1, 0.1) * max(1, abs(yl[1]))
      
      plot(a$DeductibleRatio_pct, a[[metric]],
           type = "l", lwd = 2, col = "steelblue4", log = "x",
           xlim = ratio_xlim_pct, ylim = yl,
           xlab = "Deductible / Coverage (%) [log scale]",
           ylab = metric,
           main = paste0("2010 ", sample_tag, ": XR", k,
                         " vs continuous ER", k))
      lines(b$DeductibleRatio_pct, b[[metric]],
            lwd = 2, col = "firebrick3", lty = 2)
      if (metric == "OP_ratio")
        abline(h = 1, lty = 3, col = "gray50")
      legend("topright",
             legend = c(paste0("XR", k), paste0("continuous ER", k)),
             col = c("steelblue4", "firebrick3"),
             lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
    }
    invisible(file)
  }
  invisible(lapply(c("MAE", "RMSE", "MPD", "OP_ratio"), plot_metric))
  
  writeLines(c(
    paste0("2010 OOS deductible-ratio diagnostic: ", sample_tag),
    "Training fits fixed at 2006-2009.",
    "Deduct in dollars; Coverage in millions of dollars.",
    paste0("Coverage conversion to dollars: ", coverage_to_dollars),
    "Dimensionless ratio = Deduct / (Coverage * coverage_to_dollars).",
    "Horizontal plot axis is 100 * ratio, in percent (logarithmic).",
    paste0("Shared full-sample x-axis percentage range: ",
           paste(signif(ratio_xlim_pct, 6), collapse = " to ")),
    "Full and ExcludeMax use identical smoothing grids and bandwidths.",
    "No quantile bins or arbitrary deductible-ratio groups.",
    paste0("Gaussian bandwidth in log10 units: ", signif(bandwidth, 6)),
    paste0("Minimum effective sample size displayed: ", min_ess),
    "Local O/P = weighted observed / weighted predicted.",
    "Zero-loss policies are retained. Curves are descriptive, not causal."
  ), oos_file(paste0("DeductibleRatio_README_2010_", sample_tag, ".txt")))
  
  cat("  2010 deductible-ratio diagnostics written (", sample_tag,
      "; n=", length(sel),
      "; scale=", coverage_to_dollars,
      "; x-axis %=(", paste(signif(ratio_xlim_pct, 4), collapse = ", "),
      "); bandwidth=", signif(bandwidth, 4), ")\n", sep = "")
  invisible(list(policy = policy, exact = exact, smooth = smooth,
                 bandwidth = bandwidth, ratio_xlim_pct = ratio_xlim_pct))
}


# ================================================================
# 2010 OOS: metric curves over D/C WITHIN EACH DEDUCTIBLE LEVEL
# Model fits/predictions are supplied from the pre-existing OOS run.
# No deductible-ratio bins are created.
# ================================================================

oos_ratio_by_deductible_level <- function(
    oos_pol, full_oos_pol, oos_predictions, model_pairs,
    oos_valid, oos_file, coverage_to_dollars,
    exclude_max = FALSE, bandwidth = NULL, min_ess = 5L,
    n_grid = 120L) {
  
  stopifnot(length(oos_valid) == nrow(oos_pol),
            is.finite(coverage_to_dollars), coverage_to_dollars > 0)
  ix <- which(oos_valid)
  if (!length(ix)) stop("No complete 2010 OOS observations")
  
  # Full-2010 deductible levels define the panels and each panel's x-axis.
  # Coverage is in MILLIONS of USD, while Deduct is already in USD.
  full_ratio <- full_oos_pol$Deduct /
    (full_oos_pol$Coverage * coverage_to_dollars)
  ratio <- oos_pol$Deduct[ix] /
    (oos_pol$Coverage[ix] * coverage_to_dollars)
  if (any(!is.finite(c(full_ratio, ratio))) ||
      any(c(full_ratio, ratio) <= 0))
    stop("Invalid D/C ratios")
  
  levels <- sort(unique(full_oos_pol$Deduct))
  if (is.null(bandwidth))
    bandwidth <- max(0.08, 0.20 * stats::IQR(log10(full_ratio)))
  if (!is.finite(bandwidth) || bandwidth <= 0)
    stop("Invalid kernel bandwidth")
  
  keep_models <- model_pairs[
    model_pairs$Framework %in% c("XR_marginal", "Case4_both_cont"),
    , drop = FALSE
  ]
  if (nrow(keep_models) != 6L)
    stop("Expected XR1-XR3 and continuous ER1-ER3")
  
  y <- oos_pol$S_limited[ix]
  data_rows <- lapply(seq_len(nrow(keep_models)), function(j) {
    pp <- keep_models[j, , drop = FALSE]
    p <- as.numeric(oos_predictions[[pp$Model]]$Limited[ix])
    if (length(p) != length(y) || any(!is.finite(p) | p < 0))
      stop("Invalid premiums: ", pp$Model)
    dev <- rep(Inf, length(y))
    dev[y == 0] <- 2 * p[y == 0]
    ok <- y > 0 & p > 0
    dev[ok] <- 2 * (y[ok] * log(y[ok] / p[ok]) - y[ok] + p[ok])
    dev[is.finite(dev)] <- pmax(dev[is.finite(dev)], 0)
    data.frame(
      Model = pp$Model, Framework = pp$Framework,
      Scheme = pp$Scheme, Deductible = oos_pol$Deduct[ix],
      Ratio_pct = 100 * ratio, Observed = y, Predicted = p,
      AE = abs(y - p), SE = (y - p)^2, Deviance = dev,
      stringsAsFactors = FALSE
    )
  })
  dat <- do.call(rbind, data_rows)
  tag <- if (exclude_max) "ExcludeMax" else "Full"
  write.csv(dat, oos_file(paste0(
    "DeductibleRatio_ByLevel_PolicyLevel_2010_", tag, ".csv")),
    row.names = FALSE)
  
  # These are actual deductible values, not ratio-quantile buckets.
  # Same grids and bandwidth for Full and ExcludeMax.
  grid_list <- lapply(levels, function(d) {
    fr <- full_ratio[full_oos_pol$Deduct == d]
    rr <- range(fr)
    if (rr[1] == rr[2]) return(100 * rr[1])
    100 * 10^seq(log10(rr[1]), log10(rr[2]), length.out = n_grid)
  })
  names(grid_list) <- as.character(levels)
  
  count <- do.call(rbind, lapply(levels, function(d) {
    data.frame(
      Deductible = d,
      FullN = sum(full_oos_pol$Deduct == d),
      EvaluationN = sum(oos_pol$Deduct[ix] == d),
      UniqueRatios = length(unique(ratio[oos_pol$Deduct[ix] == d])),
      RatioMin_pct = min(grid_list[[as.character(d)]]),
      RatioMax_pct = max(grid_list[[as.character(d)]])
    )
  }))
  write.csv(count, oos_file(paste0(
    "DeductibleRatio_ByLevel_Counts_2010_", tag, ".csv")),
    row.names = FALSE)
  
  smooth_parts <- list()
  for (d in levels) {
    grid <- grid_list[[as.character(d)]]
    dd <- dat[dat$Deductible == d, , drop = FALSE]
    if (!nrow(dd)) next
    for (model in unique(dd$Model)) {
      z <- dd[dd$Model == model, , drop = FALSE]
      logz <- log10(z$Ratio_pct)
      for (g in grid) {
        logw <- -0.5 * ((logz - log10(g)) / bandwidth)^2
        w <- exp(logw - max(logw))
        w <- w / sum(w)
        wp <- sum(w * z$Predicted)
        wy <- sum(w * z$Observed)
        mpd <- if (any(is.infinite(z$Deviance) & w > 0)) Inf
        else sum(w * z$Deviance)
        smooth_parts[[length(smooth_parts) + 1L]] <- data.frame(
          Deductible = d, Model = model,
          Framework = z$Framework[1], Scheme = z$Scheme[1],
          Ratio_pct = g, N = nrow(z), EffectiveN = 1 / sum(w^2),
          MAE = sum(w * z$AE),
          RMSE = sqrt(sum(w * z$SE)), MPD = mpd,
          OP_ratio = if (wp > 0) wy / wp else NA_real_,
          WeightedObserved = wy, WeightedPredicted = wp
        )
      }
    }
  }
  curves <- if (length(smooth_parts)) do.call(rbind, smooth_parts)
  else data.frame()
  write.csv(curves, oos_file(paste0(
    "DeductibleRatio_ByLevel_Smoothed_2010_", tag, ".csv")),
    row.names = FALSE)
  
  # One PDF per metric: deductible levels = rows; schemes = columns.
  # Each panel overlays XRk and continuous ERk.
  # For a given deductible level, the x-axis is identical across
  # models, metrics, and Full/ExcludeMax runs.
  draw_metric <- function(metric) {
    outfile <- oos_file(paste0(
      "DeductibleRatio_ByLevel_", metric, "_2010_", tag, ".pdf"))
    grDevices::pdf(outfile, width = 15, height = max(7, 2.8 * length(levels)))
    oldpar <- par(no.readonly = TRUE)
    on.exit({par(oldpar); grDevices::dev.off()}, add = TRUE)
    par(mfrow = c(length(levels), 3L),
        mar = c(3.9, 4.4, 2.8, 1), oma = c(0, 0, 1.5, 0))
    
    for (d in levels) {
      xfull <- grid_list[[as.character(d)]]
      for (k in 1:3) {
        zr <- curves[curves$Deductible == d &
                       curves$Framework == "XR_marginal" &
                       curves$Scheme == paste0("XR", k), , drop = FALSE]
        ze <- curves[curves$Deductible == d &
                       curves$Framework == "Case4_both_cont" &
                       curves$Scheme == paste0("ER", k), , drop = FALSE]
        zr <- zr[is.finite(zr[[metric]]) & zr$EffectiveN >= min_ess,
                 , drop = FALSE]
        ze <- ze[is.finite(ze[[metric]]) & ze$EffectiveN >= min_ess,
                 , drop = FALSE]
        good <- c(zr[[metric]], ze[[metric]])
        good <- good[is.finite(good)]
        label <- paste0("D=$", format(d, big.mark = ","),
                        " | XR", k, " vs cont. ER", k)
        if (!length(good) || length(unique(xfull)) < 2L) {
          plot.new()
          title(main = label)
          text(.5, .5, "Insufficient ratio variation / support", cex = .8)
          next
        }
        yl <- range(good)
        if (diff(yl) < 1e-12)
          yl <- yl + c(-.1, .1) * max(1, abs(yl[1]))
        plot(NA_real_, NA_real_, type = "n", log = "x",
             xlim = range(xfull), ylim = yl,
             xlab = "Deductible / Coverage (%)", ylab = metric,
             main = label)
        if (nrow(zr)) lines(zr$Ratio_pct, zr[[metric]],
                            col = "steelblue4", lwd = 2)
        if (nrow(ze)) lines(ze$Ratio_pct, ze[[metric]],
                            col = "firebrick3", lwd = 2, lty = 2)
        if (metric == "OP_ratio") abline(h = 1, col = "gray50", lty = 3)
        legend("topright", c(paste0("XR", k), paste0("ER", k, " cont.")),
               col = c("steelblue4", "firebrick3"),
               lty = c(1, 2), lwd = 2, bty = "n", cex = .65)
      }
    }
    mtext(paste0("2010 OOS ", tag,
                 " | Continuous D/C within each deductible level"),
          side = 3, outer = TRUE, cex = 1)
    invisible(outfile)
  }
  if (nrow(curves)) {
    invisible(lapply(c("MAE", "RMSE", "MPD", "OP_ratio"), draw_metric))
  }
  invisible(list(curves = curves, counts = count, policy = dat))
}

run_oos_2010 <- function(exclude_max = FALSE, reuse_full = NULL) {
  cat("\n================ OUT-OF-SAMPLE: YEAR 2010 | ",
      if (exclude_max) "EXCLUDE MAX" else "FULL SAMPLE",
      " ================\n", sep = "")
  
  oos_year <- 2010L
  needed <- c("cfg", "poldat", "claimdat", "pol", "clm", "X_pol", "fml",
              "model_pairs", "frequency_fits", "severity_fits", "fit_Dc", "fit_Dd",
              "prepare_eval", "conditional_N_mean", "conditional_paid_limited",
              "conditional_paid_unlimited", "gl_nodes", "raw_predictions")
  missing_items <- needed[!vapply(needed, exists, logical(1), inherits = TRUE)]
  if (length(missing_items)) stop("Run the existing empirical script first. Missing: ",
                                  paste(missing_items, collapse = ", "))
  if (!identical(as.integer(cfg$years_in), 2006:2009)) {
    stop("Expected fitted training years cfg$years_in = 2006:2009; current setting is: ",
         paste(cfg$years_in, collapse = ", "))
  }
  if (any(pol$Year == oos_year) || any(clm$Year == oos_year))
    stop("2010 data appear in training objects. Out-of-sample assessment is invalid.")
  if (length(cfg$holdout_years))
    stop("This appended analysis expects the original 2006-2009 in-sample fit, with holdout_years = integer(0).")
  
  # Keep original cfg, eval_pol, eval_dat, and in-sample predictions unchanged.
  oos_dir <- file.path(cfg$output_dir,
                       if (exclude_max) "out_of_sample_2010_exclude_max_loss"
                       else "out_of_sample_2010")
  oos_report_dir <- file.path(oos_dir, "manuscript_reporting")
  dir.create(oos_report_dir, recursive = TRUE, showWarnings = FALSE)
  oos_file <- function(name) file.path(oos_dir, name)
  oos_report_file <- function(name) file.path(oos_report_dir, name)
  
  oos2010_policies <- poldat[poldat$Year == oos_year, , drop = FALSE]
  oos2010_claims <- claimdat[claimdat$Year == oos_year, , drop = FALSE]
  if (!exclude_max) {
    # Descriptive 2010-only statistics, unmodified by the sensitivity exclusion.
    pos_2010 <- pmax(0, oos2010_claims$LossBeforeDeductible - oos2010_claims$Deduct)
    nn2010 <- as.integer(oos2010_claims$LossBeforeDeductible > oos2010_claims$Deduct)
    count_2010 <- if (nrow(oos2010_claims))
      stats::aggregate(nn2010,
                       by = list(PolicyNum = oos2010_claims$PolicyNum, Year = oos2010_claims$Year), FUN = sum)
    else data.frame(PolicyNum = character(0), Year = integer(0), x = integer(0))
    keys_pol <- paste(oos2010_policies$PolicyNum, oos2010_policies$Year)
    keys_claim <- paste(count_2010$PolicyNum, count_2010$Year)
    nstar <- count_2010$x[match(keys_pol, keys_claim)]
    nstar[is.na(nstar)] <- 0
    stat2010 <- function(x) c(Mean = if (length(x)) mean(x) else NA_real_,
                              SD = if (length(x) > 1L) sd(x) else NA_real_)
    summary_2010 <- t(rbind(
      N = stat2010(oos2010_policies$Freq), Nstar = stat2010(nstar),
      Y = stat2010(oos2010_claims$LossBeforeDeductible),
      Y_censored = stat2010(pos_2010),
      Y_truncated = stat2010(pos_2010[pos_2010 > 0])))
    dir.create(file.path(cfg$output_dir, "out_of_sample_2010"), showWarnings = FALSE)
    write.csv(summary_2010, file.path(cfg$output_dir, "out_of_sample_2010", "Summary_Statistics_2010.csv"))
  }
  if (!nrow(oos2010_policies)) stop("There are no 2010 policies in poldat.")
  if (anyDuplicated(paste(oos2010_policies$PolicyNum, oos2010_policies$Year)))
    stop("Duplicate PolicyNum-Year keys in 2010 policy records.")
  if (any(!is.finite(oos2010_policies$Coverage) | oos2010_policies$Coverage <= 0))
    stop("Invalid 2010 policy coverage; cannot construct design/premium inputs.")
  if (any(!is.finite(oos2010_policies$Deduct) | oos2010_policies$Deduct <= 0))
    stop("Invalid 2010 policy deductible.")
  if (anyNA(match(oos2010_policies$Deduct, cfg$ded_levels)))
    stop("2010 contains deductible levels outside cfg$ded_levels; do not silently drop these policies.")
  if (nrow(oos2010_claims)) {
    bad_claims <- !paste(oos2010_claims$PolicyNum, oos2010_claims$Year) %in%
      paste(oos2010_policies$PolicyNum, oos2010_policies$Year)
    if (any(bad_claims)) stop(sum(bad_claims),
                              " 2010 claims have no matching 2010 policy record. Resolve linkage before scoring.")
    if (any(!is.finite(oos2010_claims$LossBeforeDeductible) |
            oos2010_claims$LossBeforeDeductible < 0))
      stop("Nonfinite or negative 2010 ground-up claim amounts.")
  } else {
    warning("There are no 2010 claim records: all observed 2010 aggregate losses will be zero.")
  }
  
  # Exact same preprocessing/payment definition used in original evaluation.
  oos_pol <- prepare_eval(oos2010_policies, oos2010_claims)
  # POST-HOC sensitivity exclusion, same fixed 2006-2009 parameter estimates.
  # Exclude all claims for the policy-year with largest limited aggregate S.
  full_oos_pol <- oos_pol
  oos_sel <- seq_len(nrow(oos_pol))
  excluded_policy <- NULL
  if (exclude_max) {
    if (nrow(oos_pol) < 2L) stop("Need at least two 2010 policy-years")
    if (any(!is.finite(oos_pol$S_limited)))
      stop("Invalid 2010 observed aggregate losses")
    imax <- which(oos_pol$S_limited == max(oos_pol$S_limited))
    if (length(imax) != 1L) stop("Maximum 2010 loss is tied; specify policy manually")
    excluded_policy <- oos_pol[imax, , drop = FALSE]
    oos_sel <- seq_len(nrow(oos_pol))[-imax]
    excluded_policy$Exclusion_reason <-
      "Largest 2010 observed aggregate coverage-limited payment (post hoc)"
    write.csv(excluded_policy,
              oos_file("Excluded_2010_Extreme_PolicyYear.csv"), row.names = FALSE)
    oos_pol <- oos_pol[-imax, , drop = FALSE]
    rownames(oos_pol) <- NULL
    cat("  Excluded one 2010 policy-year: ",
        as.character(excluded_policy$PolicyNum),
        "; observed limited S = ", excluded_policy$S_limited, "\n", sep = "")
  }
  oos_X <- model.matrix(fml, data = oos_pol)
  if (!identical(colnames(oos_X), colnames(X_pol)))
    stop("2010 design-matrix columns differ from training design-matrix columns.")
  
  oos_dat <- list(
    X = oos_X,
    d = oos_pol$Deduct,
    coverage_limit = oos_pol$CoverageLimit_Dollars,
    dratio = oos_pol$dedratio,
    dlev = oos_pol$discDed,
    n_cen = oos_pol$n_ground,
    n_tru = oos_pol$n_paid,
    unit = seq_len(nrow(oos_pol))
  )
  oos_zero_limit <- oos_pol$CoverageBelowOrEqualDeductible
  stopifnot(all(oos_pol$S_limited[oos_zero_limit] == 0),
            all(oos_pol$S_limited >= 0),
            all(oos_pol$S_limited <= oos_pol$S_raw + 1e-7))
  
  # In this step the fitted theta_D, theta_Y, theta_N, and copula families
  # are fixed at their 2006-2009 estimates. Outcomes are not input to predictors.
  oos_Dtheta <- function(fit) switch(fit$spec$dedtype,
                                     cont = fit_Dc$par, disc = fit_Dd$par, none = NULL)
  if (is.null(reuse_full)) {
    cat("Predicting 2010 conditional frequency means (fixed 2006-2009 fits)...\n")
    oos_freq_mean <- lapply(frequency_fits[premium_frequency_models], function(f)
      as.numeric(conditional_N_mean(f, oos_Dtheta(f), oos_dat)))
    cat("Predicting 2010 coverage-limited severity payments (fixed fits)...\n")
    oos_sev_limited <- lapply(severity_fits, function(f)
      as.numeric(conditional_paid_limited(f, oos_Dtheta(f), oos_dat)))
  } else {
    # Predictions come from FULL 2010 inputs; only evaluation rows are removed.
    # No re-estimation AND no repeated prediction integration for the excluded case.
    stopifnot(exclude_max,
              length(oos_sel) == nrow(oos_pol),
              identical(as.character(full_oos_pol$PolicyNum),
                        as.character(reuse_full$full_policy_data$PolicyNum)))
    oos_freq_mean <- lapply(reuse_full$freq_mean, function(z) z[oos_sel])
    oos_sev_limited <- lapply(reuse_full$sev_limited, function(z) z[oos_sel])
  }
  
  if (!all(vapply(oos_sev_limited, function(v)
    length(v) == nrow(oos_pol) && all(v[oos_zero_limit] == 0), logical(1))))
    stop("2010 limited severity payments violate the contractual C <= D zero-payment condition.")
  
  # Supplementary: retain the same unlimited calculation and GB2 moment guard.
  if (!is.null(reuse_full)) {
    oos_sev_unlimited <- lapply(reuse_full$sev_unlimited, function(z) z[oos_sel])
  } else if (isTRUE(cfg$compute_unlimited)) {
    cat("Predicting supplementary unlimited severity payments for 2010...\n")
    oos_sev_unlimited <- lapply(severity_fits, function(f)
      as.numeric(conditional_paid_unlimited(f, oos_Dtheta(f), oos_dat)))
  } else {
    oos_sev_unlimited <- lapply(severity_fits, function(f)
      rep(NA_real_, nrow(oos_pol)))
  }
  
  # Preserve EXACTLY the nine Case 1 / 4c / 4d model pairings.
  oos_main_pairs <- model_pairs  # exactly nine Case 1, 4c, 4d schemes
  stopifnot(nrow(oos_main_pairs) == 9L,
            all(oos_main_pairs$Frequency %in% names(oos_freq_mean)),
            all(oos_main_pairs$Severity %in% names(oos_sev_limited)))
  
  oos_predictions <- setNames(vector("list", nrow(model_pairs)), model_pairs$Model)
  for (ii in seq_len(nrow(model_pairs))) {
    pp <- model_pairs[ii, ]
    pred_limited <- as.numeric(oos_freq_mean[[pp$Frequency]] *
                                 oos_sev_limited[[pp$Severity]])
    pred_unlimited <- as.numeric(oos_freq_mean[[pp$Frequency]] *
                                   oos_sev_unlimited[[pp$Severity]])
    pred_limited[oos_zero_limit] <- 0
    if (length(pred_limited) != nrow(oos_pol)) stop("Incorrect prediction length: ", pp$Model)
    oos_predictions[[pp$Model]] <- list(Limited = pred_limited,
                                        Unlimited = pred_unlimited)
  }
  saveRDS(oos_predictions, oos_file("Premium_Predictions_ByPolicy_2010.rds"))
  
  # Export individual policy-year predictions for reproducible paired comparisons.
  oos_prem_wide <- data.frame(
    PolicyNum = oos_pol$PolicyNum, Year = oos_pol$Year,
    Deductible = oos_pol$Deduct, Coverage = oos_pol$Coverage,
    CoverageLimit_Dollars = oos_pol$CoverageLimit_Dollars,
    CoverageBelowOrEqualDeductible = oos_zero_limit,
    Observed_Aggregate_Limited = oos_pol$S_limited,
    Observed_Aggregate_Unlimited = oos_pol$S_raw,
    stringsAsFactors = FALSE
  )
  for (nm in names(oos_predictions))
    oos_prem_wide[[paste0("Premium_Limited__", nm)]] <- oos_predictions[[nm]]$Limited
  write.csv(oos_prem_wide, oos_file("Premium_Predictions_ByPolicy_2010.csv"), row.names = FALSE)
  oos_prem_long <- do.call(rbind, lapply(seq_len(nrow(model_pairs)), function(i) {
    pp <- model_pairs[i, , drop = FALSE]
    data.frame(Case = pp$Case, Framework = pp$Framework,
               Scheme = pp$Scheme, Model = pp$Model,
               PolicyNum = oos_pol$PolicyNum, Year = oos_pol$Year,
               Deductible = oos_pol$Deduct,
               Observed_Aggregate_Limited = oos_pol$S_limited,
               Premium_Limited = oos_predictions[[pp$Model]]$Limited,
               stringsAsFactors = FALSE)
  }))
  write.csv(oos_prem_long,
            oos_file("Premium_Predictions_ByPolicy_Long_2010.csv"), row.names = FALSE)
  
  # -------------------- IDENTICAL METRICS / ZERO CONVENTIONS --------------------
  oos_op <- function(y, p) {
    sy <- sum(y); sp <- sum(p)
    if (!is.finite(sy) || !is.finite(sp)) return(NA_real_)
    if (sp > 0) return(sy / sp)
    if (sy > 0) return(Inf)
    NA_real_  # O/P undefined when both observed and predicted totals are zero.
  }
  
  oos_poisson_deviance <- function(y, p) {
    stopifnot(length(y) == length(p),
              all(is.finite(y)), all(is.finite(p)), all(y >= 0), all(p >= 0))
    if (!length(y)) return(NA_real_)
    a <- numeric(length(y))
    a[y == 0] <- 2 * p[y == 0]
    a[y > 0 & p == 0] <- Inf
    ok <- y > 0 & p > 0
    a[ok] <- 2 * (y[ok] * log(y[ok] / p[ok]) - (y[ok] - p[ok]))
    a[is.finite(a)] <- pmax(a[is.finite(a)], 0)
    mean(a)
  }
  stopifnot(isTRUE(all.equal(oos_poisson_deviance(c(0, 12), c(0, 12)), 0)),
            isTRUE(all.equal(oos_poisson_deviance(0, 3), 6)),
            is.infinite(oos_poisson_deviance(2, 0)))
  
  oos_scores <- function(pred, obs) {
    stopifnot(length(pred) == length(obs))
    good <- is.finite(pred) & is.finite(obs) & pred >= 0 & obs >= 0
    p <- as.numeric(pred[good]); y <- as.numeric(obs[good]); n <- length(y)
    if (!n) return(data.frame(
      N = 0L, N_positive = 0L,
      Observed_total = NA_real_, Predicted_total = NA_real_,
      Observed_mean = NA_real_, Predicted_mean = NA_real_,
      OP_ratio = NA_real_, RelBias_pct = NA_real_,
      Mean_Poisson_Deviance = NA_real_, MAE_obs = NA_real_,
      RMSE_obs = NA_real_, RMSRE_positive = NA_real_))
    pos <- y > 0
    data.frame(
      N = n, N_positive = sum(pos),
      Observed_total = sum(y), Predicted_total = sum(p),
      Observed_mean = mean(y), Predicted_mean = mean(p),
      OP_ratio = oos_op(y, p),
      RelBias_pct = if (sum(y) > 0) 100 * (sum(p) - sum(y)) / sum(y) else NA_real_,
      Mean_Poisson_Deviance = oos_poisson_deviance(y, p),
      MAE_obs = mean(abs(p - y)),
      RMSE_obs = sqrt(mean((p - y)^2)),
      RMSRE_positive = if (any(pos))
        sqrt(mean(((p[pos] - y[pos]) / y[pos])^2)) else NA_real_
    )
  }
  
  # Main comparison: EXACTLY THE SAME 2010 POLICY-YEARS for all nine models.
  oos_valid <- is.finite(oos_pol$S_limited) & oos_pol$S_limited >= 0
  oos_completeness <- data.frame(Model = oos_main_pairs$Model,
                                 Invalid_predictor_rows = NA_integer_)
  for (ii in seq_len(nrow(oos_main_pairs))) {
    nm <- oos_main_pairs$Model[ii]
    p <- oos_predictions[[nm]]$Limited
    good <- is.finite(p) & p >= 0
    oos_completeness$Invalid_predictor_rows[ii] <- sum(!good)
    oos_valid <- oos_valid & good
  }
  write.csv(oos_completeness,
            oos_report_file("Appendix_Prediction_Completeness_2010.csv"), row.names = FALSE)
  n_oos <- sum(oos_valid)
  if (!n_oos) stop("No common valid 2010 policy-year observations across all nine headline models. Check completeness CSV.")
  if (n_oos < nrow(oos_pol)) warning("Headline comparison retains ", n_oos,
                                     " / ", nrow(oos_pol), " valid 2010 policy-years.")
  
  # Continuous deductible-ratio plots: executed separately for BOTH 2010 OOS runs.
  # All parameters/copula selections remain frozen at 2006-2009 training fits.
  # Use the FULL 2010 policy universe for both runs' horizontal scale
  # and Gaussian-kernel bandwidth. Exclusion changes only evaluation rows.
  full_ratio <- full_oos_pol$Deduct /
    (full_oos_pol$Coverage * cfg$coverage_to_dollars)
  if (any(!is.finite(full_ratio) | full_ratio <= 0))
    stop("Invalid deductible ratio in the full 2010 policy sample.")
  ratio_xlim_pct <- range(100 * full_ratio)
  ratio_bandwidth <- max(0.08, 0.20 * stats::IQR(log10(full_ratio)))
  
  oos_ratio_results <- oos_deductible_ratio_diagnostics(
    oos_pol = oos_pol,
    oos_predictions = oos_predictions,
    model_pairs = oos_main_pairs,
    oos_valid = oos_valid,
    oos_file = oos_file,
    exclude_max = exclude_max,
    coverage_to_dollars = cfg$coverage_to_dollars,
    ratio_xlim_pct = ratio_xlim_pct,
    bandwidth = ratio_bandwidth,
    min_ess = 15
  )
  
  # Additional diagnostic: continuous D/C metrics WITHIN each deductible amount.
  # Keeps the common full-2010 ratio axis and bandwidth across sensitivity runs.
  oos_ratio_by_level <- oos_ratio_by_deductible_level(
    oos_pol = oos_pol,
    full_oos_pol = full_oos_pol,
    oos_predictions = oos_predictions,
    model_pairs = oos_main_pairs,
    oos_valid = oos_valid,
    oos_file = oos_file,
    coverage_to_dollars = cfg$coverage_to_dollars,
    exclude_max = exclude_max,
    bandwidth = ratio_bandwidth,
    min_ess = 5L
  )
  
  # Table 1 and Table 2 are NOT reselected on the 2010 validation sample.
  # Copy the chosen TRAINING specifications solely as reference tables.
  if (exists("selected_severity")) write.csv(selected_severity,
                                             oos_report_file("Table1_TrainingSelectedSeverityCopulas_Reference.csv"), row.names = FALSE)
  if (exists("selected_frequency")) write.csv(selected_frequency,
                                              oos_report_file("Table2_TrainingSelectedFrequencyCopulas_Reference.csv"), row.names = FALSE)
  
  oos_all_scores <- do.call(rbind, lapply(seq_len(nrow(oos_main_pairs)), function(i) {
    pp <- oos_main_pairs[i, , drop = FALSE]
    cbind(pp[, c("Case", "Framework", "Scheme", "Model", "Frequency", "Severity"), drop = FALSE],
          oos_scores(oos_predictions[[pp$Model]]$Limited[oos_valid],
                     oos_pol$S_limited[oos_valid]))
  }))
  rownames(oos_all_scores) <- NULL
  oos_table3 <- oos_all_scores[, c(
    "Case", "Framework", "Scheme", "Model", "Frequency", "Severity", "N",
    "OP_ratio", "RelBias_pct", "Mean_Poisson_Deviance",
    "MAE_obs", "RMSE_obs", "Observed_total", "Predicted_total"
  )]
  write.csv(oos_table3,
            oos_report_file("Table3_2010_Limited_PurePremium_Main.csv"), row.names = FALSE)
  write.csv(oos_all_scores,
            oos_report_file("Appendix_2010_Limited_PurePremium_AllMetrics.csv"), row.names = FALSE)
  
  # No Case 3 / thinning-sensitivity ratemaking combinations are reported.
  
  # Additional output matching original raw (limited + unlimited) metrics.
  oos_metric_parts <- list(); oos_metric_k <- 0L
  oos_d_metric_parts <- list(); oos_d_metric_k <- 0L
  for (i in seq_len(nrow(model_pairs))) {
    pp <- model_pairs[i, , drop = FALSE]
    for (kind in c("Limited", "Unlimited")) {
      pred <- oos_predictions[[pp$Model]][[kind]]
      obs <- if (kind == "Limited") oos_pol$S_limited else oos_pol$S_raw
      oos_metric_k <- oos_metric_k + 1L
      oos_metric_parts[[oos_metric_k]] <- cbind(
        pp[, c("Case", "Framework", "Scheme", "Model", "Frequency", "Severity"), drop = FALSE],
        Payment_definition = kind, N_all = nrow(oos_pol),
        oos_scores(pred, obs), Sample = if (exclude_max) "Out_of_sample_2010_exclude_max" else "Out_of_sample_2010"
      )
      for (d in cfg$ded_levels) {
        ix <- which(oos_pol$Deduct == d)
        if (!length(ix)) next
        oos_d_metric_k <- oos_d_metric_k + 1L
        oos_d_metric_parts[[oos_d_metric_k]] <- cbind(
          pp[, c("Case", "Framework", "Scheme", "Model", "Frequency", "Severity"), drop = FALSE],
          Payment_definition = kind, Deductible = d, N_all = length(ix),
          oos_scores(pred[ix], obs[ix]), Sample = if (exclude_max) "Out_of_sample_2010_exclude_max" else "Out_of_sample_2010"
        )
      }
    }
  }
  oos_overall_metrics <- do.call(rbind, oos_metric_parts)
  oos_by_d_metrics <- do.call(rbind, oos_d_metric_parts)
  write.csv(oos_overall_metrics, oos_file("Premium_Performance_Overall_2010.csv"), row.names = FALSE)
  write.csv(oos_by_d_metrics, oos_file("Premium_Performance_ByDeductible_2010.csv"), row.names = FALSE)
  write.csv(oos_overall_metrics[oos_overall_metrics$Payment_definition == "Unlimited", , drop = FALSE],
            oos_report_file("Appendix_2010_Unlimited_Premiums.csv"), row.names = FALSE)
  
  # Keep the original primary/supplementary export naming convention.
  write.csv(oos_overall_metrics[oos_overall_metrics$Payment_definition == "Limited", , drop = FALSE],
            oos_file("Premium_Performance_CoverageLimited_Primary.csv"), row.names = FALSE)
  write.csv(oos_overall_metrics[oos_overall_metrics$Payment_definition == "Unlimited", , drop = FALSE],
            oos_file("Premium_Performance_Unlimited_Supplementary.csv"), row.names = FALSE)
  
  # Reproduce diagnostics for the THREE frameworks (Cases 1/4c/4d).
  oos_plot_deductible_metric <- function(metric, filename) {
    tt <- oos_by_d_metrics[
      oos_by_d_metrics$Payment_definition == "Limited", , drop = FALSE]
    z <- tt[[metric]]
    z <- z[is.finite(z)]
    if (!length(z)) return(invisible(NULL))
    use_ylog <- metric != "RelBias_pct" && all(z > 0)
    ylim <- if (metric == "RelBias_pct") {
      b <- max(1, abs(z)); c(-1.05*b, 1.05*b)
    } else if (use_ylog) {
      r <- range(z); c(r[1]/1.2, r[2]*1.2)
    } else {
      r <- range(z); r + c(-.08, .08) * max(diff(r), 1)
    }
    col <- c("firebrick3", "steelblue4", "darkgreen")
    grDevices::pdf(oos_file(filename), width = 11, height = 8)
    op <- par(no.readonly = TRUE)
    par(mfrow = c(3, 2), mar = c(4.2, 4.4, 2.6, .8))
    for (g in unique(oos_main_pairs$Framework)) {
      plot(NA_real_, NA_real_, xlim = range(cfg$ded_levels), ylim = ylim,
           log = if (use_ylog) "xy" else "x", xaxt = "n",
           xlab = "Deductible (USD)",
           ylab = if (metric == "RMSRE_positive") "RMSRE (observed S > 0)" else metric,
           main = g)
      axis(1, at = cfg$ded_levels,
           labels = format(cfg$ded_levels, scientific = FALSE, big.mark = ","),
           las = 2, cex.axis = .6)
      if (metric == "RelBias_pct") abline(h = 0, lty = 3, col = "grey50")
      for (j in 1:3) {
        sc <- paste0(if (g == "XR_marginal") "XR" else "ER", j)
        aa <- tt[tt$Framework == g & tt$Scheme == sc, , drop = FALSE]
        aa <- aa[order(aa$Deductible), , drop = FALSE]
        if (nrow(aa)) lines(aa$Deductible, aa[[metric]], type = "b",
                            col = col[j], lty = j, pch = 15+j)
      }
      legend("topleft", legend = paste0(if (g == "XR_marginal") "XR" else "ER", 1:3),
             col = col, lty = 1:3, pch = 16:18, bty = "n", cex = .72)
    }
    plot.new()
    text(.5, .65, "2010 OUT-OF-SAMPLE", cex = 1)
    text(.5, .4, "Training 2006-2009; fixed model parameters", cex = .92)
    par(op); dev.off()
    invisible(NULL)
  }
  oos_plot_deductible_metric("MAE_obs", "Premium_ByDeductible_MAE_obs_Limited.pdf")
  oos_plot_deductible_metric("RMSRE_positive", "Premium_ByDeductible_RMSRE_obs_Limited.pdf")
  oos_plot_deductible_metric("RelBias_pct", "Premium_ByDeductible_RB_obs_pct_Limited.pdf")
  
  # -------------------- DEDUCTIBLE-LEVEL O/P AND 95% BOOTSTRAP CI -------------
  oos_ded_parts <- list(); kk <- 0L
  for (i in seq_len(nrow(oos_main_pairs))) {
    pp <- oos_main_pairs[i, , drop = FALSE]
    p <- oos_predictions[[pp$Model]]$Limited
    for (d in cfg$ded_levels) {
      ix <- which(oos_valid & oos_pol$Deduct == d)
      if (!length(ix)) next
      kk <- kk + 1L
      oos_ded_parts[[kk]] <- cbind(
        pp[, c("Case", "Framework", "Scheme", "Model"), drop = FALSE],
        Deductible = d, oos_scores(p[ix], oos_pol$S_limited[ix])
      )
    }
  }
  if (!length(oos_ded_parts)) stop("No scored deductible groups in 2010.")
  oos_ded_scores <- do.call(rbind, oos_ded_parts)
  rownames(oos_ded_scores) <- NULL
  
  oos_cluster_ci <- function(y, p, ids, B, min_clusters) {
    sums <- rowsum(cbind(y = y, p = p), as.character(ids), reorder = FALSE)
    K <- nrow(sums)
    if (K < min_clusters || B < 2L)
      return(c(N_clusters = K, OP_lower = NA_real_, OP_upper = NA_real_))
    v <- rep(NA_real_, B)
    for (bb in seq_len(B)) {
      ix <- sample.int(K, K, replace = TRUE)
      v[bb] <- oos_op(sum(sums[ix, "y"]), sum(sums[ix, "p"]))
    }
    v <- v[is.finite(v)]
    if (length(v) < max(20L, ceiling(B / 2)))
      return(c(N_clusters = K, OP_lower = NA_real_, OP_upper = NA_real_))
    limits <- as.numeric(stats::quantile(v, c(.025, .975), names = FALSE))
    c(N_clusters = K, OP_lower = limits[1L], OP_upper = limits[2L])
  }
  
  set.seed(cfg$report_bootstrap_seed + 2010L)
  oos_ded_scores$N_clusters <- NA_integer_
  oos_ded_scores$OP_lower <- NA_real_
  oos_ded_scores$OP_upper <- NA_real_
  for (i in seq_len(nrow(oos_ded_scores))) {
    d <- oos_ded_scores$Deductible[i]
    nm <- oos_ded_scores$Model[i]
    ix <- which(oos_valid & oos_pol$Deduct == d)
    ci <- oos_cluster_ci(oos_pol$S_limited[ix],
                         oos_predictions[[nm]]$Limited[ix],
                         oos_pol$PolicyNum[ix],
                         B = as.integer(cfg$report_bootstrap_reps),
                         min_clusters = as.integer(cfg$report_min_bootstrap_clusters))
    oos_ded_scores$N_clusters[i] <- as.integer(ci["N_clusters"])
    oos_ded_scores$OP_lower[i] <- ci["OP_lower"]
    oos_ded_scores$OP_upper[i] <- ci["OP_upper"]
  }
  write.csv(oos_ded_scores,
            oos_report_file("Figure1_Deductible_OP_and_95pctClusterCI.csv"), row.names = FALSE)
  write.csv(oos_ded_scores,
            oos_report_file("Figure1B_Deductible_OP_by_ObservationScheme.csv"), row.names = FALSE)
  
  # Figure 1. Each framework gets one panel with the three schemes and CIs.
  oos_groups <- unique(oos_main_pairs$Framework)
  oos_schemes <- c(1L, 2L, 3L)
  oos_scheme_colors <- c("#21618C", "#D35400", "#148F77")
  oos_cap <- cfg$report_op_plot_cap
  
  oos_plot_framework <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mfrow = c(3, 2), mar = c(4.1, 4.2, 2.5, .8),
        oma = c(.5, 0, 1.4, 0))
    for (g in oos_groups) {
      tt <- oos_ded_scores[oos_ded_scores$Framework == g, , drop = FALSE]
      plot(NA_real_, NA_real_, xlim = range(cfg$ded_levels),
           ylim = c(0, oos_cap), log = "x", xaxt = "n",
           xlab = "Deductible (USD)", ylab = "Observed / predicted",
           main = g)
      axis(1, at = cfg$ded_levels,
           labels = format(cfg$ded_levels, scientific = FALSE, big.mark = ","),
           las = 2, cex.axis = .53)
      abline(h = 1, col = "grey45", lty = 2)
      ss <- unique(tt$Scheme)
      for (j in seq_along(ss)) {
        aa <- tt[tt$Scheme == ss[j], , drop = FALSE]
        aa <- aa[order(aa$Deductible), , drop = FALSE]
        xx <- aa$Deductible * (1 + (j - 2) * .025)
        yy <- pmin(aa$OP_ratio, oos_cap)
        yy[!is.finite(yy)] <- NA_real_
        lines(xx, yy, col = oos_scheme_colors[j], type = "b", lwd = 1.5,
              pch = 15 + j)
        ci <- is.finite(aa$OP_lower) & is.finite(aa$OP_upper)
        if (any(ci)) arrows(xx[ci], pmax(0, pmin(aa$OP_lower[ci], oos_cap)),
                            xx[ci], pmin(aa$OP_upper[ci], oos_cap),
                            angle = 90, code = 3, length = .025,
                            col = oos_scheme_colors[j])
      }
      legend("topright", legend = ss, col = oos_scheme_colors[seq_along(ss)],
             pch = 15 + seq_along(ss), lty = 1, bty = "n", cex = .75)
    }
    plot.new()
    text(.5, .65, "2010 OUT-OF-SAMPLE\n95% descriptive policy-cluster intervals", cex = 1)
    text(.5, .35, sprintf("Display capped at O/P = %s; CSV is uncapped", oos_cap), cex = .88)
    mtext("2010 validation | fitted on 2006-2009", outer = TRUE, side = 3, line = .1)
  }
  grDevices::pdf(oos_report_file("Figure1_Deductible_OP_95pctCI.pdf"), width = 11, height = 8.5)
  oos_plot_framework()
  grDevices::dev.off()
  
  # Figure 1B. Each observation scheme gets one plot of THREE frameworks.
  oos_framework_cols <- c("#B22222", "#228B22", "#4169E1")[seq_along(oos_groups)]
  stopifnot(length(oos_groups) == length(oos_framework_cols))
  oos_scheme_num <- function(z) as.integer(sub(".*([123])$", "\\1", z))
  oos_ded_scores$Scheme_number <- oos_scheme_num(oos_ded_scores$Scheme)
  if (anyNA(oos_ded_scores$Scheme_number) ||
      any(!oos_ded_scores$Scheme_number %in% 1:3))
    stop("Scheme labels must end in 1, 2 or 3.")
  oos_scheme_title <- c("Complete", "Censored", "Truncated")
  oos_draw_scheme <- function(s, show_legend = TRUE) {
    tt <- oos_ded_scores[oos_ded_scores$Scheme_number == s, , drop = FALSE]
    plot(NA_real_, NA_real_, xlim = range(cfg$ded_levels),
         ylim = c(0, oos_cap), log = "x", xaxt = "n",
         xlab = "Deductible (USD)", ylab = "Observed / predicted",
         main = paste0("Scheme ", s, " - ", oos_scheme_title[s]))
    axis(1, at = cfg$ded_levels,
         labels = format(cfg$ded_levels, scientific = FALSE, big.mark = ","),
         las = 2, cex.axis = .7)
    abline(h = 1, lty = 2, col = "grey45")
    for (j in seq_along(oos_groups)) {
      aa <- tt[tt$Framework == oos_groups[j], , drop = FALSE]
      aa <- aa[order(aa$Deductible), , drop = FALSE]
      yy <- pmin(aa$OP_ratio, oos_cap); yy[!is.finite(yy)] <- NA_real_
      lines(aa$Deductible, yy, type = "b", pch = 14 + j,
            col = oos_framework_cols[j], lty = j, lwd = 1.6)
    }
    if (show_legend) legend("topright", legend = oos_groups,
                            col = oos_framework_cols, lty = seq_along(oos_groups),
                            pch = 14 + seq_along(oos_groups), bty = "n", cex = .68)
  }
  for (s in 1:3) {
    grDevices::pdf(oos_report_file(sprintf("Figure1B_Scheme%d_%s_Deductible_OP.pdf",
                                           s, oos_scheme_title[s])), width = 10, height = 6.7)
    op <- par(no.readonly = TRUE)
    par(mar = c(6.4, 4.4, 3.2, 1.0))
    oos_draw_scheme(s)
    par(op); dev.off()
  }
  grDevices::pdf(oos_report_file("Figure1B_AllSchemes_Deductible_OP.pdf"),
                 width = 17, height = 6.5)
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 3), mar = c(6, 4.3, 3.1, .8), oma = c(.4, .2, 1.2, .2))
  for (s in 1:3) oos_draw_scheme(s, show_legend = TRUE)
  mtext("2010 OUT-OF-SAMPLE | Cases 1, 4c, 4d",
        outer = TRUE, side = 3, line = .1)
  par(op); dev.off()
  
  # -------------------- FIGURE 2: RISK-GROUP CALIBRATION -----------------------
  # Define cutpoints using ONLY the 2006-2009 training XR1 predictions.
  # Then apply those FIXED cutpoints to 2010; do not refit to 2010 outcomes.
  oos_base_name <- "XR_marginal__XR1"
  train_xr1 <- raw_predictions[[oos_base_name]]$Limited
  holdout_xr1 <- oos_predictions[[oos_base_name]]$Limited
  if (is.null(train_xr1) || is.null(holdout_xr1))
    stop("XR1 premiums unavailable for training-defined risk groups.")
  train_xr1 <- train_xr1[is.finite(train_xr1) & train_xr1 >= 0]
  if (!length(train_xr1)) stop("No finite training XR1 premiums for risk cutpoints.")
  G <- as.integer(cfg$report_risk_groups)
  if (G < 2L) stop("report_risk_groups must be >= 2")
  cut_raw <- as.numeric(stats::quantile(train_xr1,
                                        probs = seq(0, 1, length.out = G + 1L), names = FALSE, na.rm = TRUE))
  cut_raw <- unique(cut_raw)
  if (length(cut_raw) < 2L) {
    warning("All training XR1 premiums are equal; using a single 2010 risk group.")
    oos_breaks <- c(-Inf, Inf)
  } else {
    oos_breaks <- cut_raw
    oos_breaks[1L] <- -Inf
    oos_breaks[length(oos_breaks)] <- Inf
  }
  oos_risk_group <- as.integer(cut(holdout_xr1,
                                   breaks = oos_breaks, include.lowest = TRUE, labels = FALSE))
  write.csv(data.frame(Lower = head(oos_breaks, -1), Upper = tail(oos_breaks, -1)),
            oos_report_file("Figure2_TrainingXR1_RiskGroup_Breaks.csv"), row.names = FALSE)
  
  oos_risk_parts <- list(); kk <- 0L
  for (i in seq_len(nrow(oos_main_pairs))) {
    pp <- oos_main_pairs[i, , drop = FALSE]
    p <- oos_predictions[[pp$Model]]$Limited
    for (z in sort(unique(oos_risk_group[oos_valid]))) {
      ix <- which(oos_valid & oos_risk_group == z)
      if (!length(ix)) next
      kk <- kk + 1L
      oos_risk_parts[[kk]] <- cbind(
        pp[, c("Case", "Framework", "Scheme", "Model"), drop = FALSE],
        Risk_group = z, Baseline_mean_premium = mean(holdout_xr1[ix]),
        oos_scores(p[ix], oos_pol$S_limited[ix])
      )
    }
  }
  oos_risk_scores <- do.call(rbind, oos_risk_parts)
  rownames(oos_risk_scores) <- NULL
  write.csv(oos_risk_scores,
            oos_report_file("Figure2_BaselineRiskGroup_Calibration.csv"), row.names = FALSE)
  
  grDevices::pdf(oos_report_file("Figure2_Observed_vs_Predicted_ByRiskGroup.pdf"),
                 width = 11, height = 8.5)
  op <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2), mar = c(4.1, 4.6, 2.6, .9), oma = c(.9, 0, 1.1, 0))
  for (g in oos_groups) {
    gg <- oos_risk_scores[oos_risk_scores$Framework == g, , drop = FALSE]
    yr <- range(c(gg$Observed_mean, gg$Predicted_mean), finite = TRUE)
    ymax <- if (length(yr) && all(is.finite(yr))) max(1, yr[2L] * 1.12) else 1
    base <- gg[gg$Scheme == unique(gg$Scheme)[1L], , drop = FALSE]
    base <- base[order(base$Risk_group), , drop = FALSE]
    plot(base$Risk_group, base$Observed_mean, type = "b", pch = 19, lwd = 2,
         col = "#2C3E50", xlim = c(1, length(oos_breaks) - 1L),
         ylim = c(0, ymax), xlab = "Training-XR1 premium risk group",
         ylab = "Mean aggregate payment (USD)", main = g)
    ss <- unique(gg$Scheme)
    for (j in seq_along(ss)) {
      aa <- gg[gg$Scheme == ss[j], , drop = FALSE]
      aa <- aa[order(aa$Risk_group), , drop = FALSE]
      lines(aa$Risk_group, aa$Predicted_mean, type = "b", lwd = 1.5,
            col = oos_scheme_colors[j], pch = 15 + j, lty = j)
    }
    legend("topleft", legend = c("Observed", ss),
           col = c("#2C3E50", oos_scheme_colors[seq_along(ss)]),
           pch = c(19, 15 + seq_along(ss)),
           lty = c(1, seq_along(ss)), bty = "n", cex = .72)
  }
  plot.new()
  text(.5, .63, "2010 out-of-sample comparison", cex = 1)
  text(.5, .40, "Risk-group cutpoints are based on\n2006-2009 XR1 fitted premiums only", cex = .95)
  mtext("2010 validation | fixed training-defined risk bins", outer = TRUE, side = 3, line = .05)
  par(op); dev.off()
  
  # -------------------- APPENDIX: COMPONENT-WISE CONTRASTS ---------------------
  oos_contrasts <- do.call(rbind, lapply(c("4c", "4d"), function(ccase) {
    do.call(rbind, lapply(1:3, function(ss) {
      ba <- model_pairs[model_pairs$Case == "1" &
                          model_pairs$Scheme == paste0("XR", ss), , drop = FALSE]
      al <- model_pairs[model_pairs$Case == ccase &
                          model_pairs$Scheme == paste0("ER", ss), , drop = FALSE]
      data.frame(Case = ccase, Scheme = ss,
                 Base_Model = ba$Model, Alternative_Model = al$Model)
    }))
  }))
  oos_contrast_parts <- lapply(seq_len(nrow(oos_contrasts)), function(i) {
    cc <- oos_contrasts[i, , drop = FALSE]
    pa <- oos_predictions[[cc$Base_Model]]$Limited
    pb <- oos_predictions[[cc$Alternative_Model]]$Limited
    both <- oos_valid & is.finite(pa) & is.finite(pb) & pa >= 0 & pb >= 0
    aa <- oos_scores(pa[both], oos_pol$S_limited[both])
    bb <- oos_scores(pb[both], oos_pol$S_limited[both])
    cbind(cc, N = aa$N, Base_OP = aa$OP_ratio, Alternative_OP = bb$OP_ratio,
          Base_Poisson_Deviance = aa$Mean_Poisson_Deviance,
          Alternative_Poisson_Deviance = bb$Mean_Poisson_Deviance,
          Delta_Poisson_Deviance = bb$Mean_Poisson_Deviance - aa$Mean_Poisson_Deviance,
          Base_MAE = aa$MAE_obs, Alternative_MAE = bb$MAE_obs,
          Delta_MAE = bb$MAE_obs - aa$MAE_obs)
  })
  oos_contrast_scores <- do.call(rbind, oos_contrast_parts)
  write.csv(oos_contrast_scores,
            oos_report_file("Appendix_Case1_vs_Case4_Contrasts.csv"), row.names = FALSE)
  
  # -------------------- APPENDIX: COMMON COPULA SENSITIVITY --------------------
  # Reuses selection candidate fits FROM 2006-2009, does not fit candidates in 2010.
  # As in the original appendix, FEM1 is fixed to isolate severity dependence.
  oos_common_family_scores <- data.frame()
  if (isTRUE(cfg$report_common_copula_sensitivity)) {
    oos_key_to_model <- c(
      sev_cont_complete = "SCM1", sev_cont_censored = "SCM2",
      sev_cont_truncated = "SCM3", sev_disc_complete = "SDM1",
      sev_disc_censored = "SDM2", sev_disc_truncated = "SDM3")
    sensitivity_parts <- list(); kk <- 0L
    for (key in names(oos_key_to_model)) {
      mn <- unname(oos_key_to_model[key])
      chosen <- severity_fits[[mn]]
      candidates <- sel[[key]]$fits
      candidate_codes <- unique(c(chosen$spec$fam,
                                  as.integer(cfg$report_common_family_codes)))
      for (family_code in candidate_codes) {
        family <- fam_name(family_code)
        f <- candidates[[family]]
        available <- !is.null(f) && is.finite(f$AIC) &&
          !is.null(f$convergence) && f$convergence == 0L
        kk <- kk + 1L
        if (!available) {
          sensitivity_parts[[kk]] <- data.frame(
            Model = mn, Candidate_family = family,
            AIC_selected_family = fam_name(chosen$spec$fam),
            Used = FALSE, Reason = "Not fitted / nonconverged",
            oos_scores(numeric(0), numeric(0)))
        } else {
          p <- oos_freq_mean$FEM1 *
            conditional_paid_limited(f, oos_Dtheta(f), oos_dat)
          p[oos_zero_limit] <- 0
          ix <- oos_valid & is.finite(p) & p >= 0
          sensitivity_parts[[kk]] <- data.frame(
            Model = mn, Candidate_family = family,
            AIC_selected_family = fam_name(chosen$spec$fam),
            Used = TRUE, Reason = "FEM1 fixed; training candidate only",
            oos_scores(p[ix], oos_pol$S_limited[ix]))
        }
      }
    }
    oos_common_family_scores <- do.call(rbind, sensitivity_parts)
    write.csv(oos_common_family_scores,
              oos_report_file("Appendix_CommonCopula_SeverityOnly.csv"), row.names = FALSE)
  }
  
  # -------------------- APPENDIX: NUMERICAL AND GB2 DIAGNOSTICS --------------
  if (exists("moments"))
    write.csv(moments,
              oos_report_file("Appendix_GB2_FirstMoment_Existence.csv"), row.names = FALSE)
  if (isTRUE(cfg$capped_quad_check)) {
    id <- split(seq_len(nrow(oos_pol)), oos_pol$Deduct)
    ix <- sort(unique(unlist(lapply(id, function(z)
      z[unique(round(seq(1, length(z), length.out = min(5L, length(z)))))]),
      use.names = FALSE)))
    oos_dat_small <- lapply(oos_dat, function(v)
      if (is.matrix(v)) v[ix, , drop = FALSE] else v[ix])
    gl_lo <- gl_nodes(cfg$premium_gl_n)
    gl_hi <- gl_nodes(2L * cfg$premium_gl_n)
    oos_quadrature <- do.call(rbind, lapply(severity_fits, function(f) {
      a <- conditional_paid_limited(f, oos_Dtheta(f), oos_dat_small, gl_lo)
      b <- conditional_paid_limited(f, oos_Dtheta(f), oos_dat_small, gl_hi)
      rel <- abs(a - b) / pmax(1, abs(b))
      data.frame(Model = f$name, N_checked = length(ix),
                 Max_relative_change = max(rel, na.rm = TRUE),
                 Stable = all(is.finite(rel) & rel <= cfg$capped_quad_tol))
    }))
    write.csv(oos_quadrature,
              oos_report_file("Appendix_2010_CoverageLimited_Quadrature.csv"), row.names = FALSE)
    if (any(!oos_quadrature$Stable)) warning("Some 2010 limited premium quadratures are unstable. See appendix CSV.")
  }
  
  # -------------------- IN-SAMPLE vs OUT-OF-SAMPLE TABLE -----------------------
  if (exists("manuscript_report") && !is.null(manuscript_report$table3)) {
    tr <- manuscript_report$table3
    fields <- c("N", "OP_ratio", "RelBias_pct", "Mean_Poisson_Deviance", "MAE_obs", "RMSE_obs")
    tr <- tr[, c("Model", fields), drop = FALSE]
    names(tr)[-1L] <- paste0("Train_2006_2009_", fields)
    te <- oos_table3[, c("Model", fields), drop = FALSE]
    names(te)[-1L] <- paste0("Test_2010_", fields)
    oos_train_test <- merge(tr, te, by = "Model", sort = FALSE)
    oos_train_test <- oos_train_test[match(oos_main_pairs$Model, oos_train_test$Model), ]
    write.csv(oos_train_test, oos_report_file("Table_Train_vs_2010.csv"), row.names = FALSE)
  }
  
  # -------------------- 2010 DATA-AUDIT AND RESULTS ARCHIVE -------------------
  oos_audit <- data.frame(
    Training_years = paste(sort(unique(pol$Year)), collapse = ","),
    Test_year = oos_year,
    Train_policy_years = nrow(pol),
    Test_policy_years = nrow(oos_pol),
    Test_claim_records = if (exclude_max)
      sum(oos2010_claims$PolicyNum != excluded_policy$PolicyNum)
    else nrow(oos2010_claims),
    Excluded_policy_years = as.integer(exclude_max),
    Excluded_PolicyNum = if (exclude_max) as.character(excluded_policy$PolicyNum) else NA_character_,
    Full_2010_policy_years = nrow(full_oos_pol),
    Complete_case_test_policy_years = n_oos,
    Test_zero_limited_payment_policy_years = sum(oos_pol$S_limited == 0),
    Test_zero_coverage_above_deductible = sum(oos_zero_limit),
    Coverage_conversion_to_dollars = cfg$coverage_to_dollars,
    Modeling_coverage_interpretation = "Per-occurrence ground-up limit (assumed)",
    Copulas_reselected_on_2010 = FALSE,
    Models_refitted_on_2010 = FALSE,
    stringsAsFactors = FALSE
  )
  write.csv(oos_audit, oos_file("OutOfSample_2010_Data_Audit.csv"), row.names = FALSE)
  write.csv(model_pairs, oos_file("Premium_Model_Pairs_TrainingFixed.csv"), row.names = FALSE)
  
  notes2010 <- c(
    "OUT-OF-SAMPLE EMPIRICAL VALIDATION: YEAR 2010",
    "Training: policy years 2006-2009; validation: policy year 2010.",
    if (exclude_max) "POST-HOC SENSITIVITY: largest 2010 paid-loss policy-year excluded from evaluation only."
    else "PRIMARY OOS EVALUATION: all 2010 policy-years retained.",
    "All model and copula parameters are frozen at 2006-2009 estimates.",
    "No selection, calibration or parameter fitting occurs on the 2010 outcomes.",
    paste0("Common complete-case evaluation n=", n_oos,
           " of ", nrow(oos_pol), " policy-years (all nine headline models)."),
    "Primary loss: per-occurrence (min(Y, C_i)-D_i)_+ summed within policy-year.",
    paste0("C_i = Coverage_i * ", cfg$coverage_to_dollars,
           " dollars; contract interpretation is assumed, not verified."),
    "O/P = sum(observed losses) / sum(predicted pure premiums).",
    "Poisson deviance, MAE, RMSE and relative portfolio bias use zero-loss policies.",
    "RMSRE_positive excludes zero observed limited losses; N_positive is reported.",
    "O/P = NA for observed=predicted=0; O/P = Inf for observed>0 and predicted=0.",
    "Figure 1: policy-cluster bootstrap CIs condition on fixed fitted predictions.",
    "Figure 1B: cases 1, 4c, 4d in each observation scheme.",
    "Figure 2: test groups use only TRAINING-period XR1 risk cutpoints.",
    "Tables 1-2 are fixed training copula specifications: NOT test-year AIC fits.",
    "Appendix common copulas uses cached training candidates and FEM1 fixed.",
    "If some test predictions are invalid, a common subset is used for all nine models.",
    "Out-of-sample does not mean causally identified; frequency-severity factorization remains assumed."
  )
  writeLines(notes2010, oos_report_file("README_OutOfSample_2010.txt"))
  
  oos_2010_results <- list(Excluded_max_policy = exclude_max,
                           Excluded_policy = excluded_policy,
                           train_years = 2006:2009, test_year = 2010L,
                           audit = oos_audit, valid_common = which(oos_valid),
                           model_pairs = model_pairs, table3 = oos_table3,
                           appendix_metrics = oos_all_scores,
                           deductible_op = oos_ded_scores,
                           risk_group_calibration = oos_risk_scores,
                           component_contrasts = oos_contrast_scores,
                           common_family_sensitivity = oos_common_family_scores,
                           predictions = oos_predictions,
                           individual_policies = oos_pol
  )
  saveRDS(oos_2010_results, oos_file("OutOfSample_2010_Results.rds"))
  cat("\n2010 out-of-sample Table 3 (all nine headline models):\n")
  print(oos_table3[, c("Framework", "Scheme", "N", "OP_ratio", "RelBias_pct",
                       "Mean_Poisson_Deviance", "MAE_obs", "RMSE_obs")],
        row.names = FALSE, digits = 4)
  cat("\n2010 evaluation complete. Results saved under: ",
      normalizePath(oos_dir), "\n", sep = "")
  
  # The two OOS branches share the same frozen 2006-2009 model estimates.
  invisible(list(table3 = oos_table3, full_policy_data = full_oos_pol,
                 freq_mean = oos_freq_mean, sev_limited = oos_sev_limited,
                 sev_unlimited = oos_sev_unlimited,
                 excluded_policy = excluded_policy, directory = oos_dir))
}

if (isTRUE(cfg$run_oos_full) || isTRUE(cfg$run_oos_exclude_max)) {
  if (!isTRUE(cfg$run_oos_full) && isTRUE(cfg$run_oos_exclude_max))
    warning("Exclusion requested without full OOS; full analysis is recommended")
  if (isTRUE(cfg$run_oos_full))
    oos_full_result <- run_oos_2010(FALSE)
  if (isTRUE(cfg$run_oos_exclude_max))
    oos_excluded_result <- run_oos_2010(TRUE, reuse_full = if (isTRUE(cfg$run_oos_full)) oos_full_result else NULL)
  if (isTRUE(cfg$run_oos_full) && isTRUE(cfg$run_oos_exclude_max)) {
    aa <- oos_full_result$table3
    bb <- oos_excluded_result$table3
    keys <- c("Case", "Scheme", "Model")
    keep <- c(keys, "N", "OP_ratio", "RelBias_pct",
              "Mean_Poisson_Deviance", "MAE_obs", "RMSE_obs")
    comparison <- merge(aa[, keep], bb[, keep], by = keys,
                        suffixes = c("_Full", "_ExcludeMax"), sort = FALSE)
    comparison <- comparison[match(model_pairs$Model, comparison$Model), , drop = FALSE]
    write.csv(comparison,
              file.path(cfg$output_dir, "Table_2010_Full_vs_ExcludeMax.csv"), row.names = FALSE)
  }
  cat("\nBoth 2010 OOS configurations completed with frozen 2006-2009 fits.\n")
}