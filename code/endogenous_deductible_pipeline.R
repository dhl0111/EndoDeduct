################################################################################
##  Insurance Ratemaking with Endogenous Deductibles  --  estimation pipeline
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
  out_file      = "endogenous_deductible_results.RData"
)
if (exists("cfg_override")) cfg[names(cfg_override)] <- cfg_override   # optional overrides set before source()

## ---- 0b. data ----------------------------------------------------------------
load("poldat.RData")          # policy-year data:   PolicyNum, Year, Type*, Coverage, Deduct, Freq
load("claimdat.RData")        # claim data:         PolicyNum, Year, Type*, Coverage, Deduct, LossBeforeDeductible

pol <- poldat[poldat$Year %in% cfg$years_in, ]
clm <- claimdat[claimdat$Year %in% cfg$years_in, ]

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
## 4. marginal models
################################################################################
cat("\n==== 4. marginal models ====\n")
## 4a. deductible models ---------------------------------------------------------
st_Dc <- stage_D("cont")
init_Dc <- c(coef(lm(log(dedratio) ~ TypeCity + TypeCounty + TypeSchool + TypeTown + TypeVillage, data = pol)), 1, 1, 1)
fit_Dc  <- fit_stage(st_Dc, init_Dc)

st_Dd <- stage_D("disc")
po    <- polr(factor(discDed, levels = 1:K_ded) ~ TypeCity + TypeCounty + TypeSchool + TypeTown + TypeVillage + log(Coverage), data = pol, Hess = FALSE)
fit_Dd <- fit_stage(st_Dd, c(po$zeta, -coef(po)))    # P(D<=j) = logistic(zeta_j - x'beta)  =>  a = zeta, b = -beta

## 4b. marginal severity models (GB2) -----------------------------------------------
init_sev <- function(dat) c(coef(lm(log(dat$ystar) ~ dat$X - 1)), 1, 1, 1)
SEM1 <- fit_stage(stage_sev("SEM1", "complete",  "none"), init_sev(dat_clm))
SEM2 <- fit_stage(stage_sev("SEM2", "censored",  "none"), init_sev(dat_pos))
SEM3 <- fit_stage(stage_sev("SEM3", "truncated", "none"), init_sev(dat_pos))

## 4c. marginal frequency models (NB) ----------------------------------------------
init_freq <- c(coef(glm(pol$Freq ~ X_pol - 1, family = poisson)), 1)
sev_none  <- list(dedtype = "none", fam = 0)
FEM1  <- fit_stage(stage_freq("FEM1",  "censored",  "none"), init_freq)
FEM2  <- fit_stage(stage_freq("FEM2",  "truncated", "none", sev_spec = sev_none), init_freq, prev = list(S = SEM3$par))   # p from SEM3
FEM2G <- fit_stage(stage_freq("FEM2G", "truncated", "none", sev_spec = sev_none), init_freq, prev = list(S = SEM1$par))   # p from SEM1

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

## 4e. sandwich s.e. of the marginal models ------------------------------------------
sw <- list(); se_compare <- NULL
##  SEM1-3 and FEM1 are single-stage (sandwich = robust stage-wise s.e.);
##  FEM2 / FEM2G depend on the severity stage through p(d;x)  ->  H_32 != 0
if (cfg$sandwich) {
  cat("\n==== 4e. sandwich s.e. for the marginal models ====\n")
  sw <- list()
  for (nm in c("SEM1", "SEM2", "SEM3", "FEM1")) sw[[nm]] <- sandwich_se(list(get(nm)), verbose = FALSE)
  sw$FEM2  <- sandwich_se(list(SEM3, FEM2))
  sw$FEM2G <- sandwich_se(list(SEM1, FEM2G))
  sw$D_cont <- sandwich_se(list(fit_Dc), verbose = FALSE)
  sw$D_disc <- sandwich_se(list(fit_Dd), verbose = FALSE)
  tex_table(est_se_matrix(last_stage_table(sw$SEM1), last_stage_table(sw$SEM2), last_stage_table(sw$SEM3)),
            "Regression estimates of the marginal severity models (sandwich s.e.)", header = c(" " = 1, "SEM1" = 2, "SEM2" = 2, "SEM3" = 2))
  tex_table(est_se_matrix(last_stage_table(sw$FEM1), last_stage_table(sw$FEM2G), last_stage_table(sw$FEM2)),
            "Regression estimates of the marginal frequency models (sandwich s.e.)", header = c(" " = 1, "FEM1" = 2, "FEM2G" = 2, "FEM2" = 2))
  tex_table(est_se_matrix(last_stage_table(sw$D_cont)), "Regression estimates of the continuous deductible model")
  tex_table(est_se_matrix(last_stage_table(sw$D_disc)), "Regression estimates of the discrete deductible model")
}

################################################################################
## 5. copula selection
################################################################################
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

## checkpointing: every fitted copula is written to cfg$cache_file[1] immediately, so an interrupted selection resumes
.ckpt <- new.env(); .ckpt$sel <- list(); .ckpt$partial <- list()
ckpt_load <- function() if (file.exists(cfg$cache_file[1])) { e <- new.env(); load(cfg$cache_file[1], envir = e)
  .ckpt$sel <- if (!is.null(e$sel)) e$sel else list(); .ckpt$partial <- if (!is.null(e$partial)) e$partial else list() }
ckpt_save <- function() { sel <- .ckpt$sel; partial <- .ckpt$partial; save(sel, partial, file = cfg$cache_file[1]) }

select_copula <- function(component, scenario, dedtype, init_marg, tau0, prev, sev_spec = NULL, families = cfg$sel_families, key = NULL) {
  fams <- if (is.null(families)) cop_families$fam else families
  done <- if (!is.null(key)) .ckpt$partial[[key]] else NULL
  fits <- lapply(fams, function(fam) {
    if (!is.null(done[[fam_name(fam)]])) { cat("   (checkpoint)", component, scenario, dedtype, fam_name(fam), "\n"); return(done[[fam_name(fam)]]) }
    cb <- cop_bounds(fam)
    init <- c(init_marg, min(max(tau0, cb$lo[1] + 1e-3), cb$hi[1] - 1e-3), if (fam == 2) 5)
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
best_fit <- function(selres, name) {           # AIC-best copula for one joint model; returns the fitted stage renamed
  i <- which.min(selres$table[, "AIC"]); f <- selres$fits[[i]]; f$name <- name; f
}

## -- the whole selection for one deductible type: severity (3 scenarios), then frequency (2 scenarios, the truncated
##    one using p(d;x) from the AIC-best truncated severity model). Independent across deductible types, so the two
##    can be run in parallel (see README).
run_copula_selection <- function(dt) {
  fitD <- if (dt == "cont") fit_Dc else fit_Dd
  do_key <- function(key, ...) { if (is.null(.ckpt$sel[[key]])) { .ckpt$sel[[key]] <- select_copula(..., key = key); ckpt_save() }; .ckpt$sel[[key]] }
  for (sc in c("complete", "censored", "truncated")) {
    marg <- switch(sc, complete = SEM1, censored = SEM2, truncated = SEM3)
    dat  <- if (sc == "truncated") dat_pos else dat_clm
    tau0 <- cor(F_D_of(dt, dat, fitD$par), marg_cdf_sev(marg, dat, sc), method = "kendall")
    do_key(paste("sev", dt, sc, sep = "_"), "sev", sc, dt, marg$par, tau0, prev = list(D = fitD$par))
  }
  S3 <- best_fit(.ckpt$sel[[paste("sev", dt, "truncated", sep = "_")]], "S3")
  ss <- list(dedtype = dt, fam = S3$spec$fam)
  tau0 <- cor(F_D_of(dt, dat_pol, fitD$par), marg_cdf_freq(FEM1, "censored", NULL), method = "kendall")
  do_key(paste("freq", dt, "censored", sep = "_"), "freq", "censored", dt, FEM1$par, tau0, prev = list(D = fitD$par))
  tau0 <- cor(F_D_of(dt, dat_pol, fitD$par), marg_cdf_freq(FEM2, "truncated", SEM3$par), method = "kendall")
  do_key(paste("freq", dt, "truncated", sep = "_"), "freq", "truncated", dt, FEM2$par, tau0, prev = list(D = fitD$par, S = S3$par), sev_spec = ss)
  invisible(NULL)
}

if (cfg$run_selection) {
  cat("\n==== 5. copula selection ====\n")
  ckpt_load()                                     # resume from an interrupted run, if any
  for (dt in cfg$sel_dedtypes) run_copula_selection(dt)
  sel <- .ckpt$sel
} else {                                           # load (and merge) cached selections
  sel <- list(); for (f in cfg$cache_file) { e <- new.env(); load(f, envir = e); sel <- c(sel, e$sel) }
}
if (cfg$stop_after_selection) quit(save = "no")

################################################################################
## 6. final joint models: the AIC-best copula of each joint model
##    (cfg$auto_select = FALSE: the families fixed in cfg$paper_families instead)
################################################################################
cat("\n==== 6. final joint models ====\n")
add_naive_hessian <- function(f) {                 # selection fits were made without the stage-wise Hessian
  nll <- function(th) -sum(f$contrib(th, f$prev))
  f$hess     <- tryCatch(numDeriv::hessian(nll, f$par, method.args = list(d = cfg$hess_d, r = cfg$hess_r)), error = function(e) NULL)
  f$se_naive <- tryCatch(sqrt(pmax(diag(solve(f$hess)), 0)), error = function(e) rep(NA_real_, length(f$par)))
  f
}
final_model <- function(name, key) {
  if (cfg$auto_select) f <- best_fit(sel[[key]], name)
  else { f <- sel[[key]]$fits[[fam_name(cfg$paper_families[[name]])]]; f$name <- name }
  cat(sprintf("   %-5s <- %-16s (loglik %.2f, AIC %.2f)\n", name, fam_name(f$spec$fam), f$loglik, f$AIC))
  add_naive_hessian(f)
}
SCM1 <- final_model("SCM1", "sev_cont_complete"); SCM2 <- final_model("SCM2", "sev_cont_censored"); SCM3 <- final_model("SCM3", "sev_cont_truncated")
SDM1 <- final_model("SDM1", "sev_disc_complete"); SDM2 <- final_model("SDM2", "sev_disc_censored"); SDM3 <- final_model("SDM3", "sev_disc_truncated")
FCM1 <- final_model("FCM1", "freq_cont_censored"); FCM2 <- final_model("FCM2", "freq_cont_truncated")
FDM1 <- final_model("FDM1", "freq_disc_censored"); FDM2 <- final_model("FDM2", "freq_disc_truncated")
chosen <- sapply(list(SCM1 = SCM1, SCM2 = SCM2, SCM3 = SCM3, SDM1 = SDM1, SDM2 = SDM2, SDM3 = SDM3,
                      FCM1 = FCM1, FCM2 = FCM2, FDM1 = FDM1, FDM2 = FDM2), function(f) f$spec$fam)
stopifnot(identical(FCM2$spec$sev_spec$fam, SCM3$spec$fam), identical(FDM2$spec$sev_spec$fam, SDM3$spec$fam))

## copula selection tables (Appendix)
for (nm in names(sel)) tex_table(sel[[nm]]$table, paste("Goodness-of-fit statistics for the joint model:", nm), digits = 2)

## 6b. sandwich standard errors of the joint models ----------------------------------
##   chain = the stages whose estimates enter the model, in estimation order
if (cfg$sandwich) {
  cat("\n==== 6b. sandwich s.e. for the joint models ====\n")
  chains <- list(SCM1 = list(fit_Dc, SCM1), SCM2 = list(fit_Dc, SCM2), SCM3 = list(fit_Dc, SCM3),
                 SDM1 = list(fit_Dd, SDM1), SDM2 = list(fit_Dd, SDM2), SDM3 = list(fit_Dd, SDM3),
                 FCM1 = list(fit_Dc, FCM1), FCM2 = list(fit_Dc, SCM3, FCM2),
                 FDM1 = list(fit_Dd, FDM1), FDM2 = list(fit_Dd, SDM3, FDM2))
  if (!is.null(cfg$sw_preload)) for (f in cfg$sw_preload) { e <- new.env(); load(f, envir = e); sw <- c(sw, e$sw[setdiff(names(e$sw), names(sw))]) }  # reuse precomputed chains
  for (nm in names(chains)) { if (!is.null(sw[[nm]])) { cat(nm, "(preloaded)\n"); next }; cat(nm, "\n"); sw[[nm]] <- sandwich_se(chains[[nm]]); save(sw, file = "sw_checkpoint.RData") }

  tex_table(est_se_matrix(last_stage_table(sw$SCM1), last_stage_table(sw$SCM2), last_stage_table(sw$SCM3)),
            "Joint severity models with a continuous deductible (sandwich s.e.)",
            header = c(" " = 1, "SCM1" = 2, "SCM2" = 2, "SCM3" = 2))
  tex_table(est_se_matrix(last_stage_table(sw$SDM1), last_stage_table(sw$SDM2), last_stage_table(sw$SDM3)),
            "Joint severity models with a discrete deductible (sandwich s.e.)",
            header = c(" " = 1, "SDM1" = 2, "SDM2" = 2, "SDM3" = 2))
  tex_table(est_se_matrix(last_stage_table(sw$FCM1), last_stage_table(sw$FCM2)),
            "Joint frequency models with a continuous deductible (sandwich s.e.)",
            header = c(" " = 1, "FCM1" = 2, "FCM2" = 2))
  tex_table(est_se_matrix(last_stage_table(sw$FDM1), last_stage_table(sw$FDM2)),
            "Joint frequency models with a discrete deductible (sandwich s.e.)",
            header = c(" " = 1, "FDM1" = 2, "FDM2" = 2))

  ## comparison of the three s.e. flavours (for the response letter)
  se_compare <- do.call(rbind, lapply(names(sw), function(nm) cbind(chain = nm, last_stage_table(sw[[nm]]))))
  print(se_compare[, c("chain", "parameter", "estimate", "se_naive", "se_robust", "se_sandwich", "ratio")], digits = 4, row.names = FALSE)
}

################################################################################
## 7. goodness of fit
################################################################################
cat("\n==== 7. goodness of fit ====\n")
## 7a. expected frequency distributions -----------------------------------------------
##  freq_pmf_matrix(): m x (nmax+1) matrix of P(N = j | d_i, x_i), j = 0..nmax, reusing the likelihood code
freq_pmf_matrix <- function(stage_fit, prev, scenario = stage_fit$spec$scenario, sev_spec = stage_fit$spec$sev_spec, nmax = cfg$nmax_gof) {
  sp <- stage_fit$spec
  sapply(0:nmax, function(j) {
    d <- dat_pol; d$n_cen <- rep(j, nrow(X_pol)); d$n_tru <- rep(j, nrow(X_pol))
    exp(freq_loglik_i(stage_fit$par, prev, scenario, sp$dedtype, sp$fam, sev_spec, d))
  })
}
gof_column <- function(pm, empirical) {           # expected counts 0..nmax-1 and >= nmax, plus chi-square
  E <- c(colSums(pm)[1:cfg$nmax_gof], nrow(pm) - sum(colSums(pm)[1:cfg$nmax_gof]))
  c(E, chisq = sum((E - empirical)^2 / E))
}
emp_counts <- function(n) { e <- tabulate(n + 1, nbins = cfg$nmax_gof)[1:cfg$nmax_gof]; c(e, length(n) - sum(e), chisq = 0) }
emp_cen <- emp_counts(dat_pol$n_cen); emp_tru <- emp_counts(dat_pol$n_tru)

## "model 1 -> truncated N" columns: the fitted N* model thinned with p(d;x) from the matching complete-data severity model
gof_freq <- cbind(
  Empirical = emp_cen,
  FEM1 = gof_column(freq_pmf_matrix(FEM1, list()), emp_cen[1:(cfg$nmax_gof + 1)]),
  FCM1 = gof_column(freq_pmf_matrix(FCM1, list(D = fit_Dc$par)), emp_cen[1:(cfg$nmax_gof + 1)]),
  FDM1 = gof_column(freq_pmf_matrix(FDM1, list(D = fit_Dd$par)), emp_cen[1:(cfg$nmax_gof + 1)]),
  Empirical = emp_tru,
  FEM1 = gof_column(freq_pmf_matrix(FEM1, list(S = SEM1$par), "truncated", sev_none), emp_tru[1:(cfg$nmax_gof + 1)]),
  FEM2 = gof_column(freq_pmf_matrix(FEM2, list(S = SEM3$par)), emp_tru[1:(cfg$nmax_gof + 1)]),
  FCM1 = gof_column(freq_pmf_matrix(FCM1, list(D = fit_Dc$par, S = SCM1$par), "truncated", list(dedtype = "cont", fam = SCM1$spec$fam)), emp_tru[1:(cfg$nmax_gof + 1)]),
  FCM2 = gof_column(freq_pmf_matrix(FCM2, list(D = fit_Dc$par, S = SCM3$par)), emp_tru[1:(cfg$nmax_gof + 1)]),
  FDM1 = gof_column(freq_pmf_matrix(FDM1, list(D = fit_Dd$par, S = SDM1$par), "truncated", list(dedtype = "disc", fam = SDM1$spec$fam)), emp_tru[1:(cfg$nmax_gof + 1)]),
  FDM2 = gof_column(freq_pmf_matrix(FDM2, list(D = fit_Dd$par, S = SDM3$par)), emp_tru[1:(cfg$nmax_gof + 1)]))
rownames(gof_freq) <- c(0:(cfg$nmax_gof - 1), paste0("$>=", cfg$nmax_gof, "$"), "$\\chi^2$--stat")
tex_table(gof_freq[, c(1:2, 5:7)], "Goodness-of-fit of the marginal frequency models", digits = 2, header = c(" " = 1, "Censored" = 2, "Truncated" = 3))
tex_table(gof_freq, "Goodness-of-fit of the joint frequency models", digits = 2, header = c(" " = 1, "Censored" = 4, "Truncated" = 7))

## 7b. QQ-plots of the ground-up losses: F_{Y*|D}(y*|d) under each severity model --------
Fcond_Y <- function(fit, theta_D = NULL) {
  sp <- fit$spec; k <- ncol(X_clm); th <- fit$par
  mu <- drop(X_clm %*% th[1:k]); FY <- pGB2(dat_clm$ystar, mu, th[k + 1], th[k + 2], th[k + 3])
  if (sp$dedtype == "none") return(FY)
  par <- tau2par(sp$fam, th[k + 4]); par2 <- if (sp$fam == 2) th[k + 5] else 0
  if (sp$dedtype == "cont") return(cop("h1", D_cont_cdf(theta_D, XD(X_clm, "cont"), dat_clm$dratio), FY, sp$fam, par, par2))
  up <- D_disc_cdf(theta_D, XD(X_clm, "disc"), dat_clm$dlev); um <- D_disc_cdf(theta_D, XD(X_clm, "disc"), dat_clm$dlev - 1)
  cop_rect(um, up, FY, sp$fam, par, par2) / (up - um)     # eq. (6): discrete conditional cdf F_{Y*|D}(y|d)
}
qq_list <- list(SEM1 = Fcond_Y(SEM1), SEM2 = Fcond_Y(SEM2), SEM3 = Fcond_Y(SEM3),
                SCM1 = Fcond_Y(SCM1, fit_Dc$par), SCM2 = Fcond_Y(SCM2, fit_Dc$par), SCM3 = Fcond_Y(SCM3, fit_Dc$par),
                SDM1 = Fcond_Y(SDM1, fit_Dd$par), SDM2 = Fcond_Y(SDM2, fit_Dd$par), SDM3 = Fcond_Y(SDM3, fit_Dd$par))
pdf("qq_groundup_losses.pdf", width = 10, height = 10)
par(mfrow = c(3, 3), pty = "s", mar = c(4, 4, 2, 1))
for (nm in names(qq_list)) { qqnorm(qnorm(clamp01(qq_list[[nm]])), col = "blue", main = nm, xlim = c(-8, 6), ylim = c(-8, 6)); abline(0, 1, col = "red") }
dev.off()

## 7c. variation of the linear predictor mu with log(Coverage) by entity type (marginal severity models) ----
pdf("mu_vs_logcoverage.pdf", width = 10, height = 6)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 2), mgp = c(2, 1, 0))
types <- c("TypeMisc", "TypeCity", "TypeCounty", "TypeSchool", "TypeTown", "TypeVillage")
type_of <- ifelse(rowSums(X_clm[, 2:6]) == 0, "TypeMisc", types[-1][max.col(X_clm[, 2:6], ties.method = "first")])
mus <- sapply(list(SEM1, SEM2, SEM3), function(f) drop(X_clm %*% f$par[1:ncol(X_clm)]))
for (tp in types) {
  i <- which(type_of == tp); if (length(i) < 2) next
  o <- order(X_clm[i, "log(Coverage)"]); lc <- X_clm[i, "log(Coverage)"][o]
  matplot(lc, mus[i, , drop = FALSE][o, , drop = FALSE], type = "l", lty = 1, col = c("red", "green3", "blue"), xlab = "log(Coverage)", ylab = expression(mu), main = tp)
  legend("top", legend = c("SEM1", "SEM2", "SEM3"), col = c("red", "green3", "blue"), lty = 1, horiz = TRUE, bty = "n", cex = 0.8)
}
dev.off()

## ---- save everything ---------------------------------------------------------------
save.image(cfg$out_file)      # everything: data, functions, fits, `sel`, `sw`, `se_compare`, `gof_freq`, `sev_family_tab`, `qq_list`
cat("\nDone. Results saved to", cfg$out_file, "\n")
