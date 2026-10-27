################################################################################
##  Summaries of the discrete-deductible simulation: Mean / s.e. / RMSE / Bias tables (Appendix layout of the
##  paper), MAE tables, coverage of the naive vs Godambe s.e., box plots.
##  usage:  Rscript summarize_sim.R [results_dir] [tex_dir]
################################################################################
args <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
texdir <- if (length(args) >= 2) args[2] else "tex_sim"
dir.create(texdir, showWarnings = FALSE)
pipeline_file <- "endogenous_deductible_pipeline.R"; source("sim_functions.R")

files <- list.files(resdir, "^case[0-9]+_rep[0-9]+\\.rds$", full.names = TRUE)
res <- lapply(files, readRDS)
cases <- sort(unique(sapply(res, `[[`, "case")))
R <- sapply(cases, function(cs) sum(sapply(res, `[[`, "case") == cs)); names(R) <- cases
cat("replications per case:\n"); print(R)

## ---- true values ------------------------------------------------------------------------------
truth <- function(model, case) {
  rho <- design$cases[[as.character(case)]]
  switch(substr(model, 1, 1),
         D = c(design$kappa, design$beta_D),
         S = c(design$beta_Y, design$alpha, if (grepl("DM", model)) rho[["rho_Y"]]),
         F = c(design$beta_N, design$r,     if (grepl("DM", model)) rho[["rho_N"]]))
}
par_mat <- function(model, case) t(sapply(res[sapply(res, `[[`, "case") == case], function(r) r$par[[model]]))
stats <- function(model, case) {
  P <- par_mat(model, case); th <- truth(model, case); Rc <- nrow(P)
  m <- colMeans(P); v <- colMeans(sweep(P, 2, m)^2)
  data.frame(parameter = colnames(P), truth = th, mean = m, se = sqrt(v / Rc), rmse = sqrt(colMeans(sweep(P, 2, th)^2)),
             bias = m - th, var = v, row.names = NULL)
}
f4 <- function(x) formatC(x, format = "f", digits = 4)

## ---- Appendix-style tables -----------------------------------------------------------------------
## groups: list(label = <latex label>, idx = <parameter indices>)
tabular_block <- function(models, case, groups) {
  S <- lapply(models, stats, case = case); nm <- length(models)
  head <- c("\\toprule",
            paste0("& ", paste(sprintf("\\multicolumn{2}{c}{%s}", models), collapse = " & "), " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_len(nm), 2 * seq_len(nm) + 1), collapse = " "),
            paste0("& ", paste(rep("Mean & RMSE", nm), collapse = " & "), " \\\\"),
            paste0("& ", paste(rep("(s.e.) & (Bias)", nm), collapse = " & "), " \\\\"), "\\midrule")
  body <- c()
  for (g in seq_along(groups)) {
    idx <- groups[[g]]$idx
    for (j in seq_along(idx)) {
      i <- idx[j]
      lab <- if (j == 1) sprintf("\\multirow{%d}{*}{%s}", 2 * length(idx), groups[[g]]$label) else ""
      l1 <- paste(sapply(S, function(s) sprintf("%s & %s", f4(s$mean[i]), f4(s$rmse[i]))), collapse = " & ")
      l2 <- paste(sapply(S, function(s) sprintf("(%s) & (%s)", f4(s$se[i]), f4(s$bias[i]))), collapse = " & ")
      body <- c(body, sprintf("%s & %s \\\\", lab, l1), sprintf(" & %s \\\\", l2))
      if (j < length(idx)) body <- c(body, sprintf("\\cmidrule{2-%d}", 2 * nm + 1))
    }
    body <- c(body, if (g < length(groups)) "\\midrule" else paste0("\\hhline{", strrep("=", 2 * nm + 1), "}"))
  }
  c(sprintf("\\begin{tabular}{l%s}", strrep("c", 2 * nm)), head, body, "\\end{tabular}")
}
two_panel_table <- function(left, right, caption, label, widths = c(0.65, 0.65)) {
  c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\resizebox{\\textwidth}{!}{%",
    sprintf("    \\begin{minipage}{%s\\textwidth}", widths[1]), "        \\centering", paste0("        ", left), "    \\end{minipage}",
    "    \\hspace{0.01\\textwidth}",
    sprintf("    \\begin{minipage}{%s\\textwidth}", widths[2]), "        \\centering", paste0("        ", right), "    \\end{minipage}",
    "}", "}", sprintf("\\caption{%s}\\label{%s}", caption, label), "\\end{table}")
}
g_sev  <- list(list(label = "$\\beta_Y$", idx = 1:3), list(label = "$\\alpha$", idx = 4))
g_sevJ <- c(g_sev, list(list(label = "$\\rho_Y$", idx = 5)))
g_frq  <- list(list(label = "$\\beta_N$", idx = 1:3), list(label = "$r$", idx = 4))
g_frqJ <- c(g_frq, list(list(label = "$\\rho_N$", idx = 5)))
g_D    <- list(list(label = "$\\kappa$", idx = 1:4), list(label = "$\\beta_D$", idx = 5:6))

for (cs in cases) {
  tab <- two_panel_table(tabular_block(c("SEM1", "SEM2", "SEM3"), cs, g_sev), tabular_block(c("FEM1", "FEM2G", "FEM2"), cs, g_frq),
                         sprintf("Simulation result of SEM and FEM in Case %d (discrete deductible)", cs), sprintf("tab:AppSimExoCase%d", cs))
  writeLines(tab, file.path(texdir, sprintf("tab_app_case%d_marginal.tex", cs)))
  tab <- two_panel_table(tabular_block(c("SDM1", "SDM2", "SDM3"), cs, g_sevJ), tabular_block(c("FDM1", "FDM2"), cs, g_frqJ),
                         sprintf("Simulation result of SDM and FDM in Case %d (discrete deductible)", cs), sprintf("tab:AppSimEndCase%d", cs),
                         widths = c(0.65, 0.5))
  writeLines(tab, file.path(texdir, sprintf("tab_app_case%d_joint.tex", cs)))
}
## deductible model: one table, one column pair per case
{
  S <- lapply(cases, function(cs) stats("D", cs)); nm <- length(cases)
  head <- c("\\toprule", paste0("& ", paste(sprintf("\\multicolumn{2}{c}{Case %d}", cases), collapse = " & "), " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_len(nm), 2 * seq_len(nm) + 1), collapse = " "),
            paste0("& ", paste(rep("Mean & RMSE", nm), collapse = " & "), " \\\\"),
            paste0("& ", paste(rep("(s.e.) & (Bias)", nm), collapse = " & "), " \\\\"), "\\midrule")
  body <- c(); labs <- c(sprintf("$\\kappa_%d$", 1:4), "$\\beta_{D,1}$", "$\\beta_{D,2}$")
  for (i in 1:6) {
    l1 <- paste(sapply(S, function(s) sprintf("%s & %s", f4(s$mean[i]), f4(s$rmse[i]))), collapse = " & ")
    l2 <- paste(sapply(S, function(s) sprintf("(%s) & (%s)", f4(s$se[i]), f4(s$bias[i]))), collapse = " & ")
    body <- c(body, sprintf("%s & %s \\\\", labs[i], l1), sprintf(" & %s \\\\", l2), if (i < 6) sprintf("\\cmidrule{2-%d}", 2 * nm + 1))
  }
  tab <- c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\resizebox{\\textwidth}{!}{%", sprintf("\\begin{tabular}{l%s}", strrep("c", 2 * nm)), head, body,
           paste0("\\hhline{", strrep("=", 2 * nm + 1), "}"), "\\end{tabular}", "}", "}",
           "\\caption{Simulation result of the ordinal deductible model (thresholds $\\kappa_1<\\dots<\\kappa_4$ and $\\bm{\\beta}_D$) in Cases 1--4}\\label{tab:AppSimDed}", "\\end{table}")
  writeLines(tab, file.path(texdir, "tab_app_D.tex"))
}

## ---- MAE tables (main text) -----------------------------------------------------------------------
metric <- function(case, key) mean(sapply(res[sapply(res, `[[`, "case") == case], function(r) r$metrics[[key]]))
mae_table <- function(cs) {
  fq <- list(cen = c("FEM1", "FDM1"), tru = c("FEM1", "FEM2", "FDM1", "FDM2"))
  sv <- list(cen = c("SEM1", "SEM2", "SDM1", "SDM2"), tru = c("SEM1", "SEM3", "SDM1", "SDM3"))
  row <- function(q, mods, type, dig) paste(sapply(mods, function(m) formatC(metric(cs, paste(q, m, type, sep = ".")), format = "f", digits = dig)), collapse = " & ")
  c("\\begin{table}[!h]", "\\renewcommand{\\tabcolsep}{0.8mm}", "\\renewcommand{\\arraystretch}{1}", "\\footnotesize", "\\centering", "{\\color{blue}",
    "\\begin{tabular}[h]{ccccccc}", "\\toprule",
    "\\multirow{2}{*}{Frequency} & \\multicolumn{2}{c}{Censored} & \\multicolumn{4}{c}{Truncated}\\\\",
    "\\cmidrule(l{3pt}r{3pt}){2-3} \\cmidrule(l{3pt}r{3pt}){4-7}",
    " & FEM1 & FDM1 & FEM1 & FEM2 & FDM1 & FDM2\\\\", "\\midrule",
    sprintf("MAE & %s & %s \\\\", row("N_cen", fq$cen, "obs", 4), row("N_tru", fq$tru, "obs", 4)),
    if (cs <= 4) sprintf("MAE$_{\\text{prem}}$ & %s & %s \\\\", row("N_cen", fq$cen, "prem", 4), row("N_tru", fq$tru, "prem", 4)),
    "\\hhline{=======}", "\\end{tabular}", "\\\\", "\\vspace{0.3cm}",
    "\\begin{tabular}[h]{ccccccccc}", "\\toprule",
    "\\multirow{2}{*}{Severity} & \\multicolumn{4}{c}{Censored} & \\multicolumn{4}{c}{Truncated}\\\\",
    "\\cmidrule(l{3pt}r{3pt}){2-5} \\cmidrule(l{3pt}r{3pt}){6-9}",
    " & SEM1 & SEM2 & SDM1 & SDM2 & SEM1 & SEM3 & SDM1 & SDM3\\\\", "\\midrule",
    sprintf("MAE & %s & %s \\\\", row("Y_cen", sv$cen, "obs", 1), row("Y_tru", sv$tru, "obs", 1)),
    if (cs <= 4) sprintf("MAE$_{\\text{prem}}$ & %s & %s \\\\", row("Y_cen", sv$cen, "prem", 1), row("Y_tru", sv$tru, "prem", 1)),
    "\\hhline{=========}", "\\end{tabular}", "}",
    sprintf("\\caption{MAE of the frequency and severity models in Case %d}\\label{tab:maecase%d}", cs, cs), "\\end{table}")
}
for (cs in cases) writeLines(mae_table(cs), file.path(texdir, sprintf("tab_mae_case%d.tex", cs)))

## ---- coverage of 95% intervals: naive vs robust vs sandwich s.e. ----------------------------------------
cov_rows <- list()
for (cs in cases) {
  rr <- res[sapply(res, `[[`, "case") == cs]; rr <- rr[!sapply(rr, function(r) is.null(r$se))]
  if (!length(rr)) next
  for (mod in names(rr[[1]]$se)) {
    th <- truth(mod, cs); pn <- rr[[1]]$se[[mod]]$parameter
    Etau <- t(sapply(rr, function(r) r$se[[mod]]$estimate)); E <- Etau
    if ("tau" %in% pn) { j <- which(pn == "tau"); E[, j] <- tau2rho(Etau[, j]) }                        # report rho = sin(pi tau / 2)
    for (type in c("se_naive", "se_robust", "se_sandwich")) {
      SE <- t(sapply(rr, function(r) r$se[[mod]][[type]]))
      if ("tau" %in% pn) SE[, j] <- SE[, j] * (pi / 2) * cos(pi * Etau[, j] / 2)                       # delta method tau -> rho
      cover <- colMeans(abs(E - matrix(th, nrow(E), ncol(E), byrow = TRUE)) <= qnorm(0.975) * SE)
      cov_rows[[length(cov_rows) + 1]] <- data.frame(case = cs, model = mod, type = sub("se_", "", type), parameter = sub("tau", "rho", pn),
                                                     mean_se = colMeans(SE), emp_sd = apply(E, 2, sd), coverage = cover, R = nrow(E))
    }
  }
}
cov <- do.call(rbind, cov_rows); write.csv(cov, file.path(texdir, "coverage.csv"), row.names = FALSE)
## LaTeX: (case, model) rows: empirical sd, mean naive / robust / sandwich s.e., coverage
cov_tex <- function(rows) {
  lines <- c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\begin{tabular}{llccccccc}", "\\toprule",
             "Model & Parameter & Emp. SD & \\multicolumn{3}{c}{Mean s.e.} & \\multicolumn{3}{c}{Coverage of 95\\% interval} \\\\",
             "\\cmidrule(lr){4-6} \\cmidrule(lr){7-9}", " & & & Naive & Robust & Sandwich & Naive & Robust & Sandwich \\\\", "\\midrule")
  Rmin <- Inf
  for (k in seq_len(nrow(rows))) {
    cs <- rows$case[k]; mod <- rows$model[k]
    d <- cov[cov$case == cs & cov$model == mod, ]; if (!nrow(d)) next; pn <- unique(d$parameter); Rmin <- min(Rmin, d$R)
    for (i in seq_along(pn)) {
      dn <- d[d$type == "naive" & d$parameter == pn[i], ]; dr <- d[d$type == "robust" & d$parameter == pn[i], ]; ds <- d[d$type == "sandwich" & d$parameter == pn[i], ]
      sub_ <- if (substr(mod, 1, 1) == "S") "Y" else "N"
      lab <- switch(pn[i], "(Intercept)" = sprintf("$\\beta_{%s,0}$", sub_), X1 = sprintf("$\\beta_{%s,1}$", sub_), X2 = sprintf("$\\beta_{%s,2}$", sub_),
                    alpha = "$\\alpha$", r = "$r$", rho = sprintf("$\\rho_%s$", sub_))
      lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %.3f & %.3f & %.3f \\\\", if (i == 1) sprintf("%s (Case %d)", mod, cs) else "", lab,
                                f4(dn$emp_sd), f4(dn$mean_se), f4(dr$mean_se), f4(ds$mean_se), dn$coverage, dr$coverage, ds$coverage))
    }
    lines <- c(lines, "\\midrule")
  }
  lines[length(lines)] <- "\\hhline{=========}"
  c(lines, "\\end{tabular}", "}",
    sprintf("\\caption{Monte Carlo standard deviation of the sequential estimates, mean estimated standard errors (naive: stage-wise Hessian; robust: stage-wise sandwich $\\bm{H}_{kk}^{-1}\\bm{J}_{kk}\\bm{H}_{kk}^{-1}$; sandwich: Godambe information $\\bm{H}^{-1}\\bm{J}\\bm{H}^{-\\top}$ of Section 2.3) and empirical coverage of the nominal 95\\%% confidence intervals ($R=%d$ replications); the dependence parameters are on the correlation scale (delta method)}\\label{tab:simcoverage}", Rmin),
    "\\end{table}")
}
rows <- data.frame(case = c(1, 4, 4), model = c("FEM2", "SDM3", "FDM2"))
if (nrow(cov)) writeLines(cov_tex(rows), file.path(texdir, "tab_coverage.tex"))


## ---- relative error of the premium components by deductible level (Case 4) -------------------------------
bylev_table <- function(cs) {
  rr <- res[sapply(res, `[[`, "case") == cs]
  rel <- function(key) { M <- sapply(rr, function(r) { v <- r$metrics[[key]]; out <- rep(NA_real_, K_ded); out[as.integer(names(v))] <- v; out }); rowMeans(M, na.rm = TRUE) }
  row <- function(key) paste(sub("^-0\\.0$", "0.0", sprintf("%+.1f", 100 * rel(key))), collapse = " & ")
  c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\begin{tabular}{llrrrrr}", "\\toprule",
    " & & \\multicolumn{5}{c}{Deductible level} \\\\", "\\cmidrule(lr){3-7}",
    sprintf(" & Model & %s \\\\", paste(format(design$levels, big.mark = ","), collapse = " & ")), "\\midrule",
    sprintf("\\multirow{2}{*}{Frequency, censored} & FEM1 & %s \\\\", row("N_cen.FEM1.bylev")), sprintf(" & FDM1 & %s \\\\", row("N_cen.FDM1.bylev")), "\\midrule",
    sprintf("\\multirow{4}{*}{Frequency, truncated} & FEM1 & %s \\\\", row("N_tru.FEM1.bylev")), sprintf(" & FEM2 & %s \\\\", row("N_tru.FEM2.bylev")),
    sprintf(" & FDM1 & %s \\\\", row("N_tru.FDM1.bylev")), sprintf(" & FDM2 & %s \\\\", row("N_tru.FDM2.bylev")), "\\midrule",
    sprintf("\\multirow{4}{*}{Severity, censored} & SEM1 & %s \\\\", row("Y_cen.SEM1.bylev")), sprintf(" & SEM2 & %s \\\\", row("Y_cen.SEM2.bylev")),
    sprintf(" & SDM1 & %s \\\\", row("Y_cen.SDM1.bylev")), sprintf(" & SDM2 & %s \\\\", row("Y_cen.SDM2.bylev")), "\\midrule",
    sprintf("\\multirow{4}{*}{Severity, truncated} & SEM1 & %s \\\\", row("Y_tru.SEM1.bylev")), sprintf(" & SEM3 & %s \\\\", row("Y_tru.SEM3.bylev")),
    sprintf(" & SDM1 & %s \\\\", row("Y_tru.SDM1.bylev")), sprintf(" & SDM3 & %s \\\\", row("Y_tru.SDM3.bylev")),
    "\\hhline{=======}", "\\end{tabular}", "}",
    sprintf("\\caption{Relative error (in \\%%) of the estimated frequency and severity components of the pure premium by deductible level in Case %d: the average over the policy-years (frequency) or claims (severity) at each level of the difference between the estimated and the true component, divided by the average true component at that level, averaged over the replications}\\label{tab:bylevel%d}", cs, cs),
    "\\end{table}")
}
if (4 %in% cases) writeLines(bylev_table(4), file.path(texdir, "tab_bylevel_case4.tex"))

## ---- figures: box plots of the estimates by model, one panel per parameter ---------------------------------
pdf_boxes <- function(cs) {
  pdf(file.path(texdir, sprintf("simcase%d.pdf", cs)), width = 11, height = 6.5)
  par(mfrow = c(2, 5), mar = c(4.2, 3.2, 2.2, 0.6), mgp = c(2, 0.7, 0), cex.axis = 0.8)
  sev <- c("SEM1", "SEM2", "SEM3", "SDM1", "SDM2", "SDM3"); frq <- c("FEM1", "FEM2G", "FEM2", "FDM1", "FDM2")
  box <- function(models, i, main, th) {
    L <- lapply(models, function(m) { P <- par_mat(m, cs); if (i <= ncol(P)) P[, i] else NA })
    boxplot(L, names = models, las = 2, main = main, outline = TRUE, pch = 20, cex = 0.5, col = ifelse(grepl("DM", models), "lightsteelblue", "grey90"))
    abline(h = th, col = "red", lty = 2)
  }
  rho <- design$cases[[as.character(cs)]]
  box(sev, 1, expression(beta[Y*","*0]), design$beta_Y[1]); box(sev, 2, expression(beta[Y*","*1]), design$beta_Y[2]); box(sev, 3, expression(beta[Y*","*2]), design$beta_Y[3])
  box(sev, 4, expression(alpha), design$alpha); box(sev[4:6], 5, expression(rho[Y]), rho[["rho_Y"]])
  box(frq, 1, expression(beta[N*","*0]), design$beta_N[1]); box(frq, 2, expression(beta[N*","*1]), design$beta_N[2]); box(frq, 3, expression(beta[N*","*2]), design$beta_N[3])
  box(frq, 4, expression(r), design$r); box(frq[4:5], 5, expression(rho[N]), rho[["rho_N"]])
  dev.off()
}
for (cs in cases) pdf_boxes(cs)

## ---- key numbers for the text ----------------------------------------------------------------------------
sink(file.path(texdir, "key_numbers.txt"))
cat("replications per case:\n"); print(R)
cat("\nconvergence failures (any model):", sum(sapply(res, function(r) any(r$convergence != 0))), "\n")
cat("mean number of claims:", mean(sapply(res, `[[`, "n_claims")), " mean number of payments (Y*>d):", mean(sapply(res, `[[`, "n_pos")), "\n")
cat("fraction of losses below the deductible by case:\n"); print(tapply(sapply(res, function(r) 1 - r$n_pos / r$n_claims), sapply(res, `[[`, "case"), mean))
cat("mean level proportions:\n"); print(round(Reduce("+", lapply(res, `[[`, "lev_prop")) / length(res), 3))
cat("mean Kendall tau (level vs N*, level vs Y*) by case:\n"); print(t(sapply(cases, function(cs) colMeans(t(sapply(res[sapply(res, `[[`, "case") == cs], `[[`, "kendall"))))))
cat("mean time per replication (s) by case:\n"); print(tapply(sapply(res, `[[`, "time_total"), sapply(res, `[[`, "case"), mean))
for (cs in cases) { cat(sprintf("\n==== Case %d ====\n", cs))
  for (m in c("D", "SEM1", "SEM2", "SEM3", "FEM1", "FEM2G", "FEM2", "SDM1", "SDM2", "SDM3", "FDM1", "FDM2")) { cat(m, "\n"); print(stats(m, cs), digits = 4) }
  cat("MAE / premium MAE:\n"); keys <- grep("obs$|prem$", names(res[[1]]$metrics), value = TRUE)
  print(round(sapply(keys, function(k) metric(cs, k)), 4))
  cat("relative premium error by deductible level (mean over replications):\n")
  keys <- grep("bylev$", names(res[[1]]$metrics), value = TRUE)
  print(round(t(sapply(keys, function(k) rowMeans(sapply(res[sapply(res, `[[`, "case") == cs], function(r) { v <- r$metrics[[k]]; out <- rep(NA_real_, K_ded); out[as.integer(names(v))] <- v; out }), na.rm = TRUE))), 4))
}
if (nrow(cov)) { cat("\ncoverage:\n"); print(cov, digits = 4) }
sink()
cat("written to", texdir, "\n")
