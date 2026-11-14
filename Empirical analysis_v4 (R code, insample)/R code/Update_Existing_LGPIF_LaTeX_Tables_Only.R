################################################################################
# UPDATE BOTH LGPIF MANUSCRIPT TABLES FROM NEWLY FITTED RESULTS
# Run AFTER: Empirical_2006_2009_Training_2010_OOS_Cases1_4_Combined.R
# No refitting. Base R only. Original templates are never overwritten.
# This script runs separately from your completed empirical analysis.
################################################################################

# -------- USER CONFIGURATION ---------------------------------------------------
results_dir <- "empirical_merged_outputs"
main_tex_file <- "Empirical_Main_Template.tex"
appendix_tex_file <- "Empirical_Appendix_Template.tex"
# If your .tex files have other names, change these two paths above.
# ------------------------------------------------------------------------------

stop_if_missing <- function(path, description = "Required file") {
  if (!file.exists(path)) stop(description, " not found: ", path,
                               "\nRun the merged empirical pipeline first, or edit the paths above.")
  path
}
read_csv <- function(path) utils::read.csv(stop_if_missing(path),
                                             stringsAsFactors = FALSE,
                                             check.names = FALSE)
res_path <- file.path(results_dir, "Empirical_Revised_Results.rds")
res <- readRDS(stop_if_missing(res_path, "Training results"))
stopifnot(is.data.frame(res$premium_model_pairs))
stopifnot(nrow(res$premium_model_pairs) == 9L,
          identical(as.character(unique(res$premium_model_pairs$Case)),
                    c("1", "4c", "4d")))

output_dir <- file.path(results_dir, "manuscript_tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
BS <- "\\"
NL <- paste0(" ", BS, BS)
latex <- function(cmd) paste0(BS, cmd)
fmt <- function(z, digits = 2L) {
  if (length(z) != 1L || is.na(z) || !is.finite(z)) return("--")
  formatC(as.numeric(z), digits = digits, format = "f", big.mark = ",")
}
tex_param <- function(p) {
  if (is.na(p) || !nzchar(p)) return("")
  switch(p,
         "(Intercept)" = "Intercept",
         "sigma" = "$\\sigma$", "alpha1" = "$\\alpha_1$",
         "alpha2" = "$\\alpha_2$", "tau" = "$\\tau$",
         "nu" = "$\\nu$", "r" = "$r$",
         "log(Coverage)" = "$\\log(\\mathrm{Coverage})$",
         p)
}
make_tab <- function(colspec, lines, caption, label, wide = FALSE) {
  content <- c(
    paste0(latex("begin{table}[h!]")),
    latex("centering"),
    if (wide) paste0(latex("resizebox{"), latex("linewidth}{!}{")),
    paste0(latex("begin{tabular}{"), colspec, "}"),
    latex("toprule"),
    lines,
    latex("bottomrule"),
    latex("end{tabular}"),
    if (wide) "}",
    paste0(latex("caption{"), caption, "}"),
    paste0(latex("label{"), label, "}"),
    latex("end{table}")
  )
  paste(content, collapse = "\n")
}
tr <- function(...) paste0(paste(c(...), collapse = " & "), NL)

# -------- Training data (only descriptive statistics) -------------------------
data_dir <- res$config$data_dir
penv <- new.env(parent = baseenv())
cenv <- new.env(parent = baseenv())
load(stop_if_missing(file.path(data_dir, "poldat.RData")), envir = penv)
load(stop_if_missing(file.path(data_dir, "claimdat.RData")), envir = cenv)
if (!exists("poldat", envir = penv, inherits = FALSE) ||
    !exists("claimdat", envir = cenv, inherits = FALSE))
  stop("Expected poldat and claimdat objects in the source .RData files.")
pol <- get("poldat", envir = penv)
clm <- get("claimdat", envir = cenv)
train_years <- res$config$years_in
pol <- pol[pol$Year %in% train_years, , drop = FALSE]
clm <- clm[clm$Year %in% train_years, , drop = FALSE]
if (!nrow(pol) || !nrow(clm)) stop("Empty training policy/claim sample")
if (anyDuplicated(paste(pol$PolicyNum, pol$Year)))
  stop("Duplicate policy-year key in original training data")
clm$Payment <- pmax(0, clm$LossBeforeDeductible - clm$Deduct)
if (!all(paste(clm$PolicyNum, clm$Year) %in% paste(pol$PolicyNum, pol$Year)))
  warning("Some training claims do not match training policy-year records")
claims_paid <- aggregate(as.integer(clm$Payment > 0),
                        by = list(PolicyNum = clm$PolicyNum, Year = clm$Year),
                        FUN = sum)
pol_key <- paste(pol$PolicyNum, pol$Year)
claim_key <- paste(claims_paid$PolicyNum, claims_paid$Year)
Nstar <- claims_paid$x[match(pol_key, claim_key)]
Nstar[is.na(Nstar)] <- 0
stat <- function(x) c(Mean = if (length(x)) mean(x) else NA_real_,
                      SD = if (length(x) > 1) sd(x) else NA_real_)

# -------- Fits and standard errors -------------------------------------------
fits <- c(res$marginal_severity, res$selected_severity, res$frequency,
          list(D_cont = res$deductible$continuous, D_disc = res$deductible$discrete))
se_for <- function(nm, fit) {
  sw <- res$godambe[[nm]]
  if (!is.null(sw) && !is.null(sw$idx) && !is.null(sw$table)) {
    ix <- sw$idx[[length(sw$idx)]]
    tb <- sw$table[ix, , drop = FALSE]
    vals <- tb$se_sandwich[match(names(fit$par), tb$parameter)]
    if (length(vals) == length(fit$par) && all(is.finite(vals)))
      return(setNames(vals, names(fit$par)))
  }
  if (!is.null(fit$se_naive)) {
    ss <- fit$se_naive
    return(setNames(as.numeric(ss), names(fit$par)))
  }
  setNames(rep(NA_real_, length(fit$par)), names(fit$par))
}
se_complete <- function(nms) all(vapply(nms, function(nm) {
  sw <- res$godambe[[nm]]
  !is.null(sw) && !is.null(sw$idx) && !is.null(sw$table)
}, logical(1)))
se_caption <- function(nms) {
  if (se_complete(nms)) "Godambe sandwich standard errors"
  else "standard errors (Godambe where available; otherwise stage-wise)"
}
parm_table <- function(models, caption, label, digits = 3L) {
  fs <- fits[models]
  if (any(vapply(fs, is.null, logical(1))))
    stop("Missing parameter model(s): ", paste(models[vapply(fs, is.null, logical(1))], collapse = ", "))
  par_order <- unique(unlist(lapply(fs, function(f) names(f$par))))
  hdr <- c("Parameter")
  for (nm in models) hdr <- c(hdr, paste0(nm, " estimate"), "s.e.")
  h <- tr(hdr)
  rows <- vapply(par_order, function(parameter) {
    entries <- c(tex_param(parameter))
    for (nm in models) {
      f <- fs[[nm]]
      se <- se_for(nm, f)
      entries <- c(entries,
        fmt(if (parameter %in% names(f$par)) unname(f$par[parameter]) else NA_real_, digits),
        fmt(if (parameter %in% names(se)) unname(se[parameter]) else NA_real_, digits))
    }
    do.call(tr, as.list(entries))
  }, character(1))
  make_tab(paste0("l", paste(rep("rr", length(models)), collapse = "")),
           c(h, latex("midrule"), rows), caption, label, wide = length(models) > 3L)
}

main_tables <- list()
appendix_tables <- list()

# Main: frequency and severity by observation scheme
ss <- rbind(N = stat(pol$Freq), Nstar = stat(Nstar),
            Y = stat(clm$LossBeforeDeductible),
            Ycensored = stat(clm$Payment),
            Ytruncated = stat(clm$Payment[clm$Payment > 0]))
main_tables[["tab:data summary"]] <- make_tab("lrrrrr", c(
  tr("", paste0(latex("multicolumn{2}{c}{Frequency}")),
     paste0(latex("multicolumn{3}{c}{Severity}"))),
  latex("cmidrule(lr){2-3}"), latex("cmidrule(lr){4-6}"),
  tr("", "$N$", "$N^*$", "$Y$", "$(Y-d)_+$", "$Y-d\\mid Y>d$"),
  latex("midrule"),
  do.call(tr, as.list(c("Mean", vapply(ss[, "Mean"], fmt, character(1), digits = 2L)))),
  do.call(tr, as.list(c("SD", vapply(ss[, "SD"], fmt, character(1), digits = 2L))))
), "Summary statistics of frequency and severity under different observation scenarios (2006--2009).",
"tab:data summary")

main_tables[["marginalY"]] <- parm_table(
  c("SEM1", "SEM2", "SEM3"),
  paste0("Regression estimates and ", se_caption(c("SEM1", "SEM2", "SEM3")),
         " of the marginal severity models."), "marginalY")
main_tables[["tab:MFreq_LGPIF"]] <- parm_table(
  c("FEM1", "FEM2_SEM3"),
  paste0("Regression estimates and ", se_caption(c("FEM1", "FEM2_SEM3")),
         " of the marginal frequency models (FEM2 denotes FEM2_SEM3)."),
  "tab:MFreq_LGPIF")

# Main: selected copula dependence estimates
selected_names <- c(paste0("SCM", 1:3), paste0("SDM", 1:3),
                    paste0("FCM", 1:2), paste0("FDM", 1:2))
dep_rows <- vapply(selected_names, function(nm) {
  f <- fits[[nm]]
  v <- res$selected_severity_table
  if (startsWith(nm, "F")) v <- res$selected_frequency_table
  row <- v[v$Model == nm, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing selected copula for ", nm)
  tau <- if ("tau" %in% names(f$par)) unname(f$par["tau"]) else 0
  se <- se_for(nm, f)
  tau_label <- if ("tau" %in% names(se))
    paste0(fmt(tau, 3L), " (", fmt(unname(se["tau"]), 3L), ")")
    else fmt(tau, 3L)
  tr(nm, row$Family, tau_label)
}, character(1))
main_tables[["tab:lgpif_dependence"]] <- make_tab("llr", c(
  tr("Model", "Selected family", "$\\widehat{\\tau}$ (s.e.)"),
  latex("midrule"), dep_rows),
  "Selected copula families and dependence estimates for the LGPIF dataset.",
  "tab:lgpif_dependence")

main_tables[["tab:Sevcopularesult"]] <- parm_table(
  c(paste0("SCM", 1:3), paste0("SDM", 1:3)),
  paste0("Regression estimates and ", se_caption(c(paste0("SCM", 1:3), paste0("SDM", 1:3))),
         " of the selected joint severity models."),
  "tab:Sevcopularesult", digits = 3L)
main_tables[["tab:freqcopularesult"]] <- parm_table(
  c("FCM1", "FCM2", "FDM1", "FDM2"),
  paste0("Regression estimates and ", se_caption(c("FCM1", "FCM2", "FDM1", "FDM2")),
         " of the selected joint frequency models."),
  "tab:freqcopularesult", digits = 3L)

# Main: common latent ground-up N GOF table
freq_gof_path <- file.path(results_dir, "Frequency_GOF_CommonGroundupN.csv")
if (file.exists(freq_gof_path)) {
  freq_gof <- read_csv(freq_gof_path)
  models <- c("FEM1", "FEM2", "FCM1", "FCM2", "FDM1", "FDM2")
  counts <- unique(freq_gof$Count[freq_gof$Model == "FEM1"])
  rowfun <- function(count) {
    sub <- freq_gof[freq_gof$Count == count, , drop = FALSE]
    m <- match(models, sub$Model)
    if (anyNA(m)) stop("Incomplete common N GOF results")
    tr(count, fmt(sub$Observed[m[1L]], 0L),
       vapply(sub$Expected[m], fmt, character(1), digits = 2L))
  }
  # tr() accepts vector arguments via c(); avoid argument/list column folding
  rows <- vapply(counts, rowfun, character(1))
  last <- freq_gof[match(models, freq_gof$Model), , drop = FALSE]
  chi <- do.call(tr, as.list(c("$\\chi^2$ statistic", "--",
                             vapply(last$Pearson_Chi2_descriptive, fmt, character(1), digits = 2L))))
  main_tables[["tab:Freqcopula"]] <- make_tab("lrrrrrrr", c(
    tr("$N$", "Empirical", models), latex("midrule"), rows,
    latex("midrule"), chi),
    "Goodness-of-fit for the common ground-up frequency distribution $N$ (descriptive Pearson statistics).",
    "tab:Freqcopula", wide = TRUE)
} else {
  warning("Common-N GOF output missing. Frequency GOF table is left unchanged.")
}

# Main: pure premiums, in-sample and two OOS conditions.
score_files <- c(
  "tab:insample" = file.path(results_dir, "manuscript_reporting",
                             "Table3_Limited_PurePremium_Main.csv"),
  "tab:oos1" = file.path(results_dir, "out_of_sample_2010",
                         "manuscript_reporting", "Table3_2010_Limited_PurePremium_Main.csv"),
  "tab:oos2" = file.path(results_dir, "out_of_sample_2010_exclude_max_loss",
                         "manuscript_reporting", "Table3_2010_Limited_PurePremium_Main.csv")
)
premium_tab <- function(tab, label, caption) {
  ord <- do.call(rbind, lapply(c("1", "4c", "4d"), function(ca) {
    ss <- tab[as.character(tab$Case) == ca, , drop = FALSE]
    ss[order(as.integer(sub(".*([123])$", "\\1", ss$Scheme))), , drop = FALSE]
  }))
  if (nrow(ord) != 9L || anyDuplicated(ord$Model))
    stop("Expected nine distinct Case 1/4c/4d models for ", label)
  metrics <- c("MAE_obs", "RMSE_obs", "OP_ratio", "Mean_Poisson_Deviance")
  row_names <- c("MAE", "RMSE", "O/P (\\%)", "MPD")
  rows <- vapply(seq_along(metrics), function(ii) {
    v <- ord[[metrics[ii]]]
    if (identical(metrics[ii], "OP_ratio")) v <- 100 * v
    do.call(tr, as.list(c(row_names[ii], vapply(v, fmt, character(1), digits = 2L))))
  }, character(1))
  make_tab("lrrrrrrrrr", c(
    tr("", paste0(latex("multicolumn{3}{c}{Case 1: Exogenous}")),
       paste0(latex("multicolumn{3}{c}{Case 4c: Continuous}")),
       paste0(latex("multicolumn{3}{c}{Case 4d: Discrete}"))),
    latex("cmidrule(lr){2-4}"), latex("cmidrule(lr){5-7}"), latex("cmidrule(lr){8-10}"),
    tr("", "XR1", "XR2", "XR3", "ER1", "ER2", "ER3", "ER1", "ER2", "ER3"),
    latex("midrule"), rows), caption, label, wide = TRUE)
}
score_captions <- c("tab:insample" = "In-sample coverage-limited pure premium comparison (training years 2006--2009).",
                    "tab:oos1" = "Out-of-sample coverage-limited pure premium comparison (all 2010 policy-years).",
                    "tab:oos2" = "Out-of-sample sensitivity analysis excluding the single 2010 policy-year with the largest observed aggregate limited loss.")
for (lab in names(score_files)) {
  if (file.exists(score_files[[lab]])) {
    main_tables[[lab]] <- premium_tab(read_csv(score_files[[lab]]), lab, score_captions[[lab]])
  } else {
    warning("No scores for ", lab, ": ", score_files[[lab]], "; existing table unchanged")
  }
}

# Appendix: deductible models in two-column layout
Dcont <- fits$D_cont
Ddisc <- fits$D_disc
cc <- length(Dcont$par)
dd <- length(Ddisc$par)
param_d_row <- function(j) {
  left <- if (j <= cc) c(tex_param(names(Dcont$par)[j]),
                        fmt(Dcont$par[j], 3L), fmt(se_for("D_cont", Dcont)[j], 3L))
          else c("", "", "")
  right <- if (j <= dd) c(tex_param(names(Ddisc$par)[j]),
                         fmt(Ddisc$par[j], 3L), fmt(se_for("D_disc", Ddisc)[j], 3L))
           else c("", "", "")
  do.call(tr, as.list(c(left, right)))
}
appendix_tables[["tab:Deduct"]] <- make_tab("lrr@{\\hspace{1cm}}lrr", c(
  tr(paste0(latex("multicolumn{3}{c}{Continuous deductible}")),
     paste0(latex("multicolumn{3}{c}{Discrete deductible}"))),
  latex("cmidrule(lr){1-3}"), latex("cmidrule(lr){4-6}"),
  tr("Parameter", "Estimate", "s.e.", "Parameter", "Estimate", "s.e."),
  latex("midrule"), vapply(seq_len(max(cc, dd)), param_d_row, character(1))
), "Regression estimates and stage-wise standard errors of deductible models.",
"tab:Deduct", wide = TRUE)

# Appendix: original policy/claim descriptive table
vtypes <- c("TypeCity", "TypeCounty", "TypeSchool", "TypeTown", "TypeVillage")
pstat <- list(Deduct = stat(pol$Deduct), Dedratio = stat(pol$Deduct / pol$Coverage),
              Freq = stat(pol$Freq), Claim = stat(clm$LossBeforeDeductible),
              Coverage = stat(pol$Coverage))
app_rows <- c(
  tr("Deduct", "Amount of deductible in dollars", fmt(pstat$Deduct["Mean"], 3L), fmt(pstat$Deduct["SD"], 3L)),
  tr("Dedratio", "Deductible divided by coverage (coverage in million dollars)",
     fmt(pstat$Dedratio["Mean"], 3L), fmt(pstat$Dedratio["SD"], 3L)),
  tr("Freq", "Number of ground-up claims per policy-year", fmt(pstat$Freq["Mean"], 3L), fmt(pstat$Freq["SD"], 3L)),
  tr("Claim", "Ground-up amount of individual claims (USD)", fmt(pstat$Claim["Mean"], 2L), fmt(pstat$Claim["SD"], 2L)),
  tr("Entity type", "", "", ""),
  vapply(vtypes, function(nm) {
    if (!nm %in% names(pol)) stop("Missing entity indicator: ", nm)
    tr(sub("Type", "", nm), paste("Indicator for", sub("Type", "", nm)),
       fmt(mean(pol[[nm]]), 3L), "")
  }, character(1)),
  tr("Coverage", "Coverage in millions of dollars",
     fmt(pstat$Coverage["Mean"], 3L), fmt(pstat$Coverage["SD"], 3L))
)
appendix_tables[["tab:lgpif summary"]] <- make_tab("llrr", c(
  tr("Variable", "Description", "Mean", "SD"), latex("midrule"), app_rows),
  "Description and summary statistics of response and explanatory variables (training years 2006--2009).",
  "tab:lgpif summary")

# Appendix: 4 AIC copula candidate tables, generated from full selection results
candidate_tab <- function(component, dedtype, models, label, caption) {
  scenarios <- if (component == "sev") c("complete", "censored", "truncated")
               else c("censored", "truncated")
  keys <- paste(component, dedtype, scenarios, sep = "_")
  sel <- res$copula_selection
  fitsel <- sel[keys]
  if (any(vapply(fitsel, is.null, logical(1))))
    stop("Missing copula selection candidates: ", paste(keys, collapse = ", "))
  names_by_key <- lapply(fitsel, function(z) rownames(z$table))
  families <- unique(unlist(names_by_key))
  header <- c("Family", rep(c("logLik", "AIC", "BIC", "$\\tau$"), length(keys)))
  rows <- vapply(families, function(family) {
    vals <- c(family)
    for (key in keys) {
      z <- sel[[key]]$table
      pos <- match(family, rownames(z))
      vals <- c(vals, if (is.na(pos)) rep("--", 4L)
                else vapply(z[pos, c("Likelihood", "AIC", "BIC", "tau")],
                            fmt, character(1), digits = 2L))
    }
    do.call(tr, as.list(vals))
  }, character(1))
  make_tab(paste0("l", paste(rep("rrrr", length(keys)), collapse = "")),
    c(tr(paste0("\\multicolumn{", 1 + 4L * length(keys), "}{c}{",
                paste(models, collapse = ", "), "}")),
      latex("midrule"), do.call(tr, as.list(header)),
      latex("midrule"), rows), caption, label, wide = TRUE)
}
appendix_tables[["tab:fcmcopula"]] <- candidate_tab("freq", "cont", c("FCM1", "FCM2"),
  "tab:fcmcopula", "Copula candidate likelihood and information criteria: continuous deductible frequency models.")
appendix_tables[["tab:fdmcopula"]] <- candidate_tab("freq", "disc", c("FDM1", "FDM2"),
  "tab:fdmcopula", "Copula candidate likelihood and information criteria: discrete deductible frequency models.")
appendix_tables[["tab:scmcopula"]] <- candidate_tab("sev", "cont", paste0("SCM", 1:3),
  "tab:scmcopula", "Copula candidate likelihood and information criteria: continuous deductible severity models.")
appendix_tables[["tab:sdmcopula"]] <- candidate_tab("sev", "disc", paste0("SDM", 1:3),
  "tab:sdmcopula", "Copula candidate likelihood and information criteria: discrete deductible severity models.")

# -------- Write independent LaTeX table fragments ----------------------------
for (lab in names(main_tables)) {
  fname <- paste0("MAIN_", gsub("[^A-Za-z0-9]+", "_", lab), ".tex")
  writeLines(main_tables[[lab]], file.path(output_dir, fname), useBytes = TRUE)
}
for (lab in names(appendix_tables)) {
  fname <- paste0("APPENDIX_", gsub("[^A-Za-z0-9]+", "_", lab), ".tex")
  writeLines(appendix_tables[[lab]], file.path(output_dir, fname), useBytes = TRUE)
}

# -------- Substitute by EXACT table labels; preserve original captions --------
# Replacements run from the end of the document to the beginning, so changes in
# one table cannot invalidate character positions in the remaining tables.
# Captions (including their formatting) are copied from the original table.
extract_caption <- function(block) {
  at <- regexpr("\\\\caption\\\\{", block, perl = TRUE)[1L]
  if (at < 0L) return(NULL)
  open <- at + nchar("\\caption")
  chars <- strsplit(block, "", fixed = TRUE)[[1L]]
  depth <- 0L
  for (j in seq.int(open, length(chars))) {
    if (chars[j] == "{" && (j == 1L || chars[j-1L] != "\\")) depth <- depth + 1L
    if (chars[j] == "}" && (j == 1L || chars[j-1L] != "\\")) {
      depth <- depth - 1L
      if (depth == 0L) return(substr(block, at, j))
    }
  }
  stop("Unbalanced caption braces")
}
replace_tables <- function(tex, replacement) {
  beginnings <- gregexpr("(?m)^[ \\t]*\\\\begin\\{table\\}", tex, perl = TRUE)[[1L]]
  endings <- gregexpr("(?m)^[ \\t]*\\\\end\\{table\\}", tex, perl = TRUE)[[1L]]
  if (length(beginnings) == 1L && beginnings[1] < 0L)
    stop("No active LaTeX table environments found")
  if (length(beginnings) != length(endings))
    stop("Unbalanced LaTeX table environments")
  done <- character()
  for (i in rev(seq_along(beginnings))) {
    first <- beginnings[i]
    last <- endings[i] + nchar("\\end{table}") - 1L
    old <- substr(tex, first, last)
    which_label <- names(replacement)[vapply(names(replacement), function(z)
      grepl(paste0("\\label{", z, "}"), old, fixed = TRUE), logical(1))]
    if (!length(which_label)) next
    if (length(which_label) != 1L) stop("Ambiguous table label at table ", i)
    new <- replacement[[which_label]]
    oldcap <- extract_caption(old)
    newcap <- extract_caption(new)
    if (!is.null(oldcap) && !is.null(newcap)) {
      cp <- regexpr(newcap, new, fixed = TRUE)[1L]
      if (cp > 0L)
        new <- paste0(substr(new, 1L, cp - 1L), oldcap,
                      substr(new, cp + nchar(newcap), nchar(new)))
    }
    tex <- paste0(substr(tex, 1L, first - 1L), new,
                  substr(tex, last + 1L, nchar(tex)))
    done <- c(which_label, done)
  }
  list(text = tex, updated = done)
}
run_replace <- function(path, replacement, target_name) {
  if (!file.exists(path)) {
    warning("Template not found: ", path, ". Table fragments still created.")
    return(NULL)
  }
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  out <- replace_tables(txt, replacement)
  outpath <- file.path(output_dir, target_name)
  writeLines(out$text, outpath, useBytes = TRUE)
  message("Updated ", length(out$updated), " LaTeX tables in ", outpath)
  message("  Labels: ", paste(out$updated, collapse = ", "))
  missing <- setdiff(names(replacement), out$updated)
  if (length(missing)) warning("No matching tables found in template: ", paste(missing, collapse = ", "))
  invisible(outpath)
}
main_updated <- run_replace(main_tex_file, main_tables, "Empirical_Main_UPDATED.tex")
appendix_updated <- run_replace(appendix_tex_file, appendix_tables,
                               "Empirical_Appendix_UPDATED.tex")

# Data provenance and change log
writeLines(c(
  "LGPIF TABLE UPDATE REPORT",
  paste("Original main:", main_tex_file),
  paste("Original appendix:", appendix_tex_file),
  paste("Results:", res_path),
  paste("Training years:", paste(train_years, collapse = ",")),
  "Models: 9 only; Case = 1 (XR1-3), 4c (ER1-3 continuous), 4d (ER1-3 discrete).",
  "OOS: frozen 2006-2009 fits evaluated in 2010, full vs exclusion of largest limited-loss policy-year.",
  "No numerical values have been copied from the old manuscripts.",
  "Original manuscript prose/figures and table captions are not rewritten; review them for numerical consistency.",
  "The table called marginal frequency GOF was commented out and is not replaced.",
  "Missing Godambe standard errors fall back to stage-wise SE and are labeled accordingly.",
  "Outlier-exclusion is post-hoc sensitivity analysis, not an alternative primary OOS result."
), file.path(output_dir, "README_TABLE_UPDATES.txt"))
cat("\nTABLE UPDATE COMPLETE. Output folder:", normalizePath(output_dir), "\n")




################################################################################
## Update severity distribution comparison table
## Uses existing empirical analysis output
################################################################################

results_dir <- "empirical_merged_outputs"

sev_tab <- readRDS(
  file.path(results_dir, "Severity_Family_Selection.rds")
)

# Distribution names and order
dist_names <- c(
  "GB2", "GP", "Burr", "IBurr", "Pareto",
  "Loglogi", "GG", "Weibull", "Gamma"
)

# Format numerical values
fmt <- function(x) {
  formatC(x, format = "f", digits = 1,
          big.mark = "", drop0trailing = FALSE)
}

# Generate LaTeX table rows
make_rows <- function(mat, title) {
  mat <- mat[c("Likelihood", "AIC", "BIC"),
             dist_names, drop = FALSE]
  
  c(
    paste0("\\textbf{", title, "} \\\\"),
    vapply(rownames(mat), function(metric) {
      paste0(
        metric, " & ",
        paste(fmt(mat[metric, ]), collapse = " & "),
        " \\\\"
      )
    }, character(1))
  )
}

latex <- c(
  "\\begin{table}[h!]",
  "\\renewcommand{\\tabcolsep}{2mm}",
  "\\renewcommand{\\arraystretch}{0.75}",
  "\\centering",
  "\\footnotesize",
  "\\begin{tabular}{lrrrrrrrrr}",
  "\\toprule",
  paste0(
    "& ",
    paste(paste0("\\textbf{", dist_names, "}"),
          collapse = " & "),
    " \\\\"
  ),
  "\\midrule",
  make_rows(sev_tab$complete, "Ground"),
  "\\midrule",
  make_rows(sev_tab$censored, "Censored"),
  "\\midrule",
  make_rows(sev_tab$truncated, "Truncated"),
  "\\bottomrule",
  "\\end{tabular}",
  paste0(
    "\\caption{In-sample Model Selection Criteria ",
    "for Ground, Censored, and Truncated Data}"
  ),
  "\\label{tab:model_comparison}",
  "\\end{table}"
)

# Export LaTeX file
output_tex <- file.path(
  results_dir,
  "Severity_Family_Model_Comparison_UPDATED.tex"
)

writeLines(latex, output_tex)

# Display in R console
cat(paste(latex, collapse = "\n"), "\n")

cat("\nSaved to:", output_tex, "\n")
