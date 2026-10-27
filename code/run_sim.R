## Driver:  Rscript run_sim.R <max_rep> [cases] [outdir]
##   Replications are processed in the order (rep 1: all cases), (rep 2: all cases), ... so that every case accumulates
##   at the same rate. Any number of workers can run this script concurrently on the same outdir (one per core):
##   a replication is claimed by creating a lock directory (atomic), finished replications are skipped.
##   e.g. on an 8-core machine:   for w in 1 2 3 4 5 6 7 8; do nohup Rscript run_sim.R 500 > log$w.txt 2>&1 & done
args <- commandArgs(trailingOnly = TRUE)
max_rep <- as.integer(args[1])
cases   <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else 1:4
outdir  <- if (length(args) >= 3) args[3] else "results"
source("sim_functions.R")
dir.create(outdir, showWarnings = FALSE); lockdir <- file.path(outdir, "locks"); dir.create(lockdir, showWarnings = FALSE)
grid <- expand.grid(case = cases, rep = seq_len(max_rep))            # case varies fastest
for (i in seq_len(nrow(grid))) {
  cs <- grid$case[i]; rp <- grid$rep[i]
  if (file.exists(file.path(outdir, sprintf("case%s_rep%04d.rds", cs, rp)))) next
  if (!dir.create(file.path(lockdir, sprintf("case%s_rep%04d", cs, rp)), showWarnings = FALSE)) next   # claimed by another worker
  sw <- if (cs == 4) c("FEM2", "SDM3", "FDM2") else "FEM2"          # Godambe s.e. for the coverage study (full set in Case 4)
  res <- tryCatch(run_replication(cs, rp, outdir = outdir, sandwich = sw),
                  error = function(e) { cat(sprintf("case %d rep %d FAILED: %s\n", cs, rp, conditionMessage(e))); NULL })
  if (!is.null(res)) cat(sprintf("[pid %d] case %d rep %4d done in %6.1f s  (claims %d, conv %s)\n", Sys.getpid(), cs, rp, res$time_total,
                                 res$n_claims, paste(res$convergence, collapse = "")))
  flush.console()
}
