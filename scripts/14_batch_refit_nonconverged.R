############################
# Batch re-fit of non-converged mvSuSiE loci
# Target: the 9 locus/trait-pair fits (7 unique loci, all
# chr5:~134.50-134.51Mb) that hit max_iter=1000 without
# converging in scripts 07, 11, 12. Re-fits with a higher
# max_iter and patches the existing results/summary files
# in place, leaving every other locus untouched.
############################

rm(list = ls())

library(data.table)
library(dplyr)
library(mvsusieR)

root <- getwd()

dmfs_summary_file <- file.path(root, "data_samples/dmfs_summary.RDS")
nteeth_summary_file <- file.path(root, "data_samples/nteeth_summary.RDS")
dfss_summary_file <- file.path(root, "data_samples/dfss_summary.RDS")

dmfs_ld_dir <- file.path(root, "LDMatrix/04_compute_LD_for_Shungin2019_EUR_DMFS_output")
nteeth_ld_dir <- file.path(root, "LDMatrix/04_compute_LD_for_Shungin2019_EUR_Nteeth_output")
dfss_ld_dir <- file.path(root, "LDMatrix/04_compute_LD_for_Shungin2019_EUR_DFSS_output")

results_dir <- file.path(root, "results")

mvsusie_prior <- matrix(
  c(
    0.2, 0.1,
    0.1, 0.2
  ),
  nrow = 2,
  byrow = TRUE
)

REFIT_MAX_ITER <- 3000

cat("Loading DMFS summary statistics...\n")
dmfs <- readRDS(dmfs_summary_file)
setDT(dmfs)
setkey(dmfs, MarkerName)

cat("Loading Nteeth summary statistics...\n")
nteeth <- readRDS(nteeth_summary_file)
setDT(nteeth)
setkey(nteeth, MarkerName)

cat("Loading DFSS summary statistics...\n")
dfss <- readRDS(dfss_summary_file)
setDT(dfss)
setkey(dfss, MarkerName)

############################
# Shared locus re-fit helper
############################

refit_locus <- function(prefix, ld_dir, trait_a_name, trait_a_data,
                         trait_b_name, trait_b_data) {

  cat("\nRe-fitting locus:", prefix, "(", trait_a_name, "+", trait_b_name, ")\n")

  ld_file <- file.path(ld_dir, paste0(prefix, ".ld"))
  snplist_file <- file.path(ld_dir, paste0(prefix, ".snplist"))
  id_map_file <- file.path(ld_dir, paste0(prefix, ".id_map"))

  snps_rsid <- fread(snplist_file, header = FALSE)[[1]]
  snps_rsid <- as.character(snps_rsid)

  id_map <- fread(id_map_file, header = FALSE)
  setnames(id_map, c("rsid", "MarkerName"))

  snp_map <- data.table(
    rsid = snps_rsid,
    ld_index = seq_along(snps_rsid)
  )

  snp_map <- merge(
    snp_map,
    id_map,
    by = "rsid",
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )

  snp_map <- snp_map[!is.na(MarkerName)]

  snp_map <- snp_map[
    MarkerName %in% trait_a_data$MarkerName &
      MarkerName %in% trait_b_data$MarkerName
  ]

  snp_map <- snp_map[order(ld_index)]

  R <- as.matrix(fread(ld_file, header = FALSE))
  storage.mode(R) <- "numeric"

  R <- R[snp_map$ld_index, snp_map$ld_index, drop = FALSE]

  finite_ld <- rowSums(!is.finite(R)) == 0 & colSums(!is.finite(R)) == 0

  R <- R[finite_ld, finite_ld, drop = FALSE]
  snp_map <- snp_map[finite_ld]

  R <- (R + t(R)) / 2
  diag(R) <- 1

  ev_min <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)

  if (ev_min < 1e-6) {
    lambda <- min(0.5, ((1e-6 - ev_min) / (1 - ev_min)) * 1.05)
    R <- (1 - lambda) * R + lambda * diag(nrow(R))
    diag(R) <- 1
    cat(
      "Regularized non-PSD LD matrix for", prefix,
      "- min eig was", round(ev_min, 4),
      ", shrinkage lambda =", round(lambda, 4), "\n"
    )
  }

  a_locus <- trait_a_data[J(snp_map$MarkerName)]
  b_locus <- trait_b_data[J(snp_map$MarkerName)]

  Z <- cbind(
    A = a_locus$Effect / a_locus$StdErr,
    B = b_locus$Effect / b_locus$StdErr
  )
  colnames(Z) <- c(trait_a_name, trait_b_name)

  N_vec <- c(
    A = median(a_locus$N, na.rm = TRUE),
    B = median(b_locus$N, na.rm = TRUE)
  )

  fit <- tryCatch(
    {
      mvsusie_rss(
        Z = Z,
        R = R,
        N = median(N_vec, na.rm = TRUE),
        prior_variance = mvsusie_prior,
        L = 10,
        max_iter = REFIT_MAX_ITER,
        tol = 1e-04
      )
    },
    error = function(e) {
      cat("Error re-fitting", prefix, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(fit)) {
    return(NULL)
  }

  cat(
    "Result for", prefix, "- converged:", fit$converged,
    "| niter:", fit$niter, "\n"
  )

  overall_pip <- tryCatch(
    {
      if (!is.null(fit$pip)) {
        pip <- as.numeric(fit$pip)
      } else if (!is.null(fit$alpha)) {
        pip <- 1 - apply(1 - fit$alpha, 2, prod)
        pip <- as.numeric(pip)
      } else {
        pip <- rep(NA_real_, nrow(snp_map))
      }

      if (length(pip) != nrow(snp_map)) {
        pip <- rep(NA_real_, nrow(snp_map))
      }

      pip[pip < 0] <- 0
      pip[pip > 1] <- 1

      pip
    },
    error = function(e) rep(NA_real_, nrow(snp_map))
  )

  cs_index <- rep(NA_integer_, nrow(snp_map))

  if (!is.null(fit$sets$cs)) {
    for (i in seq_along(fit$sets$cs)) {
      cs_index[fit$sets$cs[[i]]] <- i
    }
  }

  purity <- rep(NA_real_, nrow(snp_map))

  if (!is.null(fit$sets$purity)) {
    for (i in seq_along(fit$sets$cs)) {
      if ("min.abs.corr" %in% colnames(fit$sets$purity)) {
        purity[fit$sets$cs[[i]]] <- fit$sets$purity[i, "min.abs.corr"]
      }
    }
  }

  list(
    fit = fit,
    snp_map = snp_map,
    Z = Z,
    a_locus = a_locus,
    b_locus = b_locus,
    overall_pip = overall_pip,
    cs_index = cs_index,
    purity = purity
  )
}

############################
# Patch a saved mvSuSiE results/summary pair for one script
############################

patch_results <- function(results_rds, summary_csv, prefix, refit,
                           trait_a_name, trait_b_name,
                           trait_a_pcol, trait_b_pcol) {

  all_results <- readRDS(results_rds)

  all_results[[prefix]] <- refit$fit

  saveRDS(all_results, results_rds)

  final_summary <- read.csv(summary_csv, stringsAsFactors = FALSE)

  final_summary <- final_summary[final_summary$Locus != prefix, ]

  locus_summary <- data.frame(
    Locus = prefix,
    rsid = refit$snp_map$rsid,
    MarkerName = refit$snp_map$MarkerName,
    CS = refit$cs_index,
    PIP = refit$overall_pip,
    Purity = refit$purity,
    stringsAsFactors = FALSE
  )

  locus_summary[[paste0("Z_", trait_a_name)]] <- refit$Z[, trait_a_name]
  locus_summary[[paste0("Z_", trait_b_name)]] <- refit$Z[, trait_b_name]
  locus_summary[[paste0("Effect_", trait_a_name)]] <- refit$a_locus$Effect
  locus_summary[[paste0("StdErr_", trait_a_name)]] <- refit$a_locus$StdErr
  locus_summary[[paste0("Pvalue_", trait_a_name)]] <- refit$a_locus[[trait_a_pcol]]
  locus_summary[[paste0("N_", trait_a_name)]] <- refit$a_locus$N
  locus_summary[[paste0("Effect_", trait_b_name)]] <- refit$b_locus$Effect
  locus_summary[[paste0("StdErr_", trait_b_name)]] <- refit$b_locus$StdErr
  locus_summary[[paste0("Pvalue_", trait_b_name)]] <- refit$b_locus[[trait_b_pcol]]
  locus_summary[[paste0("N_", trait_b_name)]] <- refit$b_locus$N

  locus_summary <- locus_summary[names(final_summary)]

  final_summary <- rbind(final_summary, locus_summary)

  final_summary <- final_summary[order(-final_summary$PIP), ]

  write.csv(final_summary, summary_csv, row.names = FALSE)

  cat("Patched", results_rds, "and", summary_csv, "for", prefix, "\n")
}

############################
# Job list: 9 (locus, trait-pair) re-fits
############################

dmfs_pcol <- "P.value"
nteeth_pcol <- "P.value"
dfss_pcol <- "P-value"

jobs <- list(
  list(prefix = "LD_5_134504407", ld_dir = dmfs_ld_dir,
       results_rds = file.path(results_dir, "dmfs_nteeth_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "dmfs_nteeth_mvsusie_summary.csv"),
       trait_a_name = "DMFS", trait_a_data = dmfs, trait_a_pcol = dmfs_pcol,
       trait_b_name = "Nteeth", trait_b_data = nteeth, trait_b_pcol = nteeth_pcol),

  list(prefix = "LD_5_134504407", ld_dir = dmfs_ld_dir,
       results_rds = file.path(results_dir, "dmfs_dfss_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "dmfs_dfss_mvsusie_summary.csv"),
       trait_a_name = "DMFS", trait_a_data = dmfs, trait_a_pcol = dmfs_pcol,
       trait_b_name = "DFSS", trait_b_data = dfss, trait_b_pcol = dfss_pcol),
  list(prefix = "LD_5_134507139", ld_dir = dmfs_ld_dir,
       results_rds = file.path(results_dir, "dmfs_dfss_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "dmfs_dfss_mvsusie_summary.csv"),
       trait_a_name = "DMFS", trait_a_data = dmfs, trait_a_pcol = dmfs_pcol,
       trait_b_name = "DFSS", trait_b_data = dfss, trait_b_pcol = dfss_pcol),
  list(prefix = "LD_5_134507859", ld_dir = dmfs_ld_dir,
       results_rds = file.path(results_dir, "dmfs_dfss_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "dmfs_dfss_mvsusie_summary.csv"),
       trait_a_name = "DMFS", trait_a_data = dmfs, trait_a_pcol = dmfs_pcol,
       trait_b_name = "DFSS", trait_b_data = dfss, trait_b_pcol = dfss_pcol),
  list(prefix = "LD_5_134508559", ld_dir = dmfs_ld_dir,
       results_rds = file.path(results_dir, "dmfs_dfss_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "dmfs_dfss_mvsusie_summary.csv"),
       trait_a_name = "DMFS", trait_a_data = dmfs, trait_a_pcol = dmfs_pcol,
       trait_b_name = "DFSS", trait_b_data = dfss, trait_b_pcol = dfss_pcol),
  list(prefix = "LD_5_134509677", ld_dir = dmfs_ld_dir,
       results_rds = file.path(results_dir, "dmfs_dfss_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "dmfs_dfss_mvsusie_summary.csv"),
       trait_a_name = "DMFS", trait_a_data = dmfs, trait_a_pcol = dmfs_pcol,
       trait_b_name = "DFSS", trait_b_data = dfss, trait_b_pcol = dfss_pcol),
  list(prefix = "LD_5_134509987", ld_dir = dmfs_ld_dir,
       results_rds = file.path(results_dir, "dmfs_dfss_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "dmfs_dfss_mvsusie_summary.csv"),
       trait_a_name = "DMFS", trait_a_data = dmfs, trait_a_pcol = dmfs_pcol,
       trait_b_name = "DFSS", trait_b_data = dfss, trait_b_pcol = dfss_pcol),
  list(prefix = "LD_5_134510772", ld_dir = dmfs_ld_dir,
       results_rds = file.path(results_dir, "dmfs_dfss_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "dmfs_dfss_mvsusie_summary.csv"),
       trait_a_name = "DMFS", trait_a_data = dmfs, trait_a_pcol = dmfs_pcol,
       trait_b_name = "DFSS", trait_b_data = dfss, trait_b_pcol = dfss_pcol),

  list(prefix = "LD_5_134504407", ld_dir = dfss_ld_dir,
       results_rds = file.path(results_dir, "nteeth_dfss_mvsusie_results.RDS"),
       summary_csv = file.path(results_dir, "nteeth_dfss_mvsusie_summary.csv"),
       trait_a_name = "Nteeth", trait_a_data = nteeth, trait_a_pcol = nteeth_pcol,
       trait_b_name = "DFSS", trait_b_data = dfss, trait_b_pcol = dfss_pcol)
)

cat("\nTotal re-fit jobs:", length(jobs), "\n")

for (job in jobs) {

  refit <- refit_locus(
    prefix = job$prefix,
    ld_dir = job$ld_dir,
    trait_a_name = job$trait_a_name,
    trait_a_data = job$trait_a_data,
    trait_b_name = job$trait_b_name,
    trait_b_data = job$trait_b_data
  )

  if (is.null(refit)) {
    cat("Skipping patch for", job$prefix, "- re-fit failed\n")
    next
  }

  patch_results(
    results_rds = job$results_rds,
    summary_csv = job$summary_csv,
    prefix = job$prefix,
    refit = refit,
    trait_a_name = job$trait_a_name,
    trait_b_name = job$trait_b_name,
    trait_a_pcol = job$trait_a_pcol,
    trait_b_pcol = job$trait_b_pcol
  )
}

cat("\nBatch re-fit complete.\n")
