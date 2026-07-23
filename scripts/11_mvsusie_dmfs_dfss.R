############################
# mvSuSiE multi-trait fine-mapping
# Traits: DMFS + DFSS
# Test run on selected strongest loci
############################

rm(list = ls())

library(data.table)
library(dplyr)
library(mvsusieR)

############################
# Paths
############################

root <- "/Users/divyapandey/Documents/GitHub/dental-genetics-susie"

dmfs_summary_file <- file.path(root, "data_samples/dmfs_summary.RDS")
dfss_summary_file <- file.path(root, "data_samples/dfss_summary.RDS")

dmfs_ld_dir <- file.path(root, "LDMatrix/04_compute_LD_for_Shungin2019_EUR_DMFS_output")
dfss_ld_dir <- file.path(root, "LDMatrix/04_compute_LD_for_Shungin2019_EUR_DFSS_output")

results_dir <- file.path(root, "results")
dir.create(results_dir, showWarnings = FALSE)

output_rds <- file.path(results_dir, "dmfs_dfss_mvsusie_results.RDS")
output_csv <- file.path(results_dir, "dmfs_dfss_mvsusie_summary.csv")
output_log <- file.path(results_dir, "dmfs_dfss_mvsusie_log.txt")

############################
# mvSuSiE prior variance
############################

mvsusie_prior <- matrix(
  c(
    0.2, 0.1,
    0.1, 0.2
  ),
  nrow = 2,
  byrow = TRUE
)

############################
# Load GWAS summary statistics
############################

cat("Loading DMFS summary statistics...\n")
dmfs <- readRDS(dmfs_summary_file)
setDT(dmfs)

cat("Loading DFSS summary statistics...\n")
dfss <- readRDS(dfss_summary_file)
setDT(dfss)

setkey(dmfs, MarkerName)
setkey(dfss, MarkerName)

############################
# Helper function: run one locus
############################

run_mvsusie_locus <- function(prefix, ld_dir) {
  
  cat("\nRunning locus:", prefix, "\n")
  
  ld_file <- file.path(ld_dir, paste0(prefix, ".ld"))
  snplist_file <- file.path(ld_dir, paste0(prefix, ".snplist"))
  id_map_file <- file.path(ld_dir, paste0(prefix, ".id_map"))
  
  if (!file.exists(ld_file) || !file.exists(snplist_file) || !file.exists(id_map_file)) {
    cat("Skipping", prefix, "- missing LD, snplist, or id_map file\n")
    return(NULL)
  }
  
  ############################
  # Load snplist and id_map
  ############################
  
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
  
  if (nrow(snp_map) < 2) {
    cat("Skipping", prefix, "- fewer than 2 SNPs mapped from rsID to MarkerName\n")
    return(NULL)
  }
  
  ############################
  # Keep SNPs present in both GWAS files
  ############################
  
  snp_map <- snp_map[
    MarkerName %in% dmfs$MarkerName &
      MarkerName %in% dfss$MarkerName
  ]
  
  snp_map <- snp_map[order(ld_index)]
  
  if (nrow(snp_map) < 2) {
    cat("Skipping", prefix, "- fewer than 2 SNPs shared across LD, DMFS, and DFSS\n")
    return(NULL)
  }
  
  ############################
  # Load and subset LD matrix
  ############################
  
  R <- as.matrix(fread(ld_file, header = FALSE))
  storage.mode(R) <- "numeric"
  
  if (nrow(R) != length(snps_rsid) || ncol(R) != length(snps_rsid)) {
    cat("Skipping", prefix, "- LD dimensions do not match snplist\n")
    return(NULL)
  }
  
  R <- R[snp_map$ld_index, snp_map$ld_index, drop = FALSE]
  
  ############################
  # Remove SNPs with missing LD values
  ############################
  
  finite_ld <- rowSums(!is.finite(R)) == 0 & colSums(!is.finite(R)) == 0
  
  if (sum(finite_ld) < 2) {
    cat("Skipping", prefix, "- fewer than 2 SNPs after removing missing LD values\n")
    return(NULL)
  }
  
  if (sum(!finite_ld) > 0) {
    cat("Removing", sum(!finite_ld), "SNPs with missing LD values from", prefix, "\n")
  }
  
  R <- R[finite_ld, finite_ld, drop = FALSE]
  snp_map <- snp_map[finite_ld]
  
  R <- (R + t(R)) / 2
  diag(R) <- 1
  
  ############################
  # Align GWAS data
  ############################
  
  dmfs_locus <- dmfs[J(snp_map$MarkerName)]
  dfss_locus <- dfss[J(snp_map$MarkerName)]
  
  if (any(is.na(dmfs_locus$Effect)) || any(is.na(dfss_locus$Effect))) {
    cat("Skipping", prefix, "- missing effect sizes after alignment\n")
    return(NULL)
  }
  
  ############################
  # Build mvSuSiE inputs
  ############################
  
  Z <- cbind(
    DMFS = dmfs_locus$Effect / dmfs_locus$StdErr,
    DFSS = dfss_locus$Effect / dfss_locus$StdErr
  )
  
  if (any(!is.finite(Z))) {
    cat("Skipping", prefix, "- non-finite Z-scores\n")
    return(NULL)
  }
  
  N_vec <- c(
    DMFS = median(dmfs_locus$N, na.rm = TRUE),
    DFSS = median(dfss_locus$N, na.rm = TRUE)
  )
  
  ############################
  # Run mvSuSiE-RSS
  ############################
  
  fit <- tryCatch(
    {
      mvsusie_rss(
        Z = Z,
        R = R,
        N = median(N_vec, na.rm = TRUE),
        prior_variance = mvsusie_prior,
        L = 10,
        max_iter = 1000,
        tol = 1e-04
      )
    },
    error = function(e) {
      cat("Error in", prefix, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  ############################
  # Extract joint mvSuSiE PIPs
  ############################
  
  overall_pip <- tryCatch(
    {
      if (!is.null(fit$pip)) {
        
        pip <- as.numeric(fit$pip)
        
      } else if (!is.null(fit$alpha)) {
        
        pip <- 1 - apply(
          1 - fit$alpha,
          2,
          prod
        )
        
        pip <- as.numeric(pip)
        
      } else {
        
        pip <- rep(
          NA_real_,
          nrow(snp_map)
        )
      }
      
      if (length(pip) != nrow(snp_map)) {
        warning(
          "PIP length does not match number of SNPs in ",
          prefix,
          "; setting PIP to NA."
        )
        
        pip <- rep(
          NA_real_,
          nrow(snp_map)
        )
      }
      
      pip[pip < 0] <- 0
      pip[pip > 1] <- 1
      
      pip
    },
    error = function(e) {
      rep(
        NA_real_,
        nrow(snp_map)
      )
    }
  )
  
  if (length(overall_pip) != nrow(snp_map)) {
    overall_pip <- rep(NA_real_, nrow(snp_map))
  }
  
  ############################
  # Credible set annotation
  ############################
  
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
  
  ############################
  # Output locus table
  ############################
  
  locus_summary <- data.frame(
    Locus = prefix,
    rsid = snp_map$rsid,
    MarkerName = snp_map$MarkerName,
    CS = cs_index,
    PIP = overall_pip,
    Purity = purity,
    Z_DMFS = Z[, "DMFS"],
    Z_DFSS = Z[, "DFSS"],
    Effect_DMFS = dmfs_locus$Effect,
    StdErr_DMFS = dmfs_locus$StdErr,
    Pvalue_DMFS = dmfs_locus$P.value,
    N_DMFS = dmfs_locus$N,
    Effect_DFSS = dfss_locus$Effect,
    StdErr_DFSS = dfss_locus$StdErr,
    Pvalue_DFSS = dfss_locus$`P-value`,
    N_DFSS = dfss_locus$N
  )
  
  cat("Completed", prefix, "- max PIP:", max(overall_pip, na.rm = TRUE), "\n")
  
  return(
    list(
      locus = prefix,
      fit = fit,
      summary = locus_summary
    )
  )
}

############################
# Identify LD loci
############################

dmfs_ld_files <- list.files(dmfs_ld_dir, pattern = "\\.ld$", full.names = FALSE)
dfss_ld_files <- list.files(dfss_ld_dir, pattern = "\\.ld$", full.names = FALSE)

dmfs_prefixes <- sub("\\.ld$", "", dmfs_ld_files)
dfss_prefixes <- sub("\\.ld$", "", dfss_ld_files)

############################
# TEST RUN: selected strongest loci only
############################

all_prefixes <- c(
  "LD_13_32370438",
  "LD_5_134504407",
  "LD_5_134507139",
  "LD_5_134507859",
  "LD_5_134508559",
  "LD_5_134509677",
  "LD_5_134509987",
  "LD_5_134510772"
)

cat("Number of LD loci selected:", length(all_prefixes), "\n")

############################
# Run mvSuSiE across selected loci
############################

all_results <- list()
all_summaries <- list()
log_messages <- c()

for (prefix in all_prefixes) {
  
  if (prefix %in% dmfs_prefixes) {
    ld_dir <- dmfs_ld_dir
  } else if (prefix %in% dfss_prefixes) {
    ld_dir <- dfss_ld_dir
  } else {
    cat("Skipping", prefix, "- not found in either LD folder\n")
    log_messages <- c(log_messages, paste("Skipped:", prefix, "- not found"))
    next
  }
  
  res <- run_mvsusie_locus(prefix, ld_dir)
  
  if (!is.null(res)) {
    all_results[[prefix]] <- res$fit
    all_summaries[[prefix]] <- res$summary
    log_messages <- c(log_messages, paste("Completed:", prefix))
  } else {
    log_messages <- c(log_messages, paste("Skipped:", prefix))
  }
}

############################
# Save results
############################

saveRDS(all_results, output_rds)

if (length(all_summaries) > 0) {
  
  final_summary <- bind_rows(all_summaries)
  
  final_summary <- final_summary %>%
    arrange(desc(PIP), desc(abs(Z_DMFS)), desc(abs(Z_DFSS)))
  
  write.csv(
    final_summary,
    output_csv,
    row.names = FALSE
  )
  
  cat("\nSaved mvSuSiE summary to:\n")
  cat(output_csv, "\n")
  
} else {
  final_summary <- data.frame()
  cat("\nNo mvSuSiE loci completed successfully.\n")
}

writeLines(log_messages, output_log)

cat("\nSaved mvSuSiE model objects to:\n")
cat(output_rds, "\n")

cat("\nSaved log to:\n")
cat(output_log, "\n")

cat("\nNumber of completed loci:", length(all_results), "\n")
cat("Number of SNP rows in final summary:", nrow(final_summary), "\n")

if (nrow(final_summary) > 0) {
  cat("\nTop mvSuSiE results:\n")
  print(
    final_summary %>%
      arrange(desc(PIP)) %>%
      head(20)
  )
}
