{\rtf1\ansi\ansicpg1252\cocoartf2820
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Roman;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs24 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 # ============================================\
# SuSiE-RSS fine-mapping pipeline\
# Traits: Nteeth, DMFS, Perio\
# ============================================\
\
rm(list = ls())\
\
library(susieR)\
library(data.table)\
library(dplyr)\
\
# Use the current working directory as the project root.\
# Before running, set working directory to:\
# ~/Documents/GitHub/dental-genetics-susie\
root <- getwd()\
\
traits <- list(\
  list(\
    trait_name = "nteeth",\
    gwas_file = file.path(root, "data_samples/nteeth_summary.RDS"),\
    ld_dir = file.path(root, "LDMatrix/04_compute_LD_for_Shungin2019_EUR_Nteeth_output"),\
    output_file = file.path(root, "results/nteeth_susie_results.RDS"),\
    summary_file = file.path(root, "results/nteeth_susie_summary.csv"),\
    n_sample = 26533\
  ),\
  list(\
    trait_name = "dmfs",\
    gwas_file = file.path(root, "data_samples/dmfs_summary.RDS"),\
    ld_dir = file.path(root, "LDMatrix/04_compute_LD_for_Shungin2019_EUR_DMFS_output"),\
    output_file = file.path(root, "results/dmfs_susie_results.RDS"),\
    summary_file = file.path(root, "results/dmfs_susie_summary.csv"),\
    n_sample = 26533\
  ),\
  list(\
    trait_name = "perio",\
    gwas_file = file.path(root, "data_samples/perio_summary.RDS"),\
    ld_dir = file.path(root, "LDMatrix/04_compute_LD_for_Shungin2019_EUR_Perio_output"),\
    output_file = file.path(root, "results/perio_susie_results.RDS"),\
    summary_file = file.path(root, "results/perio_susie_summary.csv"),\
    n_sample = 26533\
  )\
)\
\
WINDOW <- 500000\
L_MAX <- 10\
\
dir.create(file.path(root, "results"), showWarnings = FALSE, recursive = TRUE)\
\
run_susie_trait <- function(trait_name, gwas_file, ld_dir, output_file, summary_file, n_sample) \{\
  \
  cat("\\n=============================\\n")\
  cat("Running trait:", trait_name, "\\n")\
  cat("=============================\\n")\
  \
  if (!file.exists(gwas_file)) \{\
    warning("GWAS file not found: ", gwas_file)\
    return(NULL)\
  \}\
  \
  if (!dir.exists(ld_dir)) \{\
    warning("LD directory not found: ", ld_dir)\
    return(NULL)\
  \}\
  \
  gwas <- readRDS(gwas_file)\
  \
  gwas$MarkerName <- sub(":ID$", "", gwas$MarkerName)\
  gwas <- gwas %>% distinct(MarkerName, .keep_all = TRUE)\
  gwas <- as.data.table(gwas)\
  \
  gwas[, c("CHR", "BP") := tstrsplit(MarkerName, ":", fixed = TRUE, keep = 1:2)]\
  gwas[, CHR := as.integer(CHR)]\
  gwas[, BP := as.integer(BP)]\
  \
  ld_files <- list.files(ld_dir, pattern = "\\\\.ld$", full.names = TRUE)\
  \
  if (length(ld_files) == 0) \{\
    warning("No LD files found for: ", trait_name)\
    return(NULL)\
  \}\
  \
  susie_results <- vector("list", length(ld_files))\
  names(susie_results) <- basename(ld_files)\
  \
  for (i in seq_along(ld_files)) \{\
    \
    ld_file <- ld_files[i]\
    locus_id <- sub("^LD_|\\\\.ld$", "", basename(ld_file))\
    \
    cat("Running locus:", locus_id, "\\n")\
    \
    base <- sub("\\\\.ld$", "", locus_id)\
    chr <- as.integer(sub("_.*", "", base))\
    pos <- as.integer(sub(".*_", "", base))\
    \
    snplist_file <- sub("\\\\.ld$", ".snplist", ld_file)\
    map_file <- sub("\\\\.ld$", ".id_map", ld_file)\
    \
    if (!file.exists(snplist_file) || !file.exists(map_file)) \{\
      warning("Missing snplist or id_map file, skipping: ", locus_id)\
      next\
    \}\
    \
    snp_names <- readLines(snplist_file)\
    snp_names <- sub(";.*", "", snp_names)\
    \
    R_raw <- as.matrix(read.table(ld_file, header = FALSE))\
    \
    if (nrow(R_raw) != length(snp_names)) \{\
      warning("LD matrix size mismatch, skipping: ", locus_id)\
      next\
    \}\
    \
    rownames(R_raw) <- colnames(R_raw) <- snp_names\
    \
    id_lookup <- read.table(map_file, header = FALSE, stringsAsFactors = FALSE)\
    colnames(id_lookup) <- c("rsid", "coords")\
    \
    matched_indices <- match(snp_names, id_lookup$rsid)\
    chr_pos_names <- id_lookup$coords[matched_indices]\
    \
    rownames(R_raw) <- colnames(R_raw) <- chr_pos_names\
    \
    snps_in_ld <- colnames(R_raw)\
    \
    gwas_locus <- gwas[\
      CHR == chr &\
        BP >= (pos - WINDOW) &\
        BP <= (pos + WINDOW)\
    ]\
    \
    if (nrow(gwas_locus) < 10) \{\
      warning("Too few GWAS SNPs, skipping: ", locus_id)\
      next\
    \}\
    \
    common_snps <- intersect(gwas_locus$MarkerName, snps_in_ld)\
    \
    if (length(common_snps) < 10) \{\
      warning("Too few overlapping SNPs, skipping: ", locus_id)\
      next\
    \}\
    \
    gwas_locus <- gwas_locus[MarkerName %in% common_snps]\
    setkey(gwas_locus, MarkerName)\
    \
    R <- R_raw[common_snps, common_snps]\
    \
    na_snps <- colSums(is.na(R)) > 0 | rowSums(is.na(R)) > 0\
    R <- R[!na_snps, !na_snps]\
    \
    if (nrow(R) < 10) \{\
      warning("Too few SNPs after NA removal, skipping: ", locus_id)\
      next\
    \}\
    \
    R <- (R + t(R)) / 2\
    diag(R) <- 1\
    \
    gwas_locus <- gwas_locus[colnames(R)]\
    \
    if (!identical(gwas_locus$MarkerName, colnames(R))) \{\
      warning("GWAS and LD SNP mismatch, skipping: ", locus_id)\
      next\
    \}\
    \
    z <- gwas_locus$Effect / gwas_locus$StdErr\
    \
    fit <- susie_rss(\
      z = z,\
      R = R,\
      n = n_sample,\
      L = L_MAX,\
      estimate_prior_method = "EM",\
      estimate_residual_variance = FALSE\
    )\
    \
    susie_results[[i]] <- list(\
      trait = trait_name,\
      locus = paste0("chr", chr, ":", pos),\
      fit = fit,\
      gwas = gwas_locus\
    )\
  \}\
  \
  saveRDS(susie_results, file = output_file)\
  \
  summary_list <- list()\
  \
  for (i in seq_along(susie_results)) \{\
    \
    res <- susie_results[[i]]\
    if (is.null(res)) next\
    \
    fit <- res$fit\
    gwas_data <- res$gwas\
    \
    if (!is.null(fit$sets$cs)) \{\
      \
      for (cs_id in names(fit$sets$cs)) \{\
        \
        snp_idx <- fit$sets$cs[[cs_id]]\
        cs_snps <- gwas_data[snp_idx, ]\
        \
        cs_snps$Trait <- res$trait\
        cs_snps$PIP <- fit$pip[snp_idx]\
        cs_snps$CS <- cs_id\
        cs_snps$Locus <- res$locus\
        \
        if (!is.null(fit$sets$purity)) \{\
          cs_snps$Purity <- fit$sets$purity[cs_id, "min.abs.corr"]\
        \} else \{\
          cs_snps$Purity <- NA\
        \}\
        \
        summary_list[[length(summary_list) + 1]] <- cs_snps\
      \}\
    \}\
  \}\
  \
  if (length(summary_list) > 0) \{\
    susie_summary <- bind_rows(summary_list)\
    susie_summary <- susie_summary %>% arrange(desc(PIP))\
    write.csv(susie_summary, file = summary_file, row.names = FALSE)\
  \}\
  \
  cat("Finished trait:", trait_name, "\\n")\
\}\
\
for (trait in traits) \{\
  do.call(run_susie_trait, trait)\
\}\
\
cat("\\nAll traits complete.\\n")}