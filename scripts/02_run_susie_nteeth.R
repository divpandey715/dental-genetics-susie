############################
# Run SuSiE-RSS fine-mapping
# Trait: Nteeth
############################

rm(list = ls())

library(susieR)
library(data.table)
library(dplyr)

############################
# Paths
############################

root <- getwd()

gwas_file <- file.path(
  root,
  "data_samples/nteeth_summary.RDS"
)

ld_dir <- file.path(
  root,
  "LDMatrix/04_compute_LD_for_Shungin2019_EUR_Nteeth_output"
)

output_file <- file.path(
  root,
  "results/nteeth_susie_results.RDS"
)

############################
# Parameters
############################

WINDOW <- 500000
L_MAX <- 10
N_SAMPLE <- 26533

############################
# Load GWAS summary data
############################

gwas <- readRDS(gwas_file)

head(gwas)
dim(gwas)

gwas$MarkerName <- sub(":ID$", "", gwas$MarkerName)

gwas <- gwas %>%
  distinct(MarkerName, .keep_all = TRUE)

gwas <- as.data.table(gwas)

############################
# Parse chr / position
############################

gwas[, c("CHR", "BP") := tstrsplit(
  MarkerName,
  ":",
  fixed = TRUE,
  keep = 1:2
)]

gwas[, CHR := as.integer(CHR)]
gwas[, BP := as.integer(BP)]

############################
# List LD files
############################

ld_files <- list.files(
  ld_dir,
  pattern = "\\.ld$",
  full.names = TRUE
)

length(ld_files)

############################
# Container for SuSiE results
############################

susie_results <- vector("list", length(ld_files))
names(susie_results) <- basename(ld_files)

############################
# Run SuSiE-RSS by locus
############################

for (i in seq_along(ld_files)) {
  
  ############################
  # Load LD matrix
  ############################
  
  ld_file <- ld_files[i]
  
  locus_id <- sub(
    "^LD_|\\.ld$",
    "",
    basename(ld_file)
  )
  
  cat("SNP:", locus_id, "is running", "\n")
  
  ############################
  # Parse chr / position from file name
  ############################
  
  base <- sub("\\.ld$", "", locus_id)
  
  chr <- as.integer(sub("_.*", "", base))
  pos <- as.integer(sub(".*_", "", base))
  
  ############################
  # Read SNP list
  ############################
  
  snplist_file <- sub("\\.ld$", ".snplist", ld_file)
  
  if (!file.exists(snplist_file)) {
    warning("Missing .snplist file, skipping ", locus_id)
    next
  }
  
  snp_names <- readLines(snplist_file)
  
  # Keep only the first SNP name if multiple names are separated by ;
  snp_names <- sub(";.*", "", snp_names)
  
  ############################
  # Read LD square matrix
  ############################
  
  R_raw <- as.matrix(
    read.table(ld_file, header = FALSE)
  )
  
  if (nrow(R_raw) != length(snp_names)) {
    warning("LD matrix dimensions do not match SNP list, skipping ", locus_id)
    next
  }
  
  rownames(R_raw) <- colnames(R_raw) <- snp_names
  
  ############################
  # Map rsID to chr:pos coordinates
  ############################
  
  map_file <- sub("\\.ld$", ".id_map", ld_file)
  
  if (!file.exists(map_file)) {
    warning("Missing .id_map file, skipping ", locus_id)
    next
  }
  
  id_lookup <- read.table(
    map_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(id_lookup) <- c("rsid", "coords")
  
  matched_indices <- match(
    snp_names,
    id_lookup$rsid
  )
  
  chr_pos_names <- id_lookup$coords[matched_indices]
  
  rownames(R_raw) <- colnames(R_raw) <- chr_pos_names
  
  snps_in_ld <- colnames(R_raw)
  
  ############################
  # Align GWAS and LD
  ############################
  
  gwas_locus <- gwas[
    CHR == chr &
      BP >= (pos - WINDOW) &
      BP <= (pos + WINDOW)
  ]
  
  if (nrow(gwas_locus) < 10) {
    warning("Too few GWAS SNPs, skipping ", locus_id)
    next
  }
  
  common_snps <- intersect(
    gwas_locus$MarkerName,
    snps_in_ld
  )
  
  if (length(common_snps) < 10) {
    warning("Too few overlapping SNPs, skipping ", locus_id)
    next
  }
  
  ############################
  # Subset and reorder
  ############################
  
  gwas_locus <- gwas_locus[
    MarkerName %in% common_snps
  ]
  
  setkey(gwas_locus, MarkerName)
  
  R <- R_raw[
    common_snps,
    common_snps
  ]
  
  ############################
  # Remove SNPs with missing LD
  ############################
  
  na_snps <- colSums(is.na(R)) > 0 |
    rowSums(is.na(R)) > 0
  
  R <- R[
    !na_snps,
    !na_snps
  ]
  
  if (nrow(R) < 10) {
    warning("Too few SNPs after removing NA LD rows/columns, skipping ", locus_id)
    next
  }
  
  ############################
  # Clean LD matrix
  ############################
  
  R <- (R + t(R)) / 2
  diag(R) <- 1
  
  ############################
  # Final GWAS / LD alignment
  ############################
  
  gwas_locus <- gwas_locus[
    colnames(R)
  ]
  
  stopifnot(
    identical(
      gwas_locus$MarkerName,
      colnames(R)
    )
  )
  
  ############################
  # Run SuSiE-RSS
  ############################
  
  z <- gwas_locus$Effect / gwas_locus$StdErr
  
  fit <- susie_rss(
    z = z,
    R = R,
    n = N_SAMPLE,
    L = L_MAX,
    estimate_prior_method = "EM",
    estimate_residual_variance = FALSE
  )
  
  susie_results[[i]] <- list(
    locus = paste0("chr", chr, ":", pos),
    fit = fit,
    gwas = gwas_locus
  )
}

############################
# Save SuSiE results
############################

saveRDS(
  susie_results,
  file = output_file
)

cat("SuSiE results saved to:", output_file, "\n")