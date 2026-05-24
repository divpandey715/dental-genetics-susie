############################
# Summarize SuSiE-RSS results
# Trait: Nteeth
############################

rm(list = ls())

library(data.table)
library(dplyr)

############################
# Paths
############################

root <- getwd()

susie_results_file <- file.path(
  root,
  "results/nteeth_susie_results.RDS"
)

summary_output_file <- file.path(
  root,
  "results/nteeth_susie_final_summary.csv"
)

############################
# Load SuSiE results
############################

susie_results <- readRDS(susie_results_file)

############################
# Extract credible sets
############################

summary_list <- list()

for (i in seq_along(susie_results)) {
  
  res <- susie_results[[i]]
  
  if (is.null(res)) next
  
  fit <- res$fit
  gwas_data <- res$gwas
  
  if (!is.null(fit$sets$cs)) {
    
    for (cs_id in names(fit$sets$cs)) {
      
      snp_idx <- fit$sets$cs[[cs_id]]
      
      cs_snps <- gwas_data[snp_idx, ]
      
      cs_snps$Trait <- "Nteeth"
      cs_snps$PIP <- fit$pip[snp_idx]
      cs_snps$CS <- cs_id
      cs_snps$Locus <- res$locus
      
      if (!is.null(fit$sets$purity)) {
        cs_snps$Purity <- fit$sets$purity[cs_id, "min.abs.corr"]
      } else {
        cs_snps$Purity <- NA
      }
      
      summary_list[[length(summary_list) + 1]] <- cs_snps
    }
  }
}

############################
# Save final summary
############################

if (length(summary_list) > 0) {
  
  susie_final_summary <- bind_rows(summary_list)
  
  susie_final_summary <- susie_final_summary %>%
    arrange(desc(PIP)) %>%
    dplyr::select(
      Trait,
      Locus,
      MarkerName,
      CS,
      PIP,
      Purity,
      everything()
    )
  
  write.csv(
    susie_final_summary,
    file = summary_output_file,
    row.names = FALSE
  )
  
  cat("Final SuSiE summary saved to:", summary_output_file, "\n")
  
} else {
  
  cat("No credible sets found.\n")
}