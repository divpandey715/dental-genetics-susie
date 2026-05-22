############################
# Display final result tables
# Trait: DMFS
############################

rm(list = ls())

library(dplyr)
library(DT)

############################
# Paths
############################

root <- getwd()

susie_summary_file <- file.path(
  root,
  "results/dmfs_susie_final_summary.csv"
)

positional_mapping_file <- file.path(
  root,
  "results/dmfs_positional_mapping.csv"
)

############################
# 1. SuSiE credible set table
############################

susie_final_summary <- read.csv(susie_summary_file)

susie_final_summary %>%
  arrange(desc(PIP)) %>%
  dplyr::select(
    Locus,
    MarkerName,
    CS,
    PIP,
    Purity,
    everything()
  ) %>%
  datatable(
    extensions = "Buttons",
    caption = "",
    options = list(
      dom = "Blfrtip",
      buttons = c("copy", "csv", "excel", "pdf", "print"),
      lengthMenu = list(
        c(10, 25, 50, -1),
        c(10, 25, 50, "All")
      ),
      scrollX = TRUE
    )
  )

############################
# 2. Position mapping table
############################

positional_mapping <- read.csv(positional_mapping_file)

positional_mapping %>%
  dplyr::select(
    Locus,
    MarkerName,
    PIP,
    external_gene_name,
    ensembl_gene_id,
    everything()
  ) %>%
  arrange(desc(PIP)) %>%
  datatable(
    caption = "SuSiE Results with Nearest hg19 Genes",
    options = list(
      pageLength = 10,
      scrollX = TRUE
    )
  )

############################
# 3. eQTL colocalization display
############################

target_tissue <- "Whole_Blood"

eqtl_report_file <- file.path(
  root,
  "results/dmfs_eqtl_colocalization_results.csv"
)

if (file.exists(eqtl_report_file)) {
  
  final_coloc_report <- read.csv(eqtl_report_file)
  
  datatable(
    final_coloc_report,
    caption = paste(
      "eQTL Colocalization Summary (Tissue:",
      target_tissue,
      ")"
    ),
    options = list(
      pageLength = 10,
      scrollX = TRUE
    )
  )
  
} else {
  
  print("No eQTL matches found in the specified tissue.")
}