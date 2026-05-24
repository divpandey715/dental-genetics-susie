############################
# eQTL colocalization / GTEx lookup
# Trait: Nteeth
############################

rm(list = ls())

library(data.table)
library(dplyr)
library(gtexr)

############################
# Paths
############################

root <- getwd()

positional_mapping_file <- file.path(
  root,
  "results/nteeth_positional_mapping.csv"
)

ld_dir <- file.path(
  root,
  "LDMatrix/04_compute_LD_for_Shungin2019_EUR_Nteeth_output"
)

output_file <- file.path(
  root,
  "results/nteeth_eqtl_colocalization_results.csv"
)

############################
# Load positional mapping results
############################

positional_mapping <- read.csv(positional_mapping_file)

head(positional_mapping)
dim(positional_mapping)

############################
# Convert chr:position to rsID using LD id_map files
############################
# The id_map files contain rsID and coordinate mappings.
############################

id_map_files <- list.files(
  ld_dir,
  pattern = "\\.id_map$",
  full.names = TRUE
)

length(id_map_files)

id_map_list <- lapply(id_map_files, function(file) {
  
  id_map <- fread(
    file,
    header = FALSE
  )
  
  colnames(id_map) <- c("rsid", "MarkerName")
  
  return(id_map)
})

id_map_all <- rbindlist(
  id_map_list,
  fill = TRUE
)

id_map_all <- id_map_all %>%
  distinct(MarkerName, .keep_all = TRUE)

head(id_map_all)
dim(id_map_all)

############################
# Merge rsIDs into positional mapping
############################

eqtl_final_summary <- positional_mapping %>%
  left_join(
    id_map_all,
    by = "MarkerName"
  )

head(eqtl_final_summary)
dim(eqtl_final_summary)

############################
# Filter SNPs with rsID and assigned gene
############################

to_process <- eqtl_final_summary %>%
  filter(
    !is.na(rsid),
    !is.na(ensembl_gene_id),
    rsid != "",
    ensembl_gene_id != ""
  ) %>%
  distinct(
    rsid,
    ensembl_gene_id,
    .keep_all = TRUE
  )

head(to_process)
dim(to_process)

############################
# Fetch eQTL data from GTEx
############################
# You can change the tissue if needed.
# Examples: "Whole_Blood", "Minor_Salivary_Gland"
############################

target_tissue <- "Whole_Blood"

fetch_eqtl <- function(rsid, gene_id, tissue) {
  
  tryCatch(
    {
      res <- calculate_expression_quantitative_trait_loci(
        tissueSiteDetailId = tissue,
        geneId = gene_id,
        variantId = rsid
      )
      
      if (nrow(res) > 0) {
        return(res)
      } else {
        return(NULL)
      }
    },
    error = function(e) {
      return(NULL)
    }
  )
}

eqtl_results_list <- mapply(
  fetch_eqtl,
  to_process$rsid,
  to_process$ensembl_gene_id,
  MoreArgs = list(
    tissue = target_tissue
  ),
  SIMPLIFY = FALSE
)

names(eqtl_results_list) <- to_process$rsid

############################
# Combine eQTL results
############################

valid_results <- eqtl_results_list[
  !sapply(eqtl_results_list, is.null)
]

if (length(valid_results) > 0) {
  
  eqtl_results <- bind_rows(
    valid_results,
    .id = "rsid"
  )
  
  eqtl_results <- eqtl_results %>%
    left_join(
      to_process,
      by = "rsid"
    )
  
  write.csv(
    eqtl_results,
    file = output_file,
    row.names = FALSE
  )
  
  cat("eQTL results saved to:", output_file, "\n")
  
  head(eqtl_results)
  dim(eqtl_results)
  
} else {
  
  cat("No eQTL matches found in the specified tissue:", target_tissue, "\n")
}