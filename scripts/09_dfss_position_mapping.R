############################################################
# DFSS positional mapping
# Query Ensembl GRCh37 region-by-region to avoid 504 timeouts
############################################################

rm(list = ls())

library(biomaRt)
library(dplyr)

############################################################
# Paths
############################################################

root <- "/Users/divyapandey/Documents/GitHub/dental-genetics-susie"

input_file <- file.path(
  root,
  "results/dfss_susie_credible_sets.csv"
)

output_file <- file.path(
  root,
  "results/dfss_positional_mapping.csv"
)

############################################################
# Load credible-set results
############################################################

susie_final_summary <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)

susie_final_summary$CHR <- as.integer(
  susie_final_summary$CHR
)

susie_final_summary$BP <- as.integer(
  susie_final_summary$BP
)

cat(
  "Loaded",
  nrow(susie_final_summary),
  "credible-set rows\n"
)

############################################################
# Connect to Ensembl GRCh37
############################################################

ensembl <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "https://grch37.ensembl.org"
)

############################################################
# Unique SNP positions
############################################################

unique_snps <- susie_final_summary %>%
  dplyr::select(
    MarkerName,
    CHR,
    BP
  ) %>%
  distinct()

cat(
  "Unique SNP positions:",
  nrow(unique_snps),
  "\n"
)

############################################################
# Query one SNP region at a time
############################################################

gene_results <- list()

for (i in seq_len(nrow(unique_snps))) {
  
  snp <- unique_snps[i, ]
  
  chr_value <- snp$CHR
  bp_value <- snp$BP
  
  start_value <- max(
    bp_value - 50000,
    1
  )
  
  end_value <- bp_value + 50000
  
  cat(
    "Querying",
    i,
    "of",
    nrow(unique_snps),
    "|",
    snp$MarkerName,
    "| chr",
    chr_value,
    ":",
    start_value,
    "-",
    end_value,
    "\n"
  )
  
  gene_query <- tryCatch(
    {
      
      getBM(
        attributes = c(
          "ensembl_gene_id",
          "external_gene_name",
          "chromosome_name",
          "start_position",
          "end_position"
        ),
        filters = c(
          "chromosome_name",
          "start",
          "end"
        ),
        values = list(
          chr_value,
          start_value,
          end_value
        ),
        mart = ensembl
      )
    },
    
    error = function(e) {
      
      cat(
        "Query failed for",
        snp$MarkerName,
        ":",
        conditionMessage(e),
        "\n"
      )
      
      return(NULL)
    }
  )
  
  if (
    !is.null(gene_query) &&
    nrow(gene_query) > 0
  ) {
    
    gene_query$MarkerName <- snp$MarkerName
    gene_query$CHR <- chr_value
    gene_query$BP <- bp_value
    
    gene_results[[length(gene_results) + 1]] <- gene_query
  }
  
  # small pause to reduce repeated rapid requests
  Sys.sleep(1)
}

############################################################
# Combine gene mappings
############################################################

if (length(gene_results) == 0) {
  
  stop(
    "No positional mapping results were returned."
  )
}

gene_mapping <- bind_rows(
  gene_results
)

############################################################
# Merge back with fine-mapping results
############################################################

positional_mapping <- susie_final_summary %>%
  left_join(
    gene_mapping,
    by = c(
      "MarkerName",
      "CHR",
      "BP"
    )
  ) %>%
  arrange(
    desc(PIP),
    MarkerName,
    external_gene_name
  )

############################################################
# Save
############################################################

write.csv(
  positional_mapping,
  output_file,
  row.names = FALSE
)

cat(
  "\nSaved DFSS positional mapping to:\n",
  output_file,
  "\n"
)

cat(
  "Rows:",
  nrow(positional_mapping),
  "\n"
)

cat(
  "Unique SNPs:",
  dplyr::n_distinct(positional_mapping$MarkerName),
  "\n"
)

cat(
  "Unique mapped genes:",
  dplyr::n_distinct(
    positional_mapping$ensembl_gene_id[
      !is.na(positional_mapping$ensembl_gene_id)
    ]
  ),
  "\n"
)

print(
  positional_mapping %>%
    dplyr::select(
      Locus,
      MarkerName,
      CS,
      PIP,
      Purity,
      Converged,
      external_gene_name,
      ensembl_gene_id
    ) %>%
    head(20)
)