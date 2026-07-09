############################################################
# DFSS eQTL overlap / direction-consistency analysis
# Uses local DFSS .id_map files for rsID mapping
# Preserves query metadata directly with GTEx results
############################################################

rm(list = ls())

library(gtexr)
library(dplyr)
library(data.table)

############################################################
# Paths
############################################################

root <- "/Users/divyapandey/Documents/GitHub/dental-genetics-susie"

position_file <- file.path(
  root,
  "results/dfss_positional_mapping.csv"
)

ld_dir <- file.path(
  root,
  "LDMatrix/04_compute_LD_for_Shungin2019_EUR_DFSS_output"
)

output_file <- file.path(
  root,
  "results/dfss_eqtl_colocalization_summary.csv"
)

target_tissue <- "Whole_Blood"

############################################################
# Load positional mapping
############################################################

if (!file.exists(position_file)) {
  
  stop(
    "Position mapping file not found: ",
    position_file
  )
}

positional_mapping <- read.csv(
  position_file,
  stringsAsFactors = FALSE
)

cat(
  "Loaded",
  nrow(positional_mapping),
  "positional mapping rows\n"
)

############################################################
# Check required columns
############################################################

required_columns <- c(
  "MarkerName",
  "PIP",
  "Effect",
  "ensembl_gene_id",
  "external_gene_name"
)

missing_columns <- setdiff(
  required_columns,
  colnames(positional_mapping)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

############################################################
# Keep valid gene mappings
############################################################

positional_mapping <- positional_mapping %>%
  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  )

cat(
  "Rows with valid gene mappings:",
  nrow(positional_mapping),
  "\n"
)

if (nrow(positional_mapping) == 0) {
  
  stop(
    "No valid SNP-gene mappings are available."
  )
}

############################################################
# Read all DFSS id_map files
############################################################

map_files <- list.files(
  ld_dir,
  pattern = "\\.id_map$",
  full.names = TRUE
)

cat(
  "Number of DFSS id_map files:",
  length(map_files),
  "\n"
)

if (length(map_files) == 0) {
  
  stop(
    "No .id_map files found."
  )
}

############################################################
# Read mappings
############################################################

map_list <- vector(
  mode = "list",
  length = length(map_files)
)

for (i in seq_along(map_files)) {
  
  cat(
    "Reading id_map",
    i,
    "of",
    length(map_files),
    "\n"
  )
  
  map_data <- fread(
    map_files[i],
    header = FALSE
  )
  
  if (ncol(map_data) < 2) {
    next
  }
  
  map_data <- map_data[, 1:2]
  
  setnames(
    map_data,
    c(
      "rsid",
      "MarkerName"
    )
  )
  
  map_list[[i]] <- map_data
}

############################################################
# Combine mapping files
############################################################

valid_maps <- map_list[
  !vapply(
    map_list,
    is.null,
    logical(1)
  )
]

if (length(valid_maps) == 0) {
  
  stop(
    "No valid id_map data could be read."
  )
}

snp_lookup <- bind_rows(
  valid_maps
) %>%
  distinct(
    rsid,
    MarkerName
  )

cat(
  "Unique rsID-coordinate mappings:",
  nrow(snp_lookup),
  "\n"
)

############################################################
# Add rsIDs to positional mapping
############################################################

eqtl_final_summary <- positional_mapping %>%
  left_join(
    snp_lookup,
    by = "MarkerName"
  )

cat(
  "Rows with rsIDs:",
  sum(!is.na(eqtl_final_summary$rsid)),
  "of",
  nrow(eqtl_final_summary),
  "\n"
)

############################################################
# Prepare unique SNP-gene pairs
############################################################

to_process <- eqtl_final_summary %>%
  filter(
    !is.na(rsid),
    rsid != "",
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%
  distinct(
    rsid,
    ensembl_gene_id,
    .keep_all = TRUE
  )

cat(
  "Unique SNP-gene pairs to query:",
  nrow(to_process),
  "\n"
)

if (nrow(to_process) == 0) {
  
  stop(
    "No valid SNP-gene pairs available."
  )
}

############################################################
# GTEx query function
############################################################

fetch_eqtl <- function(
    rsid,
    gene_id,
    tissue
) {
  
  result <- tryCatch(
    {
      
      calculate_expression_quantitative_trait_loci(
        tissueSiteDetailId = tissue,
        gencodeId = gene_id,
        variantId = rsid
      )
      
    },
    error = function(e) {
      
      cat(
        "GTEx query failed:",
        rsid,
        "|",
        gene_id,
        "|",
        conditionMessage(e),
        "\n"
      )
      
      return(NULL)
    }
  )
  
  if (
    is.null(result) ||
    nrow(result) == 0
  ) {
    
    return(NULL)
  }
  
  return(
    as.data.frame(result)
  )
}

############################################################
# Query GTEx and attach query metadata immediately
############################################################

cat(
  "\nQuerying GTEx tissue:",
  target_tissue,
  "\n"
)

eqtl_results_list <- vector(
  mode = "list",
  length = nrow(to_process)
)

for (i in seq_len(nrow(to_process))) {
  
  cat(
    "GTEx query",
    i,
    "of",
    nrow(to_process),
    "|",
    to_process$rsid[i],
    "|",
    to_process$external_gene_name[i],
    "\n"
  )
  
  gtex_result <- fetch_eqtl(
    rsid = to_process$rsid[i],
    gene_id = to_process$ensembl_gene_id[i],
    tissue = target_tissue
  )
  
  ##########################################################
  # Attach original fine-mapping metadata
  ##########################################################
  
  if (
    !is.null(gtex_result) &&
    nrow(gtex_result) > 0
  ) {
    
    metadata_row <- to_process[
      rep(i, nrow(gtex_result)),
      ,
      drop = FALSE
    ]
    
    eqtl_results_list[[i]] <- bind_cols(
      metadata_row,
      gtex_result
    )
  }
  
  Sys.sleep(0.5)
}

############################################################
# Remove NULL results
############################################################

non_null_results <- eqtl_results_list[
  !vapply(
    eqtl_results_list,
    is.null,
    logical(1)
  )
]

############################################################
# Combine GTEx results
############################################################

if (length(non_null_results) > 0) {
  
  final_coloc_report <- bind_rows(
    non_null_results
  )
  
} else {
  
  final_coloc_report <- data.frame()
}

cat(
  "\nGTEx result rows returned:",
  nrow(final_coloc_report),
  "\n"
)

############################################################
# Analyze returned GTEx results
############################################################

if (nrow(final_coloc_report) > 0) {
  
  ##########################################################
  # Check GTEx output columns
  ##########################################################
  
  cat(
    "\nGTEx columns returned:\n"
  )
  
  print(
    colnames(final_coloc_report)
  )
  
  ##########################################################
  # Verify necessary GTEx result fields
  ##########################################################
  
  if (!"nes" %in% colnames(final_coloc_report)) {
    
    stop(
      "GTEx results were returned, but no 'nes' column was found."
    )
  }
  
  if (!"pValue" %in% colnames(final_coloc_report)) {
    
    stop(
      "GTEx results were returned, but no 'pValue' column was found."
    )
  }
  
  ##########################################################
  # Direction consistency and nominal significance
  ##########################################################
  
  final_coloc_report <- final_coloc_report %>%
    mutate(
      Direction_Match =
        sign(Effect) == sign(nes),
      
      eQTL_Sig =
        ifelse(
          pValue < 0.05,
          "Yes",
          "No"
        )
    )
  
  ##########################################################
  # Rename result columns
  ##########################################################
  
  final_coloc_report <- final_coloc_report %>%
    rename(
      GWAS_Effect = Effect,
      eQTL_NES = nes,
      eQTL_P = pValue
    )
  
  ##########################################################
  # Arrange useful columns first
  ##########################################################
  
  priority_columns <- c(
    "Locus",
    "MarkerName",
    "rsid",
    "CS",
    "PIP",
    "Purity",
    "Converged",
    "external_gene_name",
    "ensembl_gene_id",
    "GWAS_Effect",
    "eQTL_NES",
    "eQTL_P",
    "Direction_Match",
    "eQTL_Sig"
  )
  
  priority_columns <- priority_columns[
    priority_columns %in%
      colnames(final_coloc_report)
  ]
  
  final_coloc_report <- final_coloc_report %>%
    dplyr::select(
      all_of(priority_columns),
      everything()
    ) %>%
    arrange(
      eQTL_P
    )
  
  ##########################################################
  # Save results
  ##########################################################
  ############################################################
  # Remove list-columns before writing CSV
  ############################################################
  
  list_columns <- vapply(
    final_coloc_report,
    is.list,
    logical(1)
  )
  
  if (any(list_columns)) {
    
    cat(
      "Removing list-columns before CSV export:",
      paste(
        names(final_coloc_report)[list_columns],
        collapse = ", "
      ),
      "\n"
    )
    
    final_coloc_report <- final_coloc_report[
      ,
      !list_columns,
      drop = FALSE
    ]
  }
  
  
  write.csv(
    final_coloc_report,
    output_file,
    row.names = FALSE
  )
  
  cat(
    "\nSaved DFSS eQTL summary to:\n",
    output_file,
    "\n"
  )
  
  cat(
    "Number of eQTL matches:",
    nrow(final_coloc_report),
    "\n"
  )
  
  cat(
    "\nDFSS eQTL results:\n"
  )
  
  print(
    final_coloc_report %>%
      dplyr::select(
        MarkerName,
        rsid,
        PIP,
        external_gene_name,
        GWAS_Effect,
        eQTL_NES,
        eQTL_P,
        Direction_Match,
        eQTL_Sig
      )
  )
  
  ############################################################
  # Zero-result case
  ############################################################
  
} else {
  
  cat(
    "\nNo eQTL matches found in tissue:",
    target_tissue,
    "\n"
  )
  
  final_coloc_report <- data.frame(
    Locus = character(0),
    MarkerName = character(0),
    rsid = character(0),
    CS = character(0),
    PIP = numeric(0),
    Purity = numeric(0),
    Converged = logical(0),
    external_gene_name = character(0),
    ensembl_gene_id = character(0),
    GWAS_Effect = numeric(0),
    eQTL_NES = numeric(0),
    eQTL_P = numeric(0),
    Direction_Match = logical(0),
    eQTL_Sig = character(0)
  )
  
  write.csv(
    final_coloc_report,
    output_file,
    row.names = FALSE
  )
  
  cat(
    "Created empty eQTL summary CSV for workflowr.\n"
  )
}

############################################################
# Final report
############################################################

cat(
  "\n====================================\n"
)

cat(
  "DFSS eQTL analysis finished\n"
)

cat(
  "Target tissue:",
  target_tissue,
  "\n"
)

cat(
  "Output file:",
  output_file,
  "\n"
)

cat(
  "====================================\n"
)