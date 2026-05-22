# Dental Genetics SuSiE Fine-Mapping

This repository contains an R workflow for SuSiE-RSS fine-mapping of dental GWAS traits using summary statistics and LD matrices.

The current analysis is implemented for **DMFS**. The workflow structure can be extended to additional dental traits including **Nteeth** and **Perio**.

## Traits

- DMFS
- Nteeth
- Perio

## Workflow

The analysis is organized into numbered R scripts:

1. `01_prepare_input_data.R`  
   Converts GWAS summary statistic `.txt` files into `.RDS` files for faster loading in R.

2. `02_run_susie.R`  
   Runs SuSiE-RSS fine-mapping using GWAS summary statistics and LD matrices.

3. `03_summarize_results.R`  
   Extracts credible set SNPs, posterior inclusion probabilities (PIPs), and credible set purity.

4. `04_position_mapping.R`  
   Maps SuSiE credible set SNPs to nearest hg19 genes.

5. `05_eqtl_colocalization.R`  
   Queries GTEx eQTL evidence for mapped SNP-gene pairs.

6. `06_display_results.R`  
   Displays final result tables using interactive `DT::datatable()` tables.

## Input files

Raw GWAS summary statistics should be placed in:

```text
data_samples/