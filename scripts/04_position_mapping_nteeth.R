############################
# Positional mapping
# Trait: Nteeth
# Uses hg19 gene coordinates from Bioconductor
############################

rm(list = ls())

library(data.table)
library(dplyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)

############################
# Paths
############################

root <- getwd()

summary_file <- file.path(
  root,
  "results/nteeth_susie_final_summary.csv"
)

output_file <- file.path(
  root,
  "results/nteeth_positional_mapping.csv"
)

############################
# Load SuSiE final summary
############################

susie_final_summary <- read.csv(summary_file)

head(susie_final_summary)
dim(susie_final_summary)

############################
# Create SNP genomic ranges
############################

snp_ranges <- GRanges(
  seqnames = paste0("chr", susie_final_summary$CHR),
  ranges = IRanges(
    start = susie_final_summary$BP,
    end = susie_final_summary$BP
  )
)

############################
# Load hg19 gene annotations
############################

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

gene_ranges <- genes(txdb)

############################
# Find nearest gene
############################

nearest_gene_index <- nearest(
  snp_ranges,
  gene_ranges
)

nearest_genes <- gene_ranges[nearest_gene_index]

entrez_ids <- nearest_genes$gene_id

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = entrez_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

ensembl_ids <- mapIds(
  org.Hs.eg.db,
  keys = entrez_ids,
  column = "ENSEMBL",
  keytype = "ENTREZID",
  multiVals = "first"
)

############################
# Add gene columns
############################

susie_final_summary$external_gene_name <- as.character(gene_symbols)
susie_final_summary$ensembl_gene_id <- as.character(ensembl_ids)
susie_final_summary$entrez_gene_id <- as.character(entrez_ids)

susie_final_summary$gene_chr <- as.character(seqnames(nearest_genes))
susie_final_summary$gene_start <- start(nearest_genes)
susie_final_summary$gene_end <- end(nearest_genes)

susie_final_summary$distance_to_gene <- distance(
  snp_ranges,
  nearest_genes
)

############################
# Format output
############################

positional_mapping <- susie_final_summary %>%
  dplyr::select(
    Trait,
    Locus,
    MarkerName,
    PIP,
    external_gene_name,
    ensembl_gene_id,
    entrez_gene_id,
    distance_to_gene,
    everything()
  ) %>%
  arrange(desc(PIP), P.value)

############################
# Save results
############################

write.csv(
  positional_mapping,
  file = output_file,
  row.names = FALSE
)

cat("Positional mapping results saved to:", output_file, "\n")

head(positional_mapping)
dim(positional_mapping)