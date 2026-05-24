############################
# Prepare GWAS summary data
############################

rm(list = ls())

############################
# Paths
############################

root <- getwd()

############################
# DMFS
############################

dmfs_summary <- read.table(
  file.path(root, "data_samples/EUR_DMFS_incl_HSHSSOL.txt"),
  header = TRUE
)

head(dmfs_summary)
dim(dmfs_summary)

saveRDS(
  dmfs_summary,
  file = file.path(root, "data_samples/dmfs_summary.RDS")
)

############################
# Nteeth and Perio
# Add these once the GWAS .txt files are available
############################

# nteeth_summary <- read.table(
#   file.path(root, "data_samples/EUR_Nteeth_incl_HSHSSOL.txt"),
#   header = TRUE
# )
#
# saveRDS(
#   nteeth_summary,
#   file = file.path(root, "data_samples/nteeth_summary.RDS")
# )

# perio_summary <- read.table(
#   file.path(root, "data_samples/EUR_Perio_incl_HSHSSOL.txt"),
#   header = TRUE
# )
#
# saveRDS(
#   perio_summary,
#   file = file.path(root, "data_samples/perio_summary.RDS")
# )

############################
# Nteeth
############################

nteeth_summary <- read.table(
  file.path(root, "data_samples/EUR_nteeth_incl_HCHSSOL.txt"),
  header = TRUE
)

head(nteeth_summary)
dim(nteeth_summary)

saveRDS(
  nteeth_summary,
  file = file.path(root, "data_samples/nteeth_summary.RDS")
)
