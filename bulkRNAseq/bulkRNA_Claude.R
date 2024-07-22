# Load required libraries
library(edgeR)
library(biomaRt)
library(cqn)
library(sva)

# Set environment variables
COUNT_TABLES_DIR <- "count_tables_human/"
METADATA_FILE <- "bulkRNA_metadata.csv"
OUTPUT_DIR <- "bulkRNA/"

# Set working directory
setwd(OUTPUT_DIR)

# Read metadata
metadata <- read.csv(METADATA_FILE, header = TRUE, stringsAsFactors = FALSE)

# 1. Read count tables
read_count_tables <- function(dir_path) {
  count_files <- list.files(dir_path, pattern = "txt$", full.names = TRUE)
  
  count_list <- lapply(count_files, function(file) {
    sample_name <- tools::file_path_sans_ext(basename(file))
    table <- read.table(file, header = TRUE, row.names = 1)
    colnames(table) <- sample_name
    return(table)
  })
  
  count_table_unified <- do.call(cbind, count_list)
  return(count_table_unified[-c(1:4), ])
}

count_table_unified <- read_count_tables(COUNT_TABLES_DIR)

# Visualize data before normalization
pdf("boxplot_before_cqn.pdf")
boxplot(log2(count_table_unified), range = 0, main = "Boxplot before CQN", las = 2)
dev.off()

# 2. Normalization with CQN
# Filter data
dataset <- 3
count_table_unified <- count_table_unified[, metadata$batch != dataset]
metadata <- metadata[metadata$batch != dataset, ]

count_table_unified <- count_table_unified[, metadata$length == 76]
metadata <- metadata[metadata$length == 76, ]

# Filter genes
filter_genes <- function(count_table) {
  # Remove unmapped genes
  count_table <- count_table[rowSums(count_table) > 0, ]
  
  # Keep genes with CPM >= 1 in at least 3 samples
  cpm_table <- cpm(count_table)
  count_table <- count_table[rowSums(cpm_table >= 1) >= 3, ]
  
  return(count_table)
}

count_table_filtered <- filter_genes(count_table_unified)

# Get Ensembl gene information
get_ensembl_info <- function(gene_ids) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ensembl_info <- getBM(
    attributes = c('ensembl_gene_id', 'percentage_gene_gc_content', 'start_position', 'end_position', 'hgnc_symbol', 'description'),
    filters = 'ensembl_gene_id',
    values = gene_ids,
    mart = mart
  )
  return(ensembl_info[order(ensembl_info$ensembl_gene_id), ])
}

ensembl_info <- get_ensembl_info(rownames(count_table_filtered))

# Prepare data for CQN
count_table_ensembl <- count_table_filtered[rownames(count_table_filtered) %in% ensembl_info$ensembl_gene_id, ]
count_table_ensembl <- count_table_ensembl[order(rownames(count_table_ensembl)), ]

stopifnot(all(ensembl_info$ensembl_gene_id == rownames(count_table_ensembl)))

genes_length <- ensembl_info$end_position - ensembl_info$start_position

# Perform CQN normalization
cqn_result <- cqn(
  counts = count_table_ensembl,
  x = ensembl_info$percentage_gene_gc_content,
  lengths = genes_length,
  sizeFactors = colSums(count_table_ensembl),
  verbose = TRUE
)

RPKM_cqn <- cqn_result$y + cqn_result$offset

# Scale matrix
RPKM_zscore <- t(scale(t(RPKM_cqn)))

# Verify sample names match
stopifnot(all(colnames(RPKM_zscore) == metadata$IdRun))

# Save results
save(RPKM_zscore, file = "RPKM_zscore.RData")
write.csv(RPKM_zscore, file = "RPKM_zscore.csv")