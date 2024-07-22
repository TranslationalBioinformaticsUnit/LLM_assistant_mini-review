library("edgeR")
library("biomaRt")
library("cqn")
library("sva")

# Set environment
data_directory <- "bulkRNA/"
count_tables_dir <- paste0(data_directory, "count_tables_human/")
setwd(data_directory)

metadata <- read.csv("bulkRNA_metadata.csv")

# Function to read and process count tables
read_and_process_counts <- function(filename, directory) {
  filepath <- file.path(directory, filename)
  count_data <- read.table(filepath, header = FALSE)
  gene_ids <- count_data[, 1]
  counts <- data.frame(count_data[, -1])
  colnames(counts) <- unlist(strsplit(filename, "_"))[1]
  list(gene_ids = gene_ids, counts = counts)
}

# Read all count tables
count_files <- list.files(count_tables_dir, pattern = "txt")
count_list <- lapply(count_files, read_and_process_counts, directory = count_tables_dir)

# Combine count tables
count_table_unified <- do.call(function(...) {
  args <- list(...)
  combined <- Reduce(function(x, y) {
    merge(x, y, by = "gene_ids", all = TRUE)
  }, args)
  rownames(combined) <- combined$gene_ids
  combined[, -1]
}, count_list)

# Preprocessing before normalization
count_table_unified <- count_table_unified[-c(1, 2, 3, 4), ]
boxplot(log2(count_table_unified + 1), range = 0, main = "Boxplot before CQN", las = 2)

# Normalization with CQN
valid_samples <- metadata$batch != 3 & metadata$length == 76
filtered_metadata <- metadata[valid_samples, ]
filtered_counts <- count_table_unified[, valid_samples]

# Filtering steps
is_expressed <- rowSums(filtered_counts) > 0
cpm_values <- cpm(filtered_counts)
is_above_threshold <- rowSums(cpm_values >= 1) > 3

final_filtered_counts <- filtered_counts[is_expressed & is_above_threshold, ]

# Fetching gene information from Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c('ensembl_gene_id', 'percentage_gene_gc_content', 'start_position', 'end_position', 'hgnc_symbol', 'description')
ensembl_data <- getBM(attributes = attributes, filters = 'ensembl_gene_id', values = rownames(final_filtered_counts), mart = mart)

# Continue with normalization and further analysis...
