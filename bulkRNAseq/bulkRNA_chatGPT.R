library("edgeR")
library("biomaRt")
library("cqn")
library("sva")

# Set environment
set_environment <- function(count_tables_dir, metadata_path) {
  setwd(dirname(metadata_path))
  metadata <- read.table(metadata_path, header = TRUE, sep = ",")
  count_tables <- dir(count_tables_dir, pattern = "txt")
  return(list(metadata = metadata, count_tables = count_tables, count_tables_dir = count_tables_dir))
}

# Read count tables
read_count_tables <- function(count_tables, count_tables_dir) {
  sample_name <- unlist(strsplit(count_tables[1], "_"))
  count_table_unified <- read.table(paste0(count_tables_dir, count_tables[1]))
  rownames_count_table_unified <- count_table_unified[, 1]
  count_table_unified <- data.frame(count_table_unified[, -1])
  names(count_table_unified) <- sample_name[1]
  
  for (i in 2:length(count_tables)) {
    sample_name <- unlist(strsplit(count_tables[i], "_"))
    table <- read.table(paste0(count_tables_dir, count_tables[i]))
    table <- table[, -1]
    count_table_unified <- cbind(count_table_unified, table)
    colnames(count_table_unified)[i] <- sample_name[1]
  }
  
  rownames(count_table_unified) <- rownames_count_table_unified
  count_table_unified <- count_table_unified[-c(1, 2, 3, 4), ]
  return(count_table_unified)
}

# Plot boxplot before CQN
plot_boxplot <- function(count_table_unified) {
  boxplot(log2(count_table_unified), range = 0, main = "Boxplot before CQN", las = 2)
}

# Normalize using CQN
normalize_cqn <- function(count_table_unified, metadata) {
  dataset <- 3
  count_table_unified <- count_table_unified[, -which(metadata$batch == dataset)]
  metadata <- metadata[-which(metadata$batch == dataset), ]
  count_table_unified <- count_table_unified[, which(metadata$length == 76)]
  metadata <- metadata[which(metadata$length == 76), ]
  
  count_table_unified <- count_table_unified[which(rowSums(as.matrix(count_table_unified)) != 0), ]
  cpm <- cpm(count_table_unified)
  count_table_unified_cpm <- count_table_unified[which((rowSums(cpm >= 1)) > 3), ]
  
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ensembleID <- getBM(attributes = c('ensembl_gene_id', 'percentage_gene_gc_content', 'start_position', 'end_position', 'hgnc_symbol', 'description'), 
                      filters = c('ensembl_gene_id'), values = rownames(count_table_unified_cpm), mart = mart)
  
  count_table_unified_ensembl <- count_table_unified_cpm[which(rownames(count_table_unified_cpm) %in% ensembleID$ensembl_gene_id), ]
  ensembleID <- ensembleID[order(ensembleID$ensembl_gene_id), ]
  count_table_unified_ensembl_order <- count_table_unified_ensembl[order(rownames(count_table_unified_ensembl)), ]
  genes_length <- ensembleID$end_position - ensembleID$start_position
  
  count_table_unified.cqn <- cqn(count_table_unified_ensembl_order, lengths = genes_length, x = ensembleID$percentage_gene_gc_content, 
                                 sizeFactors = colSums(as.matrix(count_table_unified_ensembl_order)), verbose = TRUE)
  RPKM.cqn <- count_table_unified.cqn$y + count_table_unified.cqn$offset
  
  RPKM.zscore <- t(scale(t(RPKM.cqn), center = TRUE, scale = TRUE))
  stopifnot(all(colnames(RPKM.zscore) == metadata$IdRun))
  
  return(RPKM.zscore)
}

# Main function
main <- function() {
  count_tables_dir <- "count_tables_human/"
  metadata_path <- "bulkRNA/bulkRNA_metadata.csv"
  
  env <- set_environment(count_tables_dir, metadata_path)
  count_table_unified <- read_count_tables(env$count_tables, env$count_tables_dir)
  
  plot_boxplot(count_table_unified)
  
  RPKM.zscore <- normalize_cqn(count_table_unified, env$metadata)
  return(RPKM.zscore)
}

# Run the main function
RPKM.zscore <- main()
