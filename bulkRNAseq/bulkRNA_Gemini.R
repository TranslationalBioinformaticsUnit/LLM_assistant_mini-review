# Define data paths
count_tables_path <- "count_tables_human/"
metadata_path <- "bulkRNA_metadata.csv"

# Read count tables
count_tables <- list.files(count_tables_path, pattern="txt")

# Read and format metadata
metadata <- read.table(metadata_path, header = TRUE, sep=",")

# Function to read a count table
read_count_table <- function(filename) {
  table <- read.table(paste0(count_tables_path, filename), sep="\t", header=TRUE)
  rownames(table) <- table[, 1]
  table[, -1]
}

# Read all count tables
count_table_unified <- lapply(count_tables, read_count_table)

# Combine count tables
count_table_unified <- do.call(cbind, count_table_unified)
colnames(count_table_unified) <- sapply(count_tables, function(x) strsplit(x, "_")[[1]][1])

# ... rest of your code ...

  
  