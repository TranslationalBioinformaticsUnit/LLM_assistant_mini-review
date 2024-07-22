# Define data directory
data_dir <- "/brahms/mollag/practice/filtered_gene_bc_matrices/hg19/"

# Load raw expression matrix
raw_expression_matrix <- Read10X(data.dir = data_dir)

# Create Seurat object
pbmc <- CreateSeuratObject(counts = raw_expression_matrix, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)

# Calculate mitochondrial content
pbmc[["percent_mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# ... rest of your code ...
