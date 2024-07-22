# Load required libraries
library(dplyr)
library(Seurat)
library(patchwork)

# Set random seed for reproducibility
set.seed(42)

# Define constants
DATA_DIR <- "/brahms/mollag/practice/filtered_gene_bc_matrices/hg19/"
PROJECT_NAME <- "pbmc3k"
MIN_CELLS <- 3
MIN_FEATURES <- 200
MT_PATTERN <- "^MT-"
MAX_FEATURES <- 2500
MAX_MT_PERCENT <- 5
VARIABLE_FEATURES_N <- 2000

# Load the PBMC dataset
pbmc_data <- Read10X(data.dir = DATA_DIR)

# Initialize the Seurat object
pbmc <- CreateSeuratObject(
  counts = pbmc_data, 
  project = PROJECT_NAME, 
  min.cells = MIN_CELLS, 
  min.features = MIN_FEATURES
)

# Calculate mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = MT_PATTERN)

# Visualize QC metrics
qc_violin <- VlnPlot(
  pbmc, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  ncol = 3
)

qc_scatter1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
qc_scatter2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

qc_plots <- qc_violin / (qc_scatter1 + qc_scatter2)
ggsave("qc_plots.pdf", qc_plots, width = 12, height = 10)

# Filter cells
pbmc <- subset(
  pbmc, 
  subset = nFeature_RNA > MIN_FEATURES & 
    nFeature_RNA < MAX_FEATURES & 
    percent.mt < MAX_MT_PERCENT
)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(
  pbmc, 
  selection.method = "vst", 
  nfeatures = VARIABLE_FEATURES_N
)

# Visualize variable features
top10 <- head(VariableFeatures(pbmc), 10)
var_features_plot <- VariableFeaturePlot(pbmc)
labeled_var_features_plot <- LabelPoints(
  plot = var_features_plot, 
  points = top10, 
  repel = TRUE
)

var_features_combined <- var_features_plot + labeled_var_features_plot
ggsave("variable_features.pdf", var_features_combined, width = 12, height = 6)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine PCA results
pca_results <- print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_processed.rds")

# Generate a summary report
summary_report <- capture.output({
  cat("PBMC Dataset Summary\n")
  cat("--------------------\n")
  print(pbmc)
  cat("\nTop Variable Features:\n")
  print(top10)
  cat("\nPCA Results:\n")
  print(pca_results)
})

writeLines(summary_report, "pbmc_analysis_summary.txt")