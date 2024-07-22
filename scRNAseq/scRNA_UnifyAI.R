# Load necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)

# Set data directory
data_dir <- "/brahms/mollag/practice/filtered_gene_bc_matrices/hg19/"

# Load the PBMC dataset
pbmc_data <- Read10X(data.dir = data_dir)

# Initialize the Seurat object with the raw (non-normalized) data
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Calculate the percentage of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Feature-Feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_plot <- plot1 + plot2
combined_plot

# Subsetting the data based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10_features <- head(VariableFeatures(pbmc), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10_features, repel = TRUE)
combined_feature_plot <- plot1 + plot2
combined_feature_plot

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# Examine and visualize PCA results
pca_results <- pbmc[["pca"]]
print(pca_results, dims = 1:5, nfeatures = 5)