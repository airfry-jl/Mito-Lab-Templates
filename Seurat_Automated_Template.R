# Library -----------------------------------------------------------------
suppressMessages(
  {
    library(dplyr)
    library(Seurat)
    library(patchwork)
    library(tibble)
    library(DoubletFinder)
  }
)

# Load Path: Make sure the paths are right for your objets-----------------
fastq_path <-"~/Desktop/UCSF/Steven/Hum Muscle Spatial/Human Muscle/"
fastq_folder <- list.dirs(fastq_path, full.names = FALSE, recursive = FALSE, pattern = NULL)
#If you want to run a specific subset of sc, then specify a pattern, if not leave as null
obj_path <- paste0(fastq_path,"Objects/")
proj_name <- "Hum_Muscle"


# Runs all folders and outputs seurat objects -----------------------------
for (i in 1:length(fastq_folder)){

fastq <- Read10X(data.dir = paste0(fastq_path,fastq_folder[i],"/filtered_feature_bc_matrix/"))
# Initialize the Seurat object with the raw (non-normalized data).
obj <- CreateSeuratObject(counts = fastq, project = proj_name, min.cells = 3, min.features = 200)

# QC ----------------------------------------------------------------------
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
feature_RNA <- 0.83*max(obj@meta.data$nFeature_RNA)
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < feature_RNA & percent.mt < 10)

# Normalize Data ----------------------------------------------------------
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify Variable Features ----------------------------------------------
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(obj), 10)  #identify the 10 most highly variable genes

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale Data --------------------------------------------------------------
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)

# Linear Dimensional Reduction --------------------------------------------
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 60)

#Figure out how many dims to include

ElbowPlot(obj,ndims = 50)

# For loop for resolutions 0.3 - 1.2  ------------------------------------------
for (j in 3:12){
  #cluster!
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = (j/10))
  head(Idents(obj), 5)
  
  #non-linear reduction to umap
  obj <- RunUMAP(obj, dims = 1:30)
  DimPlot(obj, reduction = "umap", label = TRUE)
}

Idents(obj)<-"RNA_snn_res.0.8"
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1) #preliminary clustering

# Save merged object ------------------------------------------------------
setwd(obj_path)
saveRDS(obj, file = paste0(fastq_folder[i],".rds"))

}
