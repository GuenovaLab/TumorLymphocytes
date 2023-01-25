# loading packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("reticulate")
BiocManager::install("multtest")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")
BiocManager::install("depth")

install.packages("Seurat")
install.packages("SingleCellExperiment")
install.packages("depth")
install.packages("rgl")
install.packages("scales")
install.packages("robustbase")

library(multtest)
library(reticulate)

library(dplyr)
library(ggplot2)
library(cowplot)
library(SingleCellExperiment)
library(scater)
library(depth)
library(rgl)
library(robustbase)



library(Seurat)


setwd("M:/DER/RECHERCHE/GUENOVA_LAB/Project_Yun-Tsan_ADCC_paper/raw data and R code/Figure_6/6A_ridge plot/")

# loading the data
data <- read.csv("all_healthy_raw_counts_TCR_without_Tc.csv", header = TRUE, row.names = NULL) 
names <- make.unique(as.character(data$Gene))
rownames(data) <- names
data <- data[, -1] # get rid of old names

data <- as.matrix(data)
    
meta <- read.csv("all_healthy_raw_counts_TCR_without_Tc_group.csv", header = TRUE, row.names = NULL)



## match the metadata and count data
rownames(meta) <- meta$labels

### Check that sample names match in both files 
all(colnames(data) %in% rownames(meta)) 
all(colnames(data) == rownames(meta)) 

my_data <- CreateSeuratObject(counts = data, meta.data = meta, project = "MF", min.cells = 3, min.features = 200)


my_data[["percent.rbp"]] <- PercentageFeatureSet(my_data, pattern = "^RP[SL]")

VlnPlot(my_data, features = c("percent.rbp"), group.by = "SampleStage", ncol = 1, cols = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2")) + NoLegend()
VlnPlot(my_data, features = c("percent.rbp"), group.by = "StageNow", ncol = 1) + NoLegend()



# normalizing the data
my_data <- NormalizeData(my_data, normalization.method = "LogNormalize", scale.factor = 10000)


my_data <- FindVariableFeatures(my_data, selection.method = "vst", nfeatures = 2000)


# Scaling the data
all.genes <- rownames(my_data)
my_data <- ScaleData(my_data, features = all.genes)

# run PCA
my_data <- RunPCA(my_data, features = VariableFeatures(object = my_data))


# tSNE
my_data <- RunTSNE(my_data, dims = 1:30) 

DimPlot(my_data, reduction = "tsne")

VlnPlot(my_data, features = c("IL32"), group.by = "SampleStage")








tumor.cells <- subset(my_data, idents = c("MainClone", "RelateMainClone"))
tumor.cells.skin <- subset(tumor.cells, subset = compartment == "skin")


tumor.cells.skin <- RunTSNE(tumor.cells.skin, dims = 1:30)

DimPlot(tumor.cells.skin, reduction = "tsne", group.by = "StageNow")

VlnPlot(tumor.cells.skin, features = c("IL32"), group.by = "StageNow")

# find cluster
clone <- FindNeighbors(tumor.cells.skin, dims = 1:30)
clone <- FindClusters(clone, resolution = 0.5)
DimPlot(clone)
table(Idents(clone))

cluster_info <- Idents(clone)


write.csv(cluster_info, file="only_skin_tumor_cluster_IDs.csv", col.names = TRUE, row.names = TRUE)

VlnPlot(clone, features = c("HLA-A", "HLA-B", "HLA-C"), ncol = 3)


tumor.markers <- FindAllMarkers(clone, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tumor.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)



tumor.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(clone, features = top10$gene) + NoLegend()





