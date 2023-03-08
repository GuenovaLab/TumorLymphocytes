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


# normalizing the data
my_data <- NormalizeData(my_data, normalization.method = "LogNormalize", scale.factor = 10000)


my_data <- FindVariableFeatures(my_data, selection.method = "vst", nfeatures = 2000)


# Scaling the data
all.genes <- rownames(my_data)
my_data <- ScaleData(my_data, features = all.genes)

# run PCA
my_data <- RunPCA(my_data, features = VariableFeatures(object = my_data))


# tSNE
my_data <- RunTSNE(my_data) 

# ordering the identity
table(my_data@meta.data$definitions)
my_data$definitions <- factor(x = my_data$definitions, levels = c("Healthy", "SingleBystanders", "BystanderGroups", "RelateMainClone", "MainClone"))

RidgePlot(my_data, features = c("IL2RA",
                                "IL2RB",
                                "IL2RG",
                                "IL4I1",
                                "IL4R",
                                "IL6R",
                                "IL7R",
                                "IL9R",
                                "IL10RA",
                                "IL10RB",
                                "IL11RA",
                                "IL16",
                                "IL17RB",
                                "IL18",
                                "IL32"), group.by = "definitions", cols = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"), ncol = 5, y.max = 5)





RidgePlot(my_data, features = c("CD36"), group.by = "definitions", cols = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"), ncol = 1, y.max = 5)



VlnPlot(my_data, features = c("CD36"), group.by = "definitions", cols = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"), ncol = 1, y.max = 1.5)




