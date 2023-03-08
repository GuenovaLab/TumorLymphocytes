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
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(SingleCellExperiment)
library(scater)
library(depth)
library(rgl)
library(robustbase)
library(stringr)


setwd("M:/DER/RECHERCHE/GUENOVA_LAB/Yun-Tsan Chang/CTCL_research/project_scRNA_seq/GeneVia/06_Raw_counts_corrected/")

# loading the data

data <- read.csv("Healthy_T_cells.csv", header = TRUE, row.names = NULL) 
names <- make.unique(as.character(data$Gene.HGNC))
rownames(data) <- names
data <- data[, -1] # get rid of old names


data_healthy <- as.matrix(data)





data_01_BE <- CreateSeuratObject(counts = data_01_BE, project = "01_BE", min.cells = 3, min.features = 200)
data_01_BL <- CreateSeuratObject(counts = data_01_BL, project = "01_BL", min.cells = 3, min.features = 200)
data_01_BLC <- CreateSeuratObject(counts = data_01_BLC, project = "01_BLC", min.cells = 3, min.features = 200)
data_01_SE <- CreateSeuratObject(counts = data_01_SE, project = "01_SE", min.cells = 3, min.features = 200)
data_01_SL <- CreateSeuratObject(counts = data_01_SL, project = "01_SL", min.cells = 3, min.features = 200)
data_02_BE <- CreateSeuratObject(counts = data_02_BE, project = "02_BE", min.cells = 3, min.features = 200)
data_02_BL <- CreateSeuratObject(counts = data_02_BL, project = "02_BL", min.cells = 3, min.features = 200)
data_02_SE <- CreateSeuratObject(counts = data_02_SE, project = "02_SE", min.cells = 3, min.features = 200)
data_02_SL <- CreateSeuratObject(counts = data_02_SL, project = "02_SL", min.cells = 3, min.features = 200)
data_03_BE <- CreateSeuratObject(counts = data_03_BE, project = "03_BE", min.cells = 3, min.features = 200)
data_03_BL <- CreateSeuratObject(counts = data_03_BL, project = "03_BL", min.cells = 3, min.features = 200)
data_03_SE <- CreateSeuratObject(counts = data_03_SE, project = "03_SE", min.cells = 3, min.features = 200)
data_03_SL <- CreateSeuratObject(counts = data_03_SL, project = "03_SL", min.cells = 3, min.features = 200)
  
data_04_S <- CreateSeuratObject(counts = data_04_S, project = "04_S", min.cells = 3, min.features = 200)
data_06_S <- CreateSeuratObject(counts = data_06_S, project = "06_S", min.cells = 3, min.features = 200)
data_07_S <- CreateSeuratObject(counts = data_07_S, project = "07_S", min.cells = 3, min.features = 200)
data_08_S <- CreateSeuratObject(counts = data_08_S, project = "08_S", min.cells = 3, min.features = 200)
data_09_S <- CreateSeuratObject(counts = data_09_S, project = "09_S", min.cells = 3, min.features = 200)
data_10_S <- CreateSeuratObject(counts = data_10_S, project = "10_S", min.cells = 3, min.features = 200)
data_11_S <- CreateSeuratObject(counts = data_11_S, project = "11_S", min.cells = 3, min.features = 200)
data_12_S <- CreateSeuratObject(counts = data_12_S, project = "12_S", min.cells = 3, min.features = 200)
data_13_S <- CreateSeuratObject(counts = data_13_S, project = "13_S", min.cells = 3, min.features = 200)
data_14_S <- CreateSeuratObject(counts = data_14_S, project = "14_S", min.cells = 3, min.features = 200)
data_healthy <- CreateSeuratObject(counts = data_healthy, project = "healthy", min.cells = 3, min.features = 200)





data_01_BE[["percent.mt"]]  <- PercentageFeatureSet(data_01_BE, pattern = "^MT-")
data_01_BE[["percent.rbp"]] <- PercentageFeatureSet(data_01_BE, pattern = "^RP[SL]")
data_01_BE[["group"]] <- "01_BE"
data_01_BE[["compartment"]] <- "blood"

data_01_BL[["percent.mt"]]  <- PercentageFeatureSet(data_01_BL, pattern = "^MT-")
data_01_BL[["percent.rbp"]] <- PercentageFeatureSet(data_01_BL, pattern = "^RP[SL]")
data_01_BL[["group"]] <- "01_BL"
data_01_BL[["compartment"]] <- "blood"

data_01_BLC[["percent.mt"]]  <- PercentageFeatureSet(data_01_BLC, pattern = "^MT-")
data_01_BLC[["percent.rbp"]] <- PercentageFeatureSet(data_01_BLC, pattern = "^RP[SL]")
data_01_BLC[["group"]] <- "01_BLC"
data_01_BLC[["compartment"]] <- "blood"

data_01_SE[["percent.mt"]]  <- PercentageFeatureSet(data_01_SE, pattern = "^MT-")
data_01_SE[["percent.rbp"]] <- PercentageFeatureSet(data_01_SE, pattern = "^RP[SL]")
data_01_SE[["group"]] <- "01_SE"
data_01_SE[["compartment"]] <- "skin"

data_01_SL[["percent.mt"]]  <- PercentageFeatureSet(data_01_SL, pattern = "^MT-")
data_01_SL[["percent.rbp"]] <- PercentageFeatureSet(data_01_SL, pattern = "^RP[SL]")
data_01_SL[["group"]] <- "01_SL"
data_01_SL[["compartment"]] <- "skin"

data_02_BE[["percent.mt"]]  <- PercentageFeatureSet(data_02_BE, pattern = "^MT-")
data_02_BE[["percent.rbp"]] <- PercentageFeatureSet(data_02_BE, pattern = "^RP[SL]")
data_02_BE[["group"]] <- "02_BE"
data_02_BE[["compartment"]] <- "blood"

data_02_BL[["percent.mt"]]  <- PercentageFeatureSet(data_02_BL, pattern = "^MT-")
data_02_BL[["percent.rbp"]] <- PercentageFeatureSet(data_02_BL, pattern = "^RP[SL]")
data_02_BL[["group"]] <- "02_BL"
data_02_BL[["compartment"]] <- "blood"

data_02_SE[["percent.mt"]]  <- PercentageFeatureSet(data_02_SE, pattern = "^MT-")
data_02_SE[["percent.rbp"]] <- PercentageFeatureSet(data_02_SE, pattern = "^RP[SL]")
data_02_SE[["group"]] <- "02_SE"
data_02_SE[["compartment"]] <- "skin"

data_02_SL[["percent.mt"]]  <- PercentageFeatureSet(data_02_SL, pattern = "^MT-")
data_02_SL[["percent.rbp"]] <- PercentageFeatureSet(data_02_SL, pattern = "^RP[SL]")
data_02_SL[["group"]] <- "02_SL"
data_02_SL[["compartment"]] <- "skin"

data_03_BE[["percent.mt"]]  <- PercentageFeatureSet(data_03_BE, pattern = "^MT-")
data_03_BE[["percent.rbp"]] <- PercentageFeatureSet(data_03_BE, pattern = "^RP[SL]")
data_03_BE[["group"]] <- "03_BE"
data_03_BE[["compartment"]] <- "blood"

data_03_BL[["percent.mt"]]  <- PercentageFeatureSet(data_03_BL, pattern = "^MT-")
data_03_BL[["percent.rbp"]] <- PercentageFeatureSet(data_03_BL, pattern = "^RP[SL]")
data_03_BL[["group"]] <- "03_BL"
data_03_BL[["compartment"]] <- "blood"

data_03_SE[["percent.mt"]]  <- PercentageFeatureSet(data_03_SE, pattern = "^MT-")
data_03_SE[["percent.rbp"]] <- PercentageFeatureSet(data_03_SE, pattern = "^RP[SL]")
data_03_SE[["group"]] <- "03_SE"
data_03_SE[["compartment"]] <- "skin"

data_03_SL[["percent.mt"]]  <- PercentageFeatureSet(data_03_SL, pattern = "^MT-")
data_03_SL[["percent.rbp"]] <- PercentageFeatureSet(data_03_SL, pattern = "^RP[SL]")
data_03_SL[["group"]] <- "03_SL"
data_03_SL[["compartment"]] <- "skin"

data_04_S[["percent.mt"]]  <- PercentageFeatureSet(data_04_S, pattern = "^MT-")
data_04_S[["percent.rbp"]] <- PercentageFeatureSet(data_04_S, pattern = "^RP[SL]")
data_04_S[["group"]] <- "04_S"
data_04_S[["compartment"]] <- "skin"

data_06_S[["percent.mt"]]  <- PercentageFeatureSet(data_06_S, pattern = "^MT-")
data_06_S[["percent.rbp"]] <- PercentageFeatureSet(data_06_S, pattern = "^RP[SL]")
data_06_S[["group"]] <- "06_S"
data_06_S[["compartment"]] <- "skin"

data_07_S[["percent.mt"]]  <- PercentageFeatureSet(data_07_S, pattern = "^MT-")
data_07_S[["percent.rbp"]] <- PercentageFeatureSet(data_07_S, pattern = "^RP[SL]")
data_07_S[["group"]] <- "07_S"
data_07_S[["compartment"]] <- "skin"

data_08_S[["percent.mt"]]  <- PercentageFeatureSet(data_08_S, pattern = "^MT-")
data_08_S[["percent.rbp"]] <- PercentageFeatureSet(data_08_S, pattern = "^RP[SL]")
data_08_S[["group"]] <- "08_S"
data_08_S[["compartment"]] <- "skin"

data_09_S[["percent.mt"]]  <- PercentageFeatureSet(data_09_S, pattern = "^MT-")
data_09_S[["percent.rbp"]] <- PercentageFeatureSet(data_09_S, pattern = "^RP[SL]")
data_09_S[["group"]] <- "09_S"
data_09_S[["compartment"]] <- "skin"

data_10_S[["percent.mt"]]  <- PercentageFeatureSet(data_10_S, pattern = "^MT-")
data_10_S[["percent.rbp"]] <- PercentageFeatureSet(data_10_S, pattern = "^RP[SL]")
data_10_S[["group"]] <- "10_S"
data_10_S[["compartment"]] <- "skin"

data_11_S[["percent.mt"]]  <- PercentageFeatureSet(data_11_S, pattern = "^MT-")
data_11_S[["percent.rbp"]] <- PercentageFeatureSet(data_11_S, pattern = "^RP[SL]")
data_11_S[["group"]] <- "11_S"
data_11_S[["compartment"]] <- "skin"

data_12_S[["percent.mt"]]  <- PercentageFeatureSet(data_12_S, pattern = "^MT-")
data_12_S[["percent.rbp"]] <- PercentageFeatureSet(data_12_S, pattern = "^RP[SL]")
data_12_S[["group"]] <- "12_S"
data_12_S[["compartment"]] <- "skin"

data_13_S[["percent.mt"]]  <- PercentageFeatureSet(data_13_S, pattern = "^MT-")
data_13_S[["percent.rbp"]] <- PercentageFeatureSet(data_13_S, pattern = "^RP[SL]")
data_13_S[["group"]] <- "13_S"
data_13_S[["compartment"]] <- "skin"

data_14_S[["percent.mt"]]  <- PercentageFeatureSet(data_14_S, pattern = "^MT-")
data_14_S[["percent.rbp"]] <- PercentageFeatureSet(data_14_S, pattern = "^RP[SL]")
data_14_S[["group"]] <- "14_S"
data_14_S[["compartment"]] <- "skin"

data_healthy[["percent.mt"]]  <- PercentageFeatureSet(data_healthy, pattern = "^MT-")
data_healthy[["percent.rbp"]] <- PercentageFeatureSet(data_healthy, pattern = "^RP[SL]")
data_healthy[["group"]] <- "healthy"
data_healthy[["compartment"]] <- "blood"






# merge into one single seurat object. Add cell ids just in case you have overlapping barcodes between the datasets.
alldata <- merge(data_01_BE, c(data_01_BL,
                               data_01_BLC,
                               data_01_SE,
                               data_01_SL,
                               data_02_BE,
                               data_02_BL,
                               data_02_SE,
                               data_02_SL,
                               data_03_BE,
                               data_03_BL,
                               data_03_SE,
                               data_03_SL,
                               data_04_S,
                               data_06_S,
                               data_07_S,
                               data_08_S,
                               data_09_S,
                               data_10_S,
                               data_11_S,
                               data_12_S,
                               data_13_S,
                               data_14_S,
                               data_healthy), add.cell.ids=c("01_BE",
                                                          "01_BL",
                                                          "01_BLC",
                                                          "01_SE",
                                                          "01_SL",
                                                          "02_BE",
                                                          "02_BL",
                                                          "02_SE",
                                                          "02_SL",
                                                          "03_BE",
                                                          "03_BL",
                                                          "03_SE",
                                                          "03_SL",
                                                          "04_S",
                                                          "06_S",
                                                          "07_S",
                                                          "08_S",
                                                          "09_S",
                                                          "10_S",
                                                          "11_S",
                                                          "12_S",
                                                          "13_S",
                                                          "14_S",
                                                          "healthy"))







metadata <- alldata@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

head(metadata)



View(alldata@meta.data)

my_data

head(my_data)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
my_data[["percent.mt"]] <- PercentageFeatureSet(my_data, pattern = "^MT-")
# my_data$mitoRatio <- my_data@meta.data$percent.mt / 100


#mt.genes <- rownames(my_data)[grep("^MT-",rownames(my_data))]
#C<-GetAssayData(object = my_data, slot = "counts")
#percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
#my_data <- AddMetaData(my_data, percent.mito, col.name = "percent.mito")


rb.genes <- rownames(my_data)[grep("^RP[SL]",rownames(my_data))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
my_data <- AddMetaData(my_data, percent.ribo, col.name = "percent.ribo")




View(alldata@meta.data)


VlnPlot(alldata, features = c("nFeature_RNA", "percent.mt", "percent.rbp"), ncol = 3, group.by = c("compartment"))

table(my_data@meta.data$compartment)


VlnPlot(my_data, features = c("nFeature_RNA", "percent.mito"), pt.size = 1, group.by = "disease") + NoLegend()

FeatureScatter(my_data, feature1="percent.ribo", feature2="percent.mito", group.by = "disease")



View(my_data@meta.data)


# count_table - calculate the number of cells in each cluster
table(my_data@meta.data$compartment)




# Add number of genes per UMI for each cell to metadata
my_data$log10GenesPerUMI <- log10(my_data$nFeature_RNA) / log10(my_data$nCount_RNA)



# Create metadata dataframe
metadata <- my_data@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
meta.data <- meta.data %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)




VlnPlot(my_data, features = c("nFeature_RNA", "percent.mt"), cols = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"))


my_data$definitions <- factor(my_data$definitions,                                    # Change ordering manually
                            levels = c("Healthy", "SingleBystanders", "BystanderGroups", "RelateMainClone", "MainClone"))








pbmc <- subset(my_data, subset = nFeature_RNA > 200 & percent.mt < 25)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(pbmc@meta.data$compartment)



# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

FeatureScatter(my_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(my_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))


# Calculate cell-cycle scores
my_data <- CellCycleScoring(object = my_data, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

# VlnPlot(my_data, features = c("S.Score","G2M.Score"))


# Scater
 sce <- as.SingleCellExperiment(my_data)
# Calture QC-metrics
 sce <- CalculateQCMetrics(sce, feature_controls = list(mito = mt.genes))

 colnames(colData(sce))

 # plotHighestExprs(sce, n = 50, colour_cells_by = "SampleStage", exprs_values = "logcounts")
 
 
# plotHighestExprs(
#   sce,
#   n = 50,
#   colour_cells_by = NULL,
#   drop_features = NULL,
#   exprs_values = "counts",
#   by_exprs_values = exprs_values,
#   feature_names_to_plot = NULL,
#   as_percentage = TRUE
# )
 
 
 # ?retrieveCellInfo
 
 
 
# plot each sample separately
# plotScater(sce, block1 = "ident", nfeatures = 1000)


sce <- runPCA(sce, use_coldata = TRUE,
              detect_outliers = TRUE)





# normalizing the data
my_data <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# extract normalized data
# data.table<- as.matrix(GetAssayData(my_data, slot = "data"))

# getwd()
# setwd("/Users/yun-tsan.chang/CTCL_research/project_scRNA_seq/GeneVia/raw_data/ALL/four_definitions/blood/")
# write.csv(data.table, file="Normalized_all_blood_healthy_raw_counts_TCR_without_Tc_MHCmodified.csv", col.names = TRUE, row.names = TRUE)


my_data <- FindVariableFeatures(my_data, selection.method = "vst", nfeatures = 2000)

# identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(my_data), 100)

# plot1 <- VariableFeaturePlot(my_data)
# LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))
# plot(plot1)


# Scaling the data
all.genes <- rownames(my_data)
my_data <- ScaleData(my_data, features = all.genes)

# run PCA
my_data <- RunPCA(my_data, features = VariableFeatures(object = my_data))

# Examine and visualize PCA results a few different ways
print(my_data[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(my_data, dims = 1:2, reduction = "pca")





DimPlot(my_data, reduction = "pca")


DimHeatmap(my_data, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(my_data, dims = 1:2, cells = 500, balanced = TRUE)


# tSNE
my_data <- RunTSNE(my_data) 
DimPlot(object = my_data, reduction = "tsne", group.by = "definitions" ,cols = c("dodgerblue2", "red2", "deeppink1", "blue3"))

View(my_data@meta.data)

DimPlot(object = my_data, reduction = "tsne", group.by = "orig.ident") + ggtitle("WaG_skin+blood_tumor")


DimPlot(object = my_data, reduction = "tsne", group.by = "StageNow")

my_data <- RunTSNE(my_data, dims = 1:20)
DimPlot(object = my_data, reduction ="tsne", group.by = "group")
DimPlot(object = my_data, reduction = "tsne", group.by = "definitions") + ggtitle("patients")


M.t.cells <- subset(my_data, idents = "")



FeaturePlot(my_data, reduction = "tsne", features = c("IL32", "HLA-A", "HLA-B", "HLA-C"))



RidgePlot(my_data, features = c("IL32"), group.by = "definitions",ncol = 2, cols = c("orchid2", "green3", "red2", "tan1", "slateblue3"))
RidgePlot(my_data, features = c("MKI67"), ncol = 2)
RidgePlot(my_data, features = c("ABCB5"), cols = c("orchid2", "green3", "red2", "tan1", "slateblue3"), ncol = 2)


RidgePlot(my_data, features = c("CD52", "PTPRC", "IL2RG", "S100A4", "IL32", "HLA-A", "HLA-B", "HLA-C"), ncol = 4)
RidgePlot(my_data, features = c("CXCR4", "CCR7", "TNFRSF8", "IKZF2"), group.by = "clone", ncol = 2)
RidgePlot(my_data, features = c("THBS1", "SPINK6", "MUC16", "PLA2G2D", "NPTX1", "AMZ1", "NTRK1", "IL6"), group.by = "definitions" ,ncol = 2)



plot1 <- RidgePlot(my_data, features = "NTRK1")
plot2 <- RidgePlot(my_data, features = "NTRK1", group.by = "definitions")
CombinePlots(plots = list (plot1, plot2))


# ordering the identity
table(my_data@meta.data$definitions)
my_data$definitions <- factor(x = my_data$definitions, levels = c("Healthy", "SingleBystanders", "BystanderGroups", "RelateMainClone", "MainClone"))

VlnPlot(my_data, features = c("HLA-A"), group.by = "definitions", cols = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"), ncol = 5)

VlnPlot(my_data, features = c("HLA-A", "IFITM1", "IFITM2", "SERF2"), group.by = "orig.ident", ncol = 4)

VlnPlot(my_data, features = c("IL32", "HLA-A", "HLA-B", "HLA-C"), group.by = "definitions", ncol = 4)



FeatureScatter(my_data, feature1 = "rna_IL32", feature2 = "HLA-A")


plot4 <- VlnPlot(my_data, features = "NTRK1", group.by = "definitions")
CombinePlots(plots = list (plot3, plot4))


FeaturePlot(my_data, features = "NTRK1")



CombinePlots(plots = list(plot1, plot2, plot3))




VlnPlot(my_data, features = c(
   "BEST1", 
   "RNF213",
   "EIF4A2",
   "HLA-A",
   "IFITM1",
   "IFITM2",
   "RGS1",
   "SERF2",
   "SQSTM1",
   "TOMM20",
   "PGK1",
   "HLA-E",
   "PLIN2",
   "RBM5",
   "ZNF500",
   "S100A11",
   "FYCO1",
   "ACTG1",
   "ANXA2",
   "AC005670.2",
   "ESYT2",
   "MT.RNR1",
   "DNAJB6",
   "KLF6",
   "CCT6A",
   "MT.RNR2",
   "SAT1",
   "C4orf51",
   "MRNIP",
   "HNRNPR",
   "FTH1",
   "MT-CO1",
   "TRBC2",
   "SLC2A3",
   "LGALS8",
   "HNRNPH1",
   "DDX5",
   "UBC",
   "VDAC3",
   "B2M.1",
   "HSPA8",
   "ABCC6",
   "LCP1",
   "TMSB4X",
   "TSPAN13",
   "EEF1A1",
   "CD96",
   "MYL6"
   
) ,group.by = "definitions", cols = c("dodgerblue2", "red2", "deeppink1", "blue3"), ncol = 10)


VlnPlot(my_data, features = c(
   "BEST1", 
   "RNF213",
   "EIF4A2",
   "HLA-A",
   "IFITM1",
   "IFITM2",
   "RGS1",
   "SERF2",
   "SQSTM1",
   "TOMM20",
   "PGK1",
   "HLA-E",
   "PLIN2",
   "RBM5",
   "ZNF500",
   "S100A11",
   "FYCO1",
   "ACTG1",
   "ANXA2",
   "AC005670.2"

) ,group.by = "orig.ident",ncol = 10)



plot2 <- RidgePlot(my_data, features = "THBS1", ncol = 2)






FeatureScatter(object = my_data, feature1 = "IL32", feature2 = "HLA.ABC", group.by = "disease")

FeatureScatter(object = my_data, feature1 = "IL32", feature2 = "HLA-A", cols = c("orchid2", "slateblue3")) + ggtitle("skin_clone 0.5")










VlnPlot(my_data, features = c("CXCR4"), group.by = "compartment_tumor")


VlnPlot(my_data, features = c("TNFAIP3", "ATF3", "IL4I1", "IRAK2", "IRAK3", "IRAK4", "SOCS1", "TRAC"), group.by = "disease", ncol = 4, cols = c("blue", "red"))

VlnPlot(my_data, features = c("ACTB", "SDHA", "YWHAZ"), group.by = "disease", ncol = 3)








FeaturePlot(my_data, reduction = "tsne", features = c("IL32", "HLA-A", "HLA-B", "HLA-C"))








VlnPlot(my_data, features = c("IL32", "HLA-A", "HLA-B", "HLA-C"), group.by = "tumor")

FeaturePlot(my_data, features = c("TNFAIP3", "ATF3", "IL4I1", "IRAK2", "IRAK3", "IRAK4", "SOCS1", "TRAC"), ncol=4)

FeaturePlot(my_data, features = c("ACTB", "SDHA", "YWHAZ"), ncol=3)




# find all markers of cluster 1
cluster1.markers <- FindMarkers(my_data, ident.1 = "MainClone", min.pct = 0.25)
head(cluster1.markers, n = 5)

# Find all markers of cluster 1 from cluster 0
cluster1.markers <- FindMarkers(my_data, ident.1 = 1, ident.2 = 0, min.pct = 0.25, group.by = "compartment")

# view results
head(cluster1.markers)

setwd("/Users/yun-tsan.chang/CTCL_research/project_scRNA_seq/GeneVia/raw_data/ALL/four_definitions/skin and blood/riched genes in each cluster compared to healthy/")
write.csv(cluster1.markers, file="skin_early_vs_late.markers.csv", col.names = TRUE, row.names = TRUE)



# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster1.markers <- FindMarkers(my_data, ident.1 = c("MainClone"), ident.2 = c("SingleBystanders"), min.pct = 0.25)
head(cluster1.markers, n = 80)

top10 <- cluster1.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(my_data, features = top10$gene) + NoLegend()




mf.markers <- FindAllMarkers(my_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 


mf.markers %>%
   group_by(cluster) %>%
   top_n(n = 2, wt = avg_log2FC)



mf.markers %>%
   group_by(cluster) %>%
   top_n(n = 30, wt = avg_log2FC) -> top10
DoHeatmap(my_data, features = top10$gene, group.colors = c("blue", "red")) + NoLegend()



?DoHeatmap







 # cluster the cells
my_data <- FindNeighbors(my_data, dims = 1:10)
my_data <- FindClusters(my_data, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(my_data), 5)

cluster_IDs <- Idents(my_data)
write.csv(cluster_IDs, file="only_early_clone_cluster_IDs.csv", col.names = TRUE, row.names = TRUE)






# UMAP
my_data <- RunUMAP(my_data, dims = 1:20)

DimPlot(my_data, reduction = "umap", group.by = "StageNow")
DimPlot(my_data, reduction = "umap")
plot_grid(p1)
plot_grid(p2)


