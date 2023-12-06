#################################################################################
# Analysis of cancerous T lymphocytes 
########################################

library(Seurat)
library(scater)
library(dplyr)
library(DoubletFinder)

sample_color = data.frame(
  sample = c(
    "CBCL"
  ),
  "sample_color" = c(
    "#326BC7"
  )
)

colors = data.frame(
  "cell_type_color" = c(
    "#765aa4ff",
    "#ad2070ff",
    "#fdd62eff",
    "#96c9e6ff",
    "#f48121ff",
    "#68c3a5ff",
    "#ef3e2cff",
    "#0d522aff",
    "#42b649ff",
    "#660a17ff",
    "#3f7fc1ff",
    "#189cd5ff",
    "#0e8342ff",
    "#f9ae1dff",
    "#552e8aff",
    "#8b87c0ff",
    "#984d9dff",
    "#fec964ff",
    "#126badff",
    "#a51d23ff",
    "#e5569fff",
    "#eae71cff"
  )
)

maindir = "/home/localadmin/Documents/scRNA_YunTsan_Revision///"

# Directories -------------------------------------------------------------
output = file.path(maindir, "output");  if(!dir.exists(output)){dir.create(output)}

#################################################################################
# Loading & QC
####################

# Load required libraries
library(Seurat)
library(harmony)
library(DoubletFinder)

# Merge all Seurat objects
# Get the study name
GEOdir_scRNA = file.path("/mnt/RECHERCHE/GUENOVA_LAB/Project_Yun-Tsan_ADCC_paper/dataset to analyze for ADCC paper/") 
studies <- list.files(GEOdir_scRNA, full.names = TRUE)
study_dir = studies[2]
study_name <- basename(study_dir)
print(study_name)

sample = unique(gsub("_clonotypes.*|_filtered_.*|_barcodes.*|_features.*|_matrix.*", "", list.files(study_dir)))

print(sample)

# Read the matrix files
barcode_file <- file.path(study_dir, paste0(sample, "_barcodes.tsv.gz"))
features_file <- file.path(study_dir, paste0(sample, "_features.tsv.gz"))
matrix_file <- file.path(study_dir, paste0(sample, "_matrix.mtx.gz"))

# Load the matrix
matrix <- ReadMtx(mtx = matrix_file, cells = barcode_file, features = features_file)

# Create a Seurat object
CBCL_seurat <- CreateSeuratObject(counts = matrix)
CBCL_seurat$study <- study_name
CBCL_seurat$sample_id <- sample

output = file.path(output, "CBCL")
if(!dir.exists(output)) dir.create(output)
gc()

# Run standard QC

# Run standard Seurat processing
CBCL_seurat <- NormalizeData(CBCL_seurat)
CBCL_seurat <- FindVariableFeatures(CBCL_seurat, selection.method = "vst", nfeatures = 2000)
CBCL_seurat <- ScaleData(CBCL_seurat)
CBCL_seurat <- RunPCA(CBCL_seurat, npcs = 30)
CBCL_seurat <- FindNeighbors(CBCL_seurat, dims = 1:30)
CBCL_seurat <- RunUMAP(CBCL_seurat, dims = 1:30)
CBCL_seurat <- FindClusters(CBCL_seurat, resolution = 0.8)


sweep.res <- paramSweep_v3(CBCL_seurat, PCs = 1:30, sct = FALSE, num.cores = 10)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(CBCL_seurat$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075* nrow(CBCL_seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
CBCL_seurat <- doubletFinder_v3(CBCL_seurat, PCs = 1:10, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$MeanBC)])),
                              nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE, )

CBCL_seurat$doublet = CBCL_seurat@meta.data[,ncol(CBCL_seurat@meta.data)]

print(DimPlot(CBCL_seurat, group.by = "doublet"))

CBCL_seurat = CBCL_seurat[, which(CBCL_seurat$doublet != "Doublet")]

FeaturePlot(CBCL_seurat, features = "CD8B")
DimPlot(CBCL_seurat, group.by = "seurat_clusters")

# Add metadata
CBCL_seurat$sample_id_color = sample_color$sample_color[match(CBCL_seurat$sample_id, sample_color$sample)]

# Perform additional downstream analysis as needed
pdf(file.path(output, paste0("scRNA.pdf")))
print(DimPlot(CBCL_seurat, reduction = "umap", group.by = "seurat_clusters"))
print(DimPlot(CBCL_seurat, reduction = "umap", group.by = "sample_id", cols = unique(CBCL_seurat$sample_id_color[order(CBCL_seurat$sample_id)])))
print(FeaturePlot(CBCL_seurat,  features = "nCount_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
print(FeaturePlot(CBCL_seurat,  features = "nFeature_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
dev.off()

pdf(file.path(output, paste0("UMAP_per_clust.pdf")), width = 6, height = 5)
for(clust in unique(CBCL_seurat$seurat_clusters)){
  CBCL_seurat. = CBCL_seurat
  CBCL_seurat.$is_here = FALSE
  CBCL_seurat.$is_here[CBCL_seurat$seurat_clusters == clust] = TRUE
  tplot = DimPlot(CBCL_seurat., reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
                  group.by = "is_here", order = TRUE) + ggtitle(clust)
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
}
dev.off()

# Save the Seurat objects
qs::qsave(CBCL_seurat, file.path(output, "CBCL_seurat.qs"))

cell_markers = c("PTPRC", "CD3E","CD4","CD8A", "FOXP3", "CD8B","TIGIT","IL2RA", 
                 "LAMP3","CD40",  "CD14", "CD16", "ITGAM", "HLA-DRA","HLA-DRQA1","HLA-DRQB1", "HLA-DRB1", "CD163", "CD1A", "CD207",
                 "MS4A2", "MS4A1",
                 "CDH5","PECAM1","CD34","KRT5","KRT10","KRT14","KRT1", "FABP5", "SPPR2G", "CDSN", 
                 "ITGA6", "ITGB1", "GRHL1", 
                 "CD200", "SOX9", "KRT19", "APOC1", "ACSL5", "ABCC3",
                 "LUM", "PDGFRB", "COL1A1",
                 "SOX10", "S100B", "DCT", "TYRP1", "MLANA")

NK_signature = c("NKG7",
                 "KLRB1",
                 "KLRC1",
                 "KLRD1",
                 "KLRK1",
                 "CD7",
                 "GZMB",
                 "GNLY",
                 "NCAM1",
                 "GZMH",
                 "CCL4",
                 "IFNG",
                 "CCL4L2",
                 "FCGR3B",
                 "FCGR3A"
)

gold_markers = c(
  "KRT6A",
  "KRT19",
  "MS4A1",
  "CD4",
  "NKG7",
  "CD8A",
  "CD14",
  "LAMP3",
  "COL1A1",
  "DCT",
  "CD34",
  "ACTA2",
  "TRBV11-2",
  "THBS1",
  "FLT4"
)

pdf(file.path(output, paste0("Markers_expression.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(cell_markers)){
  print(i)
  if(i %in% rownames(CBCL_seurat)) print( FeaturePlot(CBCL_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()

pdf(file.path(output, paste0("Markers_expression_NK.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(NK_signature)){
  print(i)
  if(i %in% rownames(CBCL_seurat)) print( FeaturePlot(CBCL_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()


png(file.path(output,paste0("UMAP_sample_All.png")), width = 1400, height = 1200, res = 250)
tplot = DimPlot(CBCL_seurat, reduction = "umap", pt.size = 1.5, 
                cols = unique(CBCL_seurat$sample_id_color[order(CBCL_seurat$sample_id)]),
                group.by = "sample_id") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

# Annotation
CBCL_seurat. = ScaleData(CBCL_seurat, features = cell_markers)

png(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.png")), width = 3200, height = 1600, res = 200)
Seurat::DoHeatmap(CBCL_seurat., features = unique(cell_markers), group.by = "seurat_clusters")
dev.off()


pdf(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.pdf")), width = 16, height = 8)
Seurat::DoHeatmap(CBCL_seurat., features =  unique(cell_markers), group.by = "seurat_clusters")
dev.off()

png(file.path(output, paste0("DotPlot_cell_markers_seurat_clusters.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(CBCL_seurat, features = unique(cell_markers), group.by = "seurat_clusters", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

png(file.path(output, paste0("DotPlot_gold_markers_seurat_clusters.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(CBCL_seurat, features = unique(gold_markers), group.by = "seurat_clusters", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()


# Ucell signature NK
CBCL_seurat = UCell::AddModuleScore_UCell(CBCL_seurat, list("NK" = NK_signature))

png(file.path(output,paste0("UMAP_NK_signature.png")), width = 1400, height = 1400, res = 150)
FeaturePlot(CBCL_seurat, features = "NK_UCell", label = F, reduction = "umap", cols = c("grey", "darkred"))
dev.off()


library(dittoSeq)
CBCL_seurat$condition = gsub("[0-9]*","",gsub("_.*", "",CBCL_seurat$sample_id))

png(file.path(output, "Sample x Cluster.png"), width = 3000, height = 1600, res = 300)
d = dittoBarPlot(CBCL_seurat, "sample_id", group.by = "seurat_clusters",
                 main = "Sample x Cluster", color.panel = unique(CBCL_seurat$sample_id_color[order(CBCL_seurat$sample_id)]))  
d
dev.off()

res = IDclust::differential_edgeR_pseudobulk_LRT(CBCL_seurat,
                                                 by = "seurat_clusters",
                                                 assay = "RNA",
                                                 biological_replicate_col = "sample_id",
                                                 logFC.th = 0.5, qval.th = 0.01, min.pct = 0.2)

res$cluster_of_origin = "Alpha"
cell_marker_db = read.csv("/media/localadmin/T7/InstitutCurie/Documents/Data/Annotation/cell_marker_db.csv", header = TRUE)

cell_marker_db = cell_marker_db[grep(c("Epithelium|Blood|Immune|Lymph"), ignore.case = T, cell_marker_db$organ),]
genes = intersect(toupper(res$gene), cell_marker_db$official.gene.symbol)

cell_marker_db = cell_marker_db[which(cell_marker_db$official.gene.symbol %in% genes),]
table(cell_marker_db$cell.type)

res$cell_type = cell_marker_db$cell.type[match(toupper(res$gene), cell_marker_db$official.gene.symbol)]
res = res %>% dplyr::mutate(
  rank_group_activation = order(pct.1),
  rank_logFC = order(avg_log2FC),
  rank_qval = order(p_val_adj),
)

cell_type_markers = res %>% dplyr::filter(!is.na(cell_type)) %>% 
  group_by(cluster_of_origin, cluster, cell_type) %>%
  summarise(total = n(),
            avg_rank_group_activation = mean(rank_group_activation),
            avg_rank_logFC = mean(rank_logFC),
            avg_rank_qval = mean(rank_qval),
            Gene = gene[1]
  ) %>% 
  mutate(
    combined_score_rank = nrow(res) / ((avg_rank_group_activation + avg_rank_logFC + avg_rank_qval) / 3)
  ) %>% slice_max(total, with_ties = TRUE) %>% slice_max(combined_score_rank)
cell_type_markers = cell_type_markers %>% mutate(Term = paste0(cell_type, " | ", total))
cell_type_markers = cell_type_markers[,c("cluster_of_origin", "cluster", "Term")]

WriteXLS::WriteXLS(res, ExcelFileName = file.path(output, "One_vs_All_clusters.xlsx"),
                   SheetNames = "One_vs_All")

# Annotate with CellTypes

map_cluster_celltype = unlist(list(
  "Epidermis_" = c(5, 13, 17) ,
  "Follicular epidermis_" = c(),
  "B cells_" = c(10, 19),
  "T helper cells_" = c(6, 0, 9),
  "Natural Killer cells_" = c(14),
  "Cytotoxic T cells_" = c(1, 7),
  "Macrophages and monocytes_" = c(11),
  "Dendritic cells_" = c(8),
  "Fibroblasts_" = c(4),
  "Melanocytes and Schwann cells_" = c(),
  "Endothelial cells_" = c(3),
  "Smooth muscle cells_" = c(16),
  "Tumour cells_" = c(2, 12),
  "Platelets_"  = c(18),
  "Erythroid cells_"  = c(15),
  "Lymphatic Endothelial cells_"  = c()
))


CBCL_seurat$celltype = gsub("_.*", "", names(map_cluster_celltype))[match(CBCL_seurat$seurat_clusters, map_cluster_celltype)]

col_to_celltype = read.csv("output/col_to_celltype.csv")
CBCL_seurat$celltype_color = col_to_celltype$x[match(CBCL_seurat$celltype, col_to_celltype$X)]

png(file.path(output,paste0("UMAP_cell_type_All.png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(CBCL_seurat, reduction = "umap", pt.size = 0.75,
                cols = c(unique(CBCL_seurat$celltype_color[order(CBCL_seurat$celltype)])),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

for(gene in NK_signature){
  png(file.path(output,paste0("Violin_cell_type_",gene,".png")), width = 1800, height = 1200, res = 250)
  print(VlnPlot(CBCL_seurat, features = gene, group.by = "celltype",
                cols =c(unique(CBCL_seurat$celltype_color[order(CBCL_seurat$celltype)])) ))
  dev.off()
}


CBCL_seurat. = CBCL_seurat
CBCL_seurat.$celltype = factor(CBCL_seurat.$celltype, levels = unique(gsub("_.*", "", names(map_cluster_celltype))))
png(file.path(output, paste0("DotPlot_gold_markers_celltype.png")), width = 2400, height = 1800, res = 250)
Seurat::DotPlot(CBCL_seurat., features = gold_markers, group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

# Differential Analysis
library(IDclust)
list_res = list()

res = IDclust::differential_edgeR_pseudobulk_LRT(CBCL_seurat,
                                                 by = "seurat_clusters",
                                                 assay = "RNA",
                                                 logFC.th = 0.5, qval.th = 0.01, min.pct = 0.2)




res$cluster_of_origin = "Alpha"
cell_marker_db = read.csv("/media/localadmin/T7/InstitutCurie/Documents/Data/Annotation/cell_marker_db.csv", header = TRUE)

cell_marker_db = cell_marker_db[grep(c("Epithelium|Blood|Immune|Lymph"), ignore.case = T, cell_marker_db$organ),]
genes = intersect(toupper(res$gene), cell_marker_db$official.gene.symbol)

cell_marker_db = cell_marker_db[which(cell_marker_db$official.gene.symbol %in% genes),]
table(cell_marker_db$cell.type)

res$cell_type = cell_marker_db$cell.type[match(toupper(res$gene), cell_marker_db$official.gene.symbol)]
res = res %>% dplyr::mutate(
  rank_group_activation = order(pct.1),
  rank_logFC = order(avg_log2FC),
  rank_qval = order(p_val_adj),
)

View(res %>% filter(cluster == 8))

WriteXLS::WriteXLS(res, ExcelFileName = file.path(output, "One_vs_All_clusters.xlsx"),
                   SheetNames = "One_vs_All")


for(gene in c("HLA-A", "HLA-B", "HLA-C", "HLA-E")){
  
  png(file.path(output,paste0("Violin_cell_type_",gene,".png")), width = 1800, height = 1200, res = 250)
  print(VlnPlot(CBCL_seurat, features = gene, group.by = "celltype",
                cols =c(unique(CBCL_seurat$celltype_color[order(CBCL_seurat$celltype)])) ))
  dev.off()
  
  png(file.path(output,paste0("UMAP_",gene,".png")), width = 1800, height = 1200, res = 250)
  print( FeaturePlot(CBCL_seurat,  cols = c("grey", "darkred"), features = gene)) 
  dev.off()
}

gene = "NK_UCell"
png(file.path(output,paste0("Violin_cell_type_",gene,".png")), width = 1800, height = 1200, res = 250)
print(VlnPlot(CBCL_seurat, features = gene, group.by = "celltype",
              cols =c(unique(CBCL_seurat$celltype_color[order(CBCL_seurat$celltype)])) ))
dev.off()

png(file.path(output, paste0("DotPlot_HLAs_celltype.png")), width = 2400, height = 1800, res = 250)
Seurat::DotPlot(CBCL_seurat., features = c("HLA-A", "HLA-B", "HLA-C", "HLA-E"), group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()


CBCL_seurat. = CBCL_seurat[, CBCL_seurat$celltype %in% c("Cytotoxic T cells", 
                                                   "Dendritic cells",
                                                   "Gamma Delta T cells",
                                                   "Macrophages and monocytes",
                                                   "Natural Killer cells",
                                                   "T helper cells",
                                                   "B cells",
                                                   "Gamma Delta T cells"
)]
png(file.path(output, "CellType x Sample x Condition Immune.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(CBCL_seurat., "celltype", group.by = "sample_id",
             color.panel = unique(CBCL_seurat.$celltype_color[order(CBCL_seurat.$celltype)]),
             main = "CellType x Sample x Condition")  + xlab("")
dev.off()


CBCL_seurat. = CBCL_seurat
CBCL_seurat.$is_NK = ifelse(CBCL_seurat.$celltype == "Natural Killer cells", "Yes", "No")
d =  dittoBarPlot(CBCL_seurat., "is_NK", group.by = "sample_id",
                  color.panel = sample_color$sample_color[match(levels(d$data$label),
                                                                sample_color$sample)],
                  main = "Natural Killer Cells",  scale = "count")  + xlab("")

png(file.path(output, "CellType x NK number.png"), width = 3000, height = 1600, res = 300)
df = d$data
df %>% filter(label == "Yes") %>%
  ggplot(aes(x = grouping, y = count, fill = grouping)) + geom_bar(stat = "identity") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = unique(
    CBCL_seurat.$sample_id_color)) + ggtitle("Natural Killer cells") + ylab("% Natural Killer cells") + xlab("")
dev.off()

png(file.path(output, "CellType x NK percent.png"), width = 3000, height = 1600, res = 300)
df = d$data %>% mutate(percent = 100 * percent) 
df %>% filter(label == "Yes") %>%
  ggplot(aes(x = grouping, y = percent, fill = grouping)) + geom_bar(stat = "identity") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = unique(
    CBCL_seurat.$sample_id_color)) + ggtitle("Natural Killer cells") + ylab("% Natural Killer cells") + xlab("")
dev.off()

d = dittoBarPlot(CBCL_seurat, "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  
png(file.path(output, "Sample x CellType x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(CBCL_seurat, "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

CBCL_seurat. = CBCL_seurat[,CBCL_seurat$celltype %in% unique(gsub("_.*","",names(map_cluster_celltype)))[c(3,4,5,6,7,8,9)]]
d = dittoBarPlot(CBCL_seurat., "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  

png(file.path(output, "Sample x CellType x Condition Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(CBCL_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "percent",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

png(file.path(output, "Sample x CellType x Condition Immune Compartment Count.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(CBCL_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()


pdf(file.path(output, paste0("Cell_type_Markers_expression_violin.pdf")), height = 10,
    width = 15, pointsize = 12)
for(i in cell_markers){
  print( VlnPlot(CBCL_seurat,  features = i,  pt.size = 1.5, group.by = "seurat_clusters"))
}
dev.off()

# Save the Seurat objects
qs::qsave(CBCL_seurat, file.path(output, "CBCL_seurat.qs"))



res = IDclust::differential_edgeR_pseudobulk_LRT(CBCL_seurat,
                                                 by = "celltype",
                                                 assay = "RNA",
                                                 logFC.th = 0.5, qval.th = 0.01, min.pct = 0.2)

WriteXLS::WriteXLS(res, ExcelFileName = file.path(output, "One_vs_All_celltypes.xlsx"),
                   SheetNames = "One_vs_All")
