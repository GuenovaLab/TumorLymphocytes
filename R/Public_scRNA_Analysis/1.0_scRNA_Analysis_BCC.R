#################################################################################
# Analysis of cancerous T lymphocytes 
########################################

library(Seurat)
library(scater)
library(dplyr)
library(DoubletFinder)

sample_color = data.frame(
  sample = c(
    "HC_P112",
    "HC_P115",
    "HC_P116",
    "HC_P121",
    "MF303_skin",
    "MF309_followup",
    "MF309_nonlesional",
    "MF309_thick",
    "MF309_thin",
    "MF311_followup",
    "MF311_thick",
    "MF311_thin",
    "MF312_followup",
    "MF312_thick",
    "MF312_thin",
    "BCC3",
    "BCC4",
    "BCC5",
    "AD74",
    "AD75",
    "AD77",
    "AD81",
    "AD96"
  ),
  "sample_color" = c(
    "grey5",
    "grey25",
    "grey60",
    "grey80",
    
    "#765aa4ff",
    "#aa3aecff",
    "#984d9dff",
    "#552e8aff",
    "#566f99ff",
    "#ad2070ff",
    
    "#fa3090ff",
    "#ca4060ff",
    "#660a17ff",   
    "#a51d23ff",
    "#651a43ff",
    
    "#3D550C",
    "#81B622",
    "#ECF87F",
    
    "#fff9ae",
    "#f8ed62",
    "#e9d700",
    "#dab600",
    "#a98600"
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
BCC_seurat = qs::qread(file.path(output, "merged_seurat.qs"))
BCC_seurat = BCC_seurat[,BCC_seurat$study == "GSE181907_BCC"]

output = file.path(output, "BCC")
if(!dir.exists(output)) dir.create(output)
gc()

# Run standard QC

# Run standard Seurat processing
BCC_seurat <- NormalizeData(BCC_seurat)
BCC_seurat <- FindVariableFeatures(BCC_seurat, selection.method = "vst", nfeatures = 2000)
BCC_seurat <- ScaleData(BCC_seurat)
BCC_seurat <- RunPCA(BCC_seurat, npcs = 30)
BCC_seurat <- FindNeighbors(BCC_seurat, dims = 1:30)
BCC_seurat <- RunUMAP(BCC_seurat, dims = 1:30)
BCC_seurat <- FindClusters(BCC_seurat, resolution = 0.8)


sweep.res <- paramSweep_v3(BCC_seurat, PCs = 1:30, sct = FALSE, num.cores = 4)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(BCC_seurat$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075* nrow(BCC_seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
BCC_seurat <- doubletFinder_v3(BCC_seurat, PCs = 1:10, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$MeanBC)])),
                              nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE, )

BCC_seurat$doublet = BCC_seurat@meta.data[,ncol(BCC_seurat@meta.data)]

print(DimPlot(BCC_seurat, group.by = "doublet"))

BCC_seurat = BCC_seurat[, which(BCC_seurat$doublet != "Doublet")]

# Recluster cluster 10
BCC_seurat_T_cells = BCC_seurat[,BCC_seurat$seurat_clusters == 9]
BCC_seurat_T_cells <- FindNeighbors(BCC_seurat_T_cells, reduction = "harmony", dims = 1:30)
BCC_seurat_T_cells <- FindClusters(BCC_seurat_T_cells, resolution = 0.1)
BCC_seurat_T_cells$seurat_clusters = as.numeric(BCC_seurat_T_cells$seurat_clusters)
BCC_seurat_T_cells$seurat_clusters[BCC_seurat_T_cells$seurat_clusters == 1] = 39
BCC_seurat_T_cells$seurat_clusters[BCC_seurat_T_cells$seurat_clusters == 2] = 40
BCC_seurat_T_cells$seurat_clusters[BCC_seurat_T_cells$seurat_clusters == 3] = 41
DimPlot(BCC_seurat_T_cells, group.by = "seurat_clusters")
FeaturePlot(BCC_seurat_T_cells, features = "CD8B")

BCC_seurat$seurat_clusters = as.numeric(BCC_seurat$seurat_clusters)
BCC_seurat$seurat_clusters[colnames(BCC_seurat) %in% colnames(BCC_seurat_T_cells)] = BCC_seurat_T_cells$seurat_clusters
BCC_seurat$seurat_clusters = as.factor(BCC_seurat$seurat_clusters )
DimPlot(BCC_seurat, group.by = "seurat_clusters")

# Add metadata
BCC_seurat$sample_id_color = sample_color$sample_color[match(BCC_seurat$sample_id, sample_color$sample)]

# Perform additional downstream analysis as needed
pdf(file.path(output, paste0("scRNA.pdf")))
print(DimPlot(BCC_seurat, reduction = "umap", group.by = "seurat_clusters"))
print(DimPlot(BCC_seurat, reduction = "umap", group.by = "sample_id", cols = unique(BCC_seurat$sample_id_color[order(BCC_seurat$sample_id)])))
print(FeaturePlot(BCC_seurat,  features = "nCount_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
print(FeaturePlot(BCC_seurat,  features = "nFeature_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
dev.off()

pdf(file.path(output, paste0("UMAP_per_clust.pdf")), width = 6, height = 5)
for(clust in unique(BCC_seurat$seurat_clusters)){
  BCC_seurat. = BCC_seurat
  BCC_seurat.$is_here = FALSE
  BCC_seurat.$is_here[BCC_seurat$seurat_clusters == clust] = TRUE
  tplot = DimPlot(BCC_seurat., reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
                  group.by = "is_here", order = TRUE) + ggtitle(clust)
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
}
dev.off()

# Save the Seurat objects
qs::qsave(BCC_seurat, file.path(output, "BCC_seurat.qs"))

cell_markers = c("PTPRC", "CD3E","CD4","CD8A", "FOXP3", "CD8B","TIGIT","IL2RA", 
                 "LAMP3","CD40",  "CD14", "CD16", "ITGAM", "HLA-DRA","HLA-DRQA1","HLA-DRQB1", "HLA-DRB1", "CD163", "CD1A", "CD207",
                 "MS4A2", "MS4A1",
                 "CDH5","PECAM1","CD34","KRT5","KRT10","KRT14","KRT1", "FABP5", "SPPR2G", "CDSN", 
                 "ITGA6", "ITGB1", "GRHL1", 
                 "CD200", "SOX9", "KRT19", "APOC1", "ACSL5", "ABCC3",
                 "LUM", "PDGFRB", "COL1A1", "ACTA2",
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
  if(i %in% rownames(BCC_seurat)) print( FeaturePlot(BCC_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()

pdf(file.path(output, paste0("Markers_expression_NK.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(NK_signature)){
  print(i)
  if(i %in% rownames(BCC_seurat)) print( FeaturePlot(BCC_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()


png(file.path(output,paste0("UMAP_sample_All.png")), width = 1400, height = 1200, res = 250)
tplot = DimPlot(BCC_seurat, reduction = "umap", pt.size = 1.5, 
                cols = unique(BCC_seurat$sample_id_color[order(BCC_seurat$sample_id)]),
                group.by = "sample_id") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

# Annotation
BCC_seurat. = ScaleData(BCC_seurat, features = cell_markers)

png(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.png")), width = 3200, height = 1600, res = 200)
Seurat::DoHeatmap(BCC_seurat., features = unique(cell_markers), group.by = "seurat_clusters")
dev.off()


pdf(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.pdf")), width = 16, height = 8)
Seurat::DoHeatmap(BCC_seurat., features =  unique(cell_markers), group.by = "seurat_clusters")
dev.off()

png(file.path(output, paste0("DotPlot_cell_markers_seurat_clusters.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(BCC_seurat, features = unique(cell_markers), group.by = "seurat_clusters", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

png(file.path(output, paste0("DotPlot_gold_markers_seurat_clusters.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(BCC_seurat, features = unique(gold_markers), group.by = "seurat_clusters", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()


# Ucell signature NK
BCC_seurat = UCell::AddModuleScore_UCell(BCC_seurat, list("NK" = NK_signature))

Idents(BCC_seurat) = BCC_seurat$celltype
png(file.path(output,paste0("UMAP_NK_signature.png")), width = 1400, height = 1400, res = 150)
FeaturePlot(BCC_seurat, features = "NK_UCell", label = F, reduction = "umap", cols = c("grey", "darkred"))
dev.off()

# Absolute number of NK cells per sample
library(dittoSeq)
BCC_seurat$condition = gsub("[0-9]*","",gsub("_.*", "",BCC_seurat$sample_id))

png(file.path(output, "Sample x Cluster.png"), width = 3000, height = 1600, res = 300)
d = dittoBarPlot(BCC_seurat, "sample_id", group.by = "seurat_clusters",
                 main = "Sample x Cluster", color.panel = unique(BCC_seurat$sample_id_color[order(BCC_seurat$sample_id)]))  
d
dev.off()

res = IDclust::differential_edgeR_pseudobulk_LRT(BCC_seurat,
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

View(res %>% filter(cluster == 20))
WriteXLS::WriteXLS(res, ExcelFileName = file.path(output, "One_vs_All_clusters.xlsx"),
                   SheetNames = "One_vs_All")

# Annotate with CellTypes

map_cluster_celltype = unlist(list(
  "Epidermis_" = c(6, 9, 2, 24, 25, 18) ,
  "Follicular epidermis_" = c(),
  "B cells_" = c(20, 21),
  "T helper cells_" = c(0, 8),
  "Mast cells_" = c(19),
  "Natural Killer cells_" = c(16),
  "Cytotoxic T cells_" = c(5),
  "Macrophages and monocytes_" = c(15, 23),
  "Dendritic cells_" = c(22),
  "Fibroblasts_" = c(12),
  "Melanocytes and Schwann cells_" = c(3),
  "Endothelial cells_" = c(11),
  "Smooth muscle cells_" = c(17),
  "Tumour cells_" = c(13,10,4,14,1,7),
  "Platelets_"  = c(),
  "Gamma Delta T cells_"  = c(),
  "Lymphatic Endothelial cells_"  = c()
))


BCC_seurat$celltype = gsub("_.*", "", names(map_cluster_celltype))[match(BCC_seurat$seurat_clusters, map_cluster_celltype)]

col_to_celltype = read.csv("output/col_to_celltype.csv")
BCC_seurat$celltype_color = col_to_celltype$x[match(BCC_seurat$celltype, col_to_celltype$X)]

png(file.path(output,paste0("UMAP_cell_type_All.png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(BCC_seurat, reduction = "umap", pt.size = 0.75,
                cols = c(unique(BCC_seurat$celltype_color[order(BCC_seurat$celltype)])),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

for(gene in NK_signature){
  png(file.path(output,paste0("Violin_cell_type_",gene,".png")), width = 1800, height = 1200, res = 250)
  print(VlnPlot(BCC_seurat, features = gene, group.by = "celltype",
                cols =c(unique(BCC_seurat$celltype_color[order(BCC_seurat$celltype)])) ))
  dev.off()
}


BCC_seurat. = BCC_seurat
BCC_seurat.$celltype = factor(BCC_seurat.$celltype, levels = unique(gsub("_.*", "", names(map_cluster_celltype))))
png(file.path(output, paste0("DotPlot_gold_markers_celltype.png")), width = 2400, height = 1800, res = 250)
Seurat::DotPlot(BCC_seurat., features = gold_markers, group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

# Differential Analysis
library(IDclust)
list_res = list()

res = IDclust::differential_edgeR_pseudobulk_LRT(BCC_seurat,
                                                 by = "celltype",
                                                 assay = "RNA",
                                                 logFC.th = 0.5, qval.th = 0.01, min.pct = 0.2)

WriteXLS::WriteXLS(res, ExcelFileName = file.path(output, "One_vs_All_celltypes.xlsx"),
                   SheetNames = "One_vs_All")


library(enrichR)
res$cluster_of_origin = "Alpha"
pathway_df = top_enriched_pathways(
  res,
  top = 25,
  gene_col = "gene",
  qval.th = 0.1, 
  max_genes_enriched = 1000
)
WriteXLS::WriteXLS(pathway_df, ExcelFileName = file.path(output, "Pathways_celltypes.xlsx"),
                   SheetNames = "Pathways")


for(gene in c("HLA-A", "HLA-B", "HLA-C", "HLA-E")){
  
  png(file.path(output,paste0("Violin_cell_type_",gene,".png")), width = 1800, height = 1200, res = 250)
  print(VlnPlot(BCC_seurat, features = gene, group.by = "celltype",
                cols =c(unique(BCC_seurat$celltype_color[order(BCC_seurat$celltype)])) ))
  dev.off()
  
  png(file.path(output,paste0("UMAP_",gene,".png")), width = 1800, height = 1200, res = 250)
  print( FeaturePlot(BCC_seurat,  cols = c("grey", "darkred"), features = gene)) 
  dev.off()
}

png(file.path(output, paste0("DotPlot_HLAs_celltype.png")), width = 2400, height = 1800, res = 250)
Seurat::DotPlot(BCC_seurat., features = c("HLA-A", "HLA-B", "HLA-C", "HLA-E"), group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()


BCC_seurat. = BCC_seurat[, BCC_seurat$celltype %in% c("Cytotoxic T cells", 
                                                      "Dendritic cells",
                                                      "Gamma Delta T cells",
                                                      "Macrophages and monocytes",
                                                      "Natural Killer cells",
                                                      "T helper cells",
                                                      "B cells",
                                                      "Gamma Delta T cells",
                                                      "Mast cells"
)]

png(file.path(output, "CellType x Sample x Condition Immune.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(BCC_seurat., "celltype", group.by = "sample_id",
             color.panel = unique(BCC_seurat.$celltype_color[order(BCC_seurat.$celltype)]),
             main = "CellType x Sample x Condition")  + xlab("")
dev.off()


BCC_seurat. = BCC_seurat
BCC_seurat.$is_NK = ifelse(BCC_seurat.$celltype == "Natural Killer cells", "Yes", "No")
d =  dittoBarPlot(BCC_seurat., "is_NK", group.by = "sample_id",
                  color.panel = sample_color$sample_color[match(levels(d$data$label),
                                                                sample_color$sample)],
                  main = "Natural Killer Cells",  scale = "count")  + xlab("")

png(file.path(output, "CellType x NK number.png"), width = 3000, height = 1600, res = 300)
df = d$data
df %>% filter(label == "Yes") %>%
  ggplot(aes(x = grouping, y = count, fill = grouping)) + geom_bar(stat = "identity") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = unique(
    BCC_seurat.$sample_id_color)) + ggtitle("Natural Killer cells") + ylab("% Natural Killer cells") + xlab("")
dev.off()

png(file.path(output, "CellType x NK percent.png"), width = 3000, height = 1600, res = 300)
df = d$data %>% mutate(percent = 100 * percent) 
df %>% filter(label == "Yes") %>%
  ggplot(aes(x = grouping, y = percent, fill = grouping)) + geom_bar(stat = "identity") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = unique(
    BCC_seurat.$sample_id_color)) + ggtitle("Natural Killer cells") + ylab("% Natural Killer cells") + xlab("")
dev.off()

d = dittoBarPlot(BCC_seurat, "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  
png(file.path(output, "Sample x CellType x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(BCC_seurat, "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()


png(file.path(output, "Sample x CellType x Condition Immune Compartment Count.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(BCC_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

BCC_seurat. = BCC_seurat[,grep(group, BCC_seurat$study)]

png(file.path(output,paste0("UMAP_celltype_", group,".png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(BCC_seurat., reduction = "umap", pt.size = 0.75,
                cols =  unique(BCC_seurat.$celltype_color[order(BCC_seurat.$celltype)]),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
print(tplot)
dev.off()


png(file.path(output, "Condition x Clonality x Cell Type Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(BCC_seurat., "is_clonal", group.by = "celltype", color.panel = c("grey", "darkred"),
             main = "Condition x Clonality x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()


png(file.path(output, "Sample x Clonality x Condition Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(BCC_seurat., "is_clonal", group.by = "sample_id", color.panel = c("grey", "darkred"),
             main = "Condition x Cell Type Immune", scale = "percent",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()




pdf(file.path(output, paste0("Cell_type_Markers_expression_violin.pdf")), height = 10,
    width = 15, pointsize = 12)
for(i in cell_markers){
  print( VlnPlot(BCC_seurat,  features = i,  pt.size = 1.5, group.by = "seurat_clusters"))
}
dev.off()

# Save the Seurat objects
qs::qsave(BCC_seurat, file.path(output, "BCC_seurat.qs"))




# 
