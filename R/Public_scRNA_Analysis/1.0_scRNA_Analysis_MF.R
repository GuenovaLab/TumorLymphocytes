#################################################################################
# Analysis of cancerous T lymphocytes 
########################################

library(Seurat)
library(scater)
library(dplyr)

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

# Merge all Seurat objects
MF_seurat = qs::qread(file.path(output, "merged_seurat.qs"))
MF_seurat = MF_seurat[,MF_seurat$study == "GSE173205__GSE165623_MF_TCR"]

output = file.path(output, "MF")
if(!dir.exists(output)) dir.create(output)
gc()

# Run standard QC

# Run standard Seurat processing
MF_seurat <- NormalizeData(MF_seurat)
MF_seurat <- FindVariableFeatures(MF_seurat, selection.method = "vst", nfeatures = 2000)
MF_seurat <- ScaleData(MF_seurat)
MF_seurat <- RunPCA(MF_seurat, npcs = 30)
MF_seurat$study = "GSE173205"
MF_seurat$study[MF_seurat$sample_id != "MF303_skin"] = "GSE165623"
MF_seurat <- RunHarmony(MF_seurat, group.by.vars = "study")
MF_seurat <- FindNeighbors(MF_seurat, reduction = "harmony", dims = 1:30)
MF_seurat <- RunUMAP(MF_seurat, reduction = "harmony", dims = 1:30)
MF_seurat <- FindClusters(MF_seurat, resolution = 0.8)

# Recluster cluster 10
MF_seurat_T_cells = MF_seurat[,MF_seurat$seurat_clusters == 9]
MF_seurat_T_cells <- FindNeighbors(MF_seurat_T_cells, reduction = "harmony", dims = 1:30)
MF_seurat_T_cells <- FindClusters(MF_seurat_T_cells, resolution = 0.1)
MF_seurat_T_cells$seurat_clusters = as.numeric(MF_seurat_T_cells$seurat_clusters)
MF_seurat_T_cells$seurat_clusters[MF_seurat_T_cells$seurat_clusters == 1] = 39
MF_seurat_T_cells$seurat_clusters[MF_seurat_T_cells$seurat_clusters == 2] = 40
MF_seurat_T_cells$seurat_clusters[MF_seurat_T_cells$seurat_clusters == 3] = 41
DimPlot(MF_seurat_T_cells, group.by = "seurat_clusters")
FeaturePlot(MF_seurat_T_cells, features = "CD8B")

MF_seurat$seurat_clusters = as.numeric(MF_seurat$seurat_clusters)
MF_seurat$seurat_clusters[colnames(MF_seurat) %in% colnames(MF_seurat_T_cells)] = MF_seurat_T_cells$seurat_clusters
MF_seurat$seurat_clusters = as.factor(MF_seurat$seurat_clusters )
DimPlot(MF_seurat, group.by = "seurat_clusters")

# Add metadata
MF_seurat$sample_id_color = sample_color$sample_color[match(MF_seurat$sample_id, sample_color$sample)]
MF_seurat$has_TCR = ifelse(is.na(MF_seurat$TCR_raw_consensus_id) | MF_seurat$TCR_raw_consensus_id == "None", FALSE, TRUE)
clones = table(MF_seurat$TCR_raw_consensus_id)
clones = clones[clones>100]
MF_seurat$is_clonal <- ifelse(MF_seurat$TCR_raw_consensus_id %in% names(clones), TRUE, FALSE)

# Perform additional downstream analysis as needed
pdf(file.path(output, paste0("scRNA.pdf")))
print(DimPlot(MF_seurat, reduction = "umap", group.by = "seurat_clusters"))
print(DimPlot(MF_seurat, reduction = "umap", group.by = "sample_id", cols = unique(MF_seurat$sample_id_color[order(MF_seurat$sample_id)])))
print(FeaturePlot(MF_seurat,  features = "nCount_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
print(FeaturePlot(MF_seurat,  features = "nFeature_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
DimPlot(MF_seurat, reduction = "umap", pt.size = 0.2, cols = c( "grey",  "darkred"),
        group.by = "has_TCR") + ggtitle("")
DimPlot(MF_seurat, reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
        group.by = "is_clonal", ) + ggtitle("")
dev.off()

png(file.path(output,paste0("UMAP_has_TCR.png")), width = 1400, height = 1200, res = 250)
DimPlot(MF_seurat, reduction = "umap", pt.size = 0.2, cols = c( "grey",  "darkred"),
        group.by = "has_TCR") + ggtitle("")
dev.off()

png(file.path(output,paste0("UMAP_is_clonal.png")), width = 1400, height = 1200, res = 250)
DimPlot(MF_seurat, reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
        group.by = "is_clonal", ) + ggtitle("")
 dev.off()


pdf(file.path(output, paste0("UMAP_per_sample.pdf")), width = 9, height = 8)
 
for(samp in unique(MF_seurat$sample_id)){
   MF_seurat. = MF_seurat
   MF_seurat.$is_here = FALSE
   MF_seurat.$is_here[MF_seurat$sample_id == samp] = TRUE
   tplot = DimPlot(MF_seurat., reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
                   group.by = "is_here", order = TRUE) + ggtitle(samp)
   tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
   print(tplot)
 }
dev.off()

pdf(file.path(output, paste0("UMAP_per_clust.pdf")), width = 6, height = 5)
for(clust in unique(MF_seurat$seurat_clusters)){
  MF_seurat. = MF_seurat
  MF_seurat.$is_here = FALSE
  MF_seurat.$is_here[MF_seurat$seurat_clusters == clust] = TRUE
  tplot = DimPlot(MF_seurat., reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
                  group.by = "is_here", order = TRUE) + ggtitle(clust)
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
}
dev.off()

# Save the Seurat objects
qs::qsave(MF_seurat, file.path(output, "MF_seurat.qs"))

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


pdf(file.path(output, paste0("Markers_expression.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(cell_markers)){
  print(i)
  if(i %in% rownames(MF_seurat)) print( FeaturePlot(MF_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()

pdf(file.path(output, paste0("Markers_expression_NK.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(NK_signature)){
  print(i)
  if(i %in% rownames(MF_seurat)) print( FeaturePlot(MF_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()


png(file.path(output,paste0("UMAP_sample_All.png")), width = 1400, height = 1200, res = 250)
tplot = DimPlot(MF_seurat, reduction = "umap", pt.size = 1.5, 
                cols = unique(MF_seurat$sample_id_color[order(MF_seurat$sample_id)]),
                group.by = "sample_id") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

png(file.path(output,paste0("UMAP_sample_Healthy.png")), width = 1400, height = 1200, res = 250)
MF_seurat. = MF_seurat[,MF_seurat$condition == "HC"]
tplot = DimPlot(MF_seurat., reduction = "umap", pt.size = 1.5, 
                cols = unique(MF_seurat.$sample_id_color[order(MF_seurat.$sample_id)]),
                group.by = "sample_id") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

png(file.path(output,paste0("UMAP_sample_MF.png")), width = 1400, height = 1200, res = 250)
MF_seurat. = MF_seurat[,MF_seurat$condition == "MF"]
tplot = DimPlot(MF_seurat, reduction = "umap", pt.size = 1.5, 
                cols = unique(MF_seurat$sample_id_color[order(MF_seurat$sample_id)]),
                group.by = "sample_id") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

# Annotation
MF_seurat. = ScaleData(MF_seurat, features = cell_markers)

png(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.png")), width = 3200, height = 1600, res = 200)
Seurat::DoHeatmap(MF_seurat., features = unique(cell_markers), cells = sample(ncol(MF_seurat.), 50000),group.by = "seurat_clusters")
dev.off()


pdf(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.pdf")), width = 16, height = 8)
Seurat::DoHeatmap(MF_seurat., features =  unique(cell_markers), cells = 1:1000, group.by = "seurat_clusters")
dev.off()

png(file.path(output, paste0("DotPlot_cell_markers_seurat_clusters.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(MF_seurat, features = unique(cell_markers), group.by = "seurat_clusters", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()


# Ucell signature NK
MF_seurat = UCell::AddModuleScore_UCell(MF_seurat, list("NK" = NK_signature))

Idents(MF_seurat) = MF_seurat$celltype
png(file.path(output,paste0("UMAP_NK_signature.png")), width = 1400, height = 1400, res = 150)
FeaturePlot(MF_seurat, features = "NK_UCell", label = F, reduction = "umap", cols = c("grey", "darkred"))
dev.off()


png(file.path(output,paste0("UMAP_NK_signature_HC.png")), width = 1400, height = 1400, res = 150)
MF_seurat. = MF_seurat[,MF_seurat$condition == "HC"]
FeaturePlot(MF_seurat., features = "NK_UCell", label = F, reduction = "umap", cols = rev(viridis::magma(100)))
dev.off()

png(file.path(output,paste0("UMAP_NK_signature_MF.png")), width = 1400, height = 1400, res = 150)
MF_seurat. = MF_seurat[,MF_seurat$condition == "HC"]
FeaturePlot(MF_seurat., features = "NK_UCell", label = F, reduction = "umap", cols = c("grey", "darkred"))
dev.off()

# Absolute number of NK cells per sample
library(dittoSeq)
MF_seurat$condition = gsub("[0-9]*","",gsub("_.*", "",MF_seurat$sample_id))

png(file.path(output, "Cluster x Clonality.png"), width = 3000, height = 1600, res = 300)
d = dittoBarPlot(MF_seurat, "is_clonal", group.by = "seurat_clusters",
                 main = "Condition x Cell Type", color.panel = c("grey", "darkred"))  
d
dev.off()

png(file.path(output, "Sample x Cluster.png"), width = 3000, height = 1600, res = 300)
d = dittoBarPlot(MF_seurat, "sample_id", group.by = "seurat_clusters",
                 main = "Sample x Cluster", color.panel = unique(MF_seurat$sample_id_color[order(MF_seurat$sample_id)]))  
d
dev.off()

# Annotate with CellTypes
map_cluster_celltype = unlist(list(
  "Epidermis_" = c(6, 24,16,36,21,37) ,
  "Follicular epidermis_" = c(34),
  "B cells_" = c(20),
  "T helper cells_" = c(3,33, 26),
  "Natural Killer cells_" = c(40),
  "Cytotoxic T cells_" = c(39, 41),
  "Macrophages and monocytes_" = c(4,7),
  "Dendritic cells_" = c(31, 22),
  "Fibroblasts_" = c(13, 12, 11, 8, 9,32),
  "Melanocytes and Schwann cells_" = c(29),
  "Endothelial cells_" = c(5, 27,25),
  "Smooth muscle cells_" = c(30, 17 ),
  "Tumour cells_" = c(1,2,28,19,18,23,35,38),
  "Platelets_"  = c(14),
  "Lymphatic Endothelial cells_"  = c(15)
))


MF_seurat$celltype = gsub("_.*", "", names(map_cluster_celltype))[match(MF_seurat$seurat_clusters, map_cluster_celltype)]

col_to_celltype = read.csv("output/col_to_celltype.csv")
MF_seurat$celltype_color = col_to_celltype$x[match(MF_seurat$celltype, col_to_celltype$X)]

png(file.path(output,paste0("UMAP_cell_type_All.png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(MF_seurat, reduction = "umap", pt.size = 0.75,
                cols = c(unique(MF_seurat$celltype_color[order(MF_seurat$celltype)])),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

for(gene in NK_signature){
  png(file.path(output,paste0("Violin_cell_type_",gene,".png")), width = 1800, height = 1200, res = 250)
  print(VlnPlot(MF_seurat, features = gene, group.by = "celltype",
                cols =c(unique(MF_seurat$celltype_color[order(MF_seurat$celltype)])) ))
  dev.off()
}

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

MF_seurat. = MF_seurat
MF_seurat.$celltype = factor(MF_seurat.$celltype, levels = unique(gsub("_.*", "", names(map_cluster_celltype))))
png(file.path(output, paste0("DotPlot_gold_markers_celltype.png")), width = 2400, height = 1800, res = 250)
Seurat::DotPlot(MF_seurat., features = gold_markers, group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()


# Differential Analysis
library(IDclust)
list_res = list()

res = IDclust::differential_edgeR_pseudobulk_LRT(MF_seurat,
                                                 by = "celltype",
                                                 assay = "RNA",
                                                 biological_replicate_col = "sample_id",
                                                 logFC.th = 0.5, qval.th = 0.01, min.pct = 0.2)

topmarkers = IDclust::top_differential_markers(
  res,
  top = 10,
  gene_col = "gene",
  logFC_col = "avg_log2FC",
  qvalue_col = "p_val_adj",
  order_by = "logFC_col",
  pseudogene_pattern = NULL
)
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
  print(VlnPlot(MF_seurat, features = gene, group.by = "celltype",
                cols =c(unique(MF_seurat$celltype_color[order(MF_seurat$celltype)])) ))
  dev.off()
  
  png(file.path(output,paste0("UMAP_",gene,".png")), width = 1800, height = 1200, res = 250)
  print( FeaturePlot(MF_seurat,  cols = c("grey", "darkred"), features = gene)) 
  dev.off()
}

png(file.path(output, paste0("DotPlot_HLAs_celltype.png")), width = 2400, height = 1800, res = 250)
Seurat::DotPlot(MF_seurat., features = c("HLA-A", "HLA-B", "HLA-C", "HLA-E"), group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

d = dittoBarPlot(MF_seurat, "sample_id", group.by = "celltype",
                 color.panel =  unique(MF_seurat$sample_id_color[order(MF_seurat$sample_id)]),
                 main = "Condition x Cell Type")  
png(file.path(output, "CellType x Sample x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(MF_seurat, "sample_id", group.by = "celltype",
             color.panel = sample_color$sample_color[match(levels(d$data$label),
                                                           sample_color$sample)],
             main = "CellType x Sample x Condition")  + xlab("")
dev.off()


MF_seurat. = MF_seurat[, MF_seurat$celltype %in% c("Cytotoxic T cells", 
                                                   "Dendritic cells",
                                                   "Gamma Delta T cells",
                                                   "Macrophages and monocytes",
                                                   "Natural Killer cells",
                                                   "T helper cells"
)]
png(file.path(output, "CellType x Sample x Condition Immune.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(MF_seurat., "celltype", group.by = "sample_id",
             color.panel = unique(MF_seurat.$celltype_color[order(MF_seurat.$celltype)]),
             main = "CellType x Sample x Condition")  + xlab("")
dev.off()


MF_seurat. = MF_seurat
MF_seurat.$is_NK = ifelse(MF_seurat.$celltype == "Natural Killer cells", "Yes", "No")
d =  dittoBarPlot(MF_seurat., "is_NK", group.by = "sample_id",
             color.panel = sample_color$sample_color[match(levels(d$data$label),
                                                           sample_color$sample)],
             main = "Natural Killer Cells",  scale = "count")  + xlab("")

png(file.path(output, "CellType x NK number.png"), width = 3000, height = 1600, res = 300)
df = d$data
df %>% filter(label == "Yes") %>%
  ggplot(aes(x = grouping, y = count, fill = grouping)) + geom_bar(stat = "identity") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = unique(
    MF_seurat.$sample_id_color)) + ggtitle("Natural Killer cells") + ylab("% Natural Killer cells") + xlab("")
dev.off()

png(file.path(output, "CellType x NK percent.png"), width = 3000, height = 1600, res = 300)
df = d$data %>% mutate(percent = 100 * percent) 
df %>% filter(label == "Yes") %>%
  ggplot(aes(x = grouping, y = percent, fill = grouping)) + geom_bar(stat = "identity") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = unique(
    MF_seurat.$sample_id_color)) + ggtitle("Natural Killer cells") + ylab("% Natural Killer cells") + xlab("")
dev.off()


table(MF_seurat$sample_id, MF_seurat$celltype)

groups = c(
  "MF309_thick",
  "MF309_followup",
  "MF311_thick",
  "MF311_followup",
  "MF312_thick",
  "MF312_followup"
)
refs = c(
  "MF309_thin",
  "MF309_thin",
  "MF311_thin",
  "MF311_thin",
  "MF312_thin",
  "MF312_thin"
)

res_list = list()
for(i in 1:length(groups)){
  cat("Comparing ", groups[i], " vs ", refs[i], "...\n")
  MF_seurat. = MF_seurat[,MF_seurat$sample_id %in% c(groups[i], refs[i])]
  MF_seurat. = MF_seurat.[,MF_seurat.$is_clonal]
  MF_seurat. = MF_seurat.[,MF_seurat.$celltype %in% c("Tumour cells", "T helper cells")]
  
  df = IDclust::differential_edgeR_pseudobulk_LRT(MF_seurat.,
                                                   by = "compare",
                                                   assay = "RNA", 
                                                   logFC.th = -1e6, qval.th = 1, min.pct = 0)
  df = df %>% filter(cluster == groups[i])
  df = df %>% mutate(signif = ifelse(p_val_adj < 0.01, ifelse(avg_log2FC > 1, "positive",
                                                              ifelse(avg_log2FC < -1,"negative", "None")), "None"))
  df = df %>% mutate(label = ifelse(signif == "None", "", gene))
  df$p_val_adj[which(df$p_val_adj == 0)] = min(df$p_val_adj[which(df$p_val_adj != 0)])
  genes = c("HLA-E", "HLA-C")
  df = df %>% mutate(label_HLA = ifelse(gene %in% genes, gene, ""))
  
  png(file.path(output, paste0(groups[i], "_vs_", refs[i], "_volcano.png")), width = 2400, height = 2000, res = 300)
  p = ggplot(df, aes(y = -log10(p_val_adj), x = avg_log2FC, color = signif, label = label)) + geom_point() +
    scale_color_manual(values = c( "darkgreen", "black","red")) + theme_classic() +
    geom_vline(xintercept = -1, lty = 5) + 
    geom_vline(xintercept = 1, lty = 5) +
    geom_hline(yintercept = 2, lty = 5) + ggtitle(paste0(groups[i], " vs ", refs[i])) + 
    ggrepel::geom_text_repel() + xlab("Log2 FC") + ylab("-log10(adjusted-pvalue)") +
    ggrepel::geom_text_repel(aes(label = label_HLA), max.overlaps = 1e6)  + theme(legend.position = "none")
  print(p)
  dev.off()
  
  pdf(file.path(output, paste0("Thick_vs_Thin_volcano.pdf")), width = 10, height = 10)
  p = ggplot(df, aes(y = -log10(p_val_adj), x = avg_log2FC, color = signif, label = label)) + geom_point() +
    scale_color_manual(values = c( "darkgreen", "black","red")) + theme_classic() +
    geom_vline(xintercept = -1, lty = 5) + 
    geom_vline(xintercept = 1, lty = 5) +
    geom_hline(yintercept = 2, lty = 5) + ggtitle("Thick vs Thin") + 
    ggrepel::geom_text_repel() + xlab("Log2 FC") + ylab("-log10(adjusted-pvalue)") +
    ggrepel::geom_text_repel(aes(label = label_HLA), max.overlaps = 1e6)  + theme(legend.position = "none")
  print(p)
  dev.off()
  
  res_list[[paste0(groups[i],"_vs_", refs[i])]] = df
}

WriteXLS::WriteXLS(res_list, ExcelFileName = file.path(output, "DA_Tumour_Late_vs_Early.xlsx"),
                   SheetNames = names(res_list))
(df)


d = dittoBarPlot(MF_seurat, "sample_id", group.by = "celltype",
                 color.panel =  unique(MF_seurat$sample_id_color[order(MF_seurat$sample_id)]),
                 main = "Condition x Cell Type")  
png(file.path(output, "CellType x Sample x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(MF_seurat, "sample_id", group.by = "celltype",
             color.panel = sample_color$sample_color[match(levels(d$data$label),
                                                           sample_color$sample)],
             main = "CellType x Sample x Condition")  + xlab("")
dev.off()



d = dittoBarPlot(MF_seurat, "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  
png(file.path(output, "Sample x CellType x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(MF_seurat, "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

MF_seurat. = MF_seurat[,MF_seurat$celltype %in% c(
  "Cytotoxic T cells", 
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
dittoBarPlot(MF_seurat., "celltype", group.by = "sample_id",
             color.panel = unique(MF_seurat.$celltype_color[order(MF_seurat.$celltype)]),
             main = "CellType x Sample x Condition")  + xlab("")
dev.off()

d = dittoBarPlot(MF_seurat., "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  

png(file.path(output, "Sample x CellType x Condition Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(MF_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "percent",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

png(file.path(output, "Sample x CellType x Condition Immune Compartment Count.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(MF_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

MF_seurat. = MF_seurat[,grep(group, MF_seurat$study)]

png(file.path(output,paste0("UMAP_celltype_", group,".png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(MF_seurat., reduction = "umap", pt.size = 0.75,
                cols =  unique(MF_seurat.$celltype_color[order(MF_seurat.$celltype)]),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
print(tplot)
dev.off()


png(file.path(output, "Condition x Clonality x Cell Type Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(MF_seurat., "is_clonal", group.by = "celltype", color.panel = c("grey", "darkred"),
             main = "Condition x Clonality x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()


png(file.path(output, "Sample x Clonality x Condition Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(MF_seurat., "is_clonal", group.by = "sample_id", color.panel = c("grey", "darkred"),
             main = "Condition x Cell Type Immune", scale = "percent",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()




pdf(file.path(seurat_dir_combined, paste0("Cell_type_Markers_expression_violin.pdf")), height = 10,
    width = 15, pointsize = 12)
for(i in cell_markers){
  print( VlnPlot(MF_seurat,  features = i,  pt.size = 1.5, group.by = "seurat_clusters"))
}
dev.off()



# Plot for each object the NK signature
for( i in 1:length(seurat_list)){
  seurat_obj = seurat_list[[i]]
  MF_seurat. = MF_seurat[,MF_seurat$sample_id == seurat_obj$sample_id[1]]
  seurat_obj$NK_Ucell = MF_seurat.$NK_UCell[match(colnames(seurat_obj),  gsub(paste0(".*",seurat_obj$sample_id[1], "_"),"",colnames(MF_seurat.)))]
  
  png(file.path(output,paste0("UMAP_NK_signature_", seurat_obj$sample_id[1],".png")), width = 1400, height = 1400, res = 150)
  print(FeaturePlot(seurat_obj, features = "NK_Ucell", label = F, reduction = "umap", cols = c("grey", "darkred")) + ggtitle(seurat_obj$sample_id[1]))
  dev.off()
  
}

# For MF study compare Tumor T cells vs Healthy T-cells

# For MF study compare Tumor T cells vs non-tumor T cells

# For MF study compare non-tumor T cells vs Healthy T-cells

# Save the Seurat objects
qs::qsave(MF_seurat, file.path(output, "MF_seurat.qs"))




# 
