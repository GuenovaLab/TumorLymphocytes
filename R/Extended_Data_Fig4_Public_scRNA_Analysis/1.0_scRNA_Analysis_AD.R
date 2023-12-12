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
AD_seurat = qs::qread(file.path(output, "merged_seurat.qs"))
AD_seurat = AD_seurat[,AD_seurat$study == "GSE222840_AD_TCR"]

output = file.path(output, "AD")
if(!dir.exists(output)) dir.create(output)
gc()

# Run standard QC

# Run standard Seurat processing
AD_seurat <- NormalizeData(AD_seurat)
AD_seurat <- FindVariableFeatures(AD_seurat, selection.method = "vst", nfeatures = 2000)
AD_seurat <- ScaleData(AD_seurat)
AD_seurat <- RunPCA(AD_seurat, npcs = 30)
AD_seurat <- FindNeighbors(AD_seurat, dims = 1:30)
AD_seurat <- RunUMAP(AD_seurat, dims = 1:30)
AD_seurat <- FindClusters(AD_seurat, resolution = 0.8)


sweep.res <- paramSweep_v3(AD_seurat, PCs = 1:30, sct = FALSE, num.cores = 10)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(AD_seurat$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075* nrow(AD_seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
AD_seurat <- doubletFinder_v3(AD_seurat, PCs = 1:10, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$MeanBC)])),
                               nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE, )

AD_seurat$doublet = AD_seurat@meta.data[,ncol(AD_seurat@meta.data)]

print(DimPlot(AD_seurat, group.by = "doublet"))

AD_seurat = AD_seurat[, which(AD_seurat$doublet != "Doublet")]

FeaturePlot(AD_seurat_T_cells, features = "CD8B")
DimPlot(AD_seurat, group.by = "seurat_clusters")

# Add metadata
AD_seurat$sample_id_color = sample_color$sample_color[match(AD_seurat$sample_id, sample_color$sample)]
AD_seurat$has_TCR = ifelse(is.na(AD_seurat$TCR_raw_consensus_id) | AD_seurat$TCR_raw_consensus_id == "None", FALSE, TRUE)
clones = table(AD_seurat$TCR_raw_consensus_id)
clones = clones[clones > 50]
AD_seurat$is_clonal <- ifelse(AD_seurat$TCR_raw_consensus_id %in% names(clones), TRUE, FALSE)

# Perform additional downstream analysis as needed
pdf(file.path(output, paste0("scRNA.pdf")))
print(DimPlot(AD_seurat, reduction = "umap", group.by = "seurat_clusters"))
print(DimPlot(AD_seurat, reduction = "umap", group.by = "sample_id", cols = unique(AD_seurat$sample_id_color[order(AD_seurat$sample_id)])))
print(FeaturePlot(AD_seurat,  features = "nCount_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
print(FeaturePlot(AD_seurat,  features = "nFeature_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
DimPlot(AD_seurat, reduction = "umap", pt.size = 0.2, cols = c( "grey",  "darkred"),
        group.by = "has_TCR") + ggtitle("")
DimPlot(AD_seurat, reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
        group.by = "is_clonal", ) + ggtitle("")
dev.off()

png(file.path(output,paste0("UMAP_has_TCR.png")), width = 1400, height = 1200, res = 250)
DimPlot(AD_seurat, reduction = "umap", pt.size = 0.1, cols = c( "grey",  "darkred"),
        group.by = "has_TCR", order = T) + ggtitle("")
dev.off()

png(file.path(output,paste0("UMAP_is_clonal.png")), width = 1400, height = 1200, res = 250)
DimPlot(AD_seurat, reduction = "umap", pt.size = 0.1, cols = c(  "grey", "darkred"),
        group.by = "is_clonal", order = T) + ggtitle("")
dev.off()


pdf(file.path(output, paste0("UMAP_per_sample.pdf")), width = 9, height = 8)

for(samp in unique(AD_seurat$sample_id)){
  AD_seurat. = AD_seurat
  AD_seurat.$is_here = FALSE
  AD_seurat.$is_here[AD_seurat$sample_id == samp] = TRUE
  tplot = DimPlot(AD_seurat., reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
                  group.by = "is_here", order = TRUE) + ggtitle(samp)
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
}
dev.off()

pdf(file.path(output, paste0("UMAP_per_clust.pdf")), width = 6, height = 5)
for(clust in unique(AD_seurat$seurat_clusters)){
  AD_seurat. = AD_seurat
  AD_seurat.$is_here = FALSE
  AD_seurat.$is_here[AD_seurat$seurat_clusters == clust] = TRUE
  tplot = DimPlot(AD_seurat., reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
                  group.by = "is_here", order = TRUE) + ggtitle(clust)
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
}
dev.off()

# Save the Seurat objects
qs::qsave(AD_seurat, file.path(output, "AD_seurat.qs"))

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
  if(i %in% rownames(AD_seurat)) print( FeaturePlot(AD_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()

pdf(file.path(output, paste0("Markers_expression_NK.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(NK_signature)){
  print(i)
  if(i %in% rownames(AD_seurat)) print( FeaturePlot(AD_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()


png(file.path(output,paste0("UMAP_sample_All.png")), width = 1400, height = 1200, res = 250)
tplot = DimPlot(AD_seurat, reduction = "umap", pt.size = 1.5, 
                cols = unique(AD_seurat$sample_id_color[order(AD_seurat$sample_id)]),
                group.by = "sample_id") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

# Annotation
AD_seurat. = ScaleData(AD_seurat, features = cell_markers)

png(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.png")), width = 3200, height = 1600, res = 200)
Seurat::DoHeatmap(AD_seurat., features = unique(cell_markers), cells = sample(ncol(AD_seurat.), 50000),group.by = "seurat_clusters")
dev.off()


pdf(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.pdf")), width = 16, height = 8)
Seurat::DoHeatmap(AD_seurat., features =  unique(cell_markers), cells = 1:1000, group.by = "seurat_clusters")
dev.off()

png(file.path(output, paste0("DotPlot_cell_markers_seurat_clusters.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(AD_seurat, features = unique(cell_markers), group.by = "seurat_clusters", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

png(file.path(output, paste0("DotPlot_gold_markers_seurat_clusters.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(AD_seurat, features = unique(gold_markers), group.by = "seurat_clusters", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()


# Ucell signature NK
AD_seurat = UCell::AddModuleScore_UCell(AD_seurat, list("NK" = NK_signature))

Idents(AD_seurat) = AD_seurat$celltype
png(file.path(output,paste0("UMAP_NK_signature.png")), width = 1400, height = 1400, res = 150)
FeaturePlot(AD_seurat, features = "NK_UCell", label = F, reduction = "umap", cols = c("grey", "darkred"))
dev.off()

# Absolute number of NK cells per sample
library(dittoSeq)
AD_seurat$condition = gsub("[0-9]*","",gsub("_.*", "",AD_seurat$sample_id))

png(file.path(output, "Cluster x Clonality.png"), width = 3000, height = 1600, res = 300)
d = dittoBarPlot(AD_seurat, "is_clonal", group.by = "seurat_clusters",
                 main = "Condition x Cell Type", color.panel = c("grey", "darkred"))  
d
dev.off()

png(file.path(output, "Sample x Cluster.png"), width = 3000, height = 1600, res = 300)
d = dittoBarPlot(AD_seurat, "sample_id", group.by = "seurat_clusters",
                 main = "Sample x Cluster", color.panel = unique(AD_seurat$sample_id_color[order(AD_seurat$sample_id)]))  
d
dev.off()

res = IDclust::differential_edgeR_pseudobulk_LRT(AD_seurat,
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
  "Epidermis_" = c(17, 5, 1, 6, 12, 16, 7, 25, 0, 7, 8, 27, 15, 31, 32, 30, 29, 19) ,
  "Follicular epidermis_" = c(),
  "B cells_" = c(),
  "T helper cells_" = c(14),
  "Natural Killer cells_" = c(21),
  "Cytotoxic T cells_" = c(3),
  "Macrophages and monocytes_" = c(10, 13),
  "Dendritic cells_" = c(22, 26),
  "Fibroblasts_" = c(2, 23, 11),
  "Melanocytes and Schwann cells_" = c(28),
  "Endothelial cells_" = c(4, 18),
  "Smooth muscle cells_" = c(9),
  "Clonal T cells_" = c(),
  "Platelets_"  = c(),
  "Gamma Delta T cells_"  = c(24),
  "Lymphatic Endothelial cells_"  = c(20)
))


AD_seurat$celltype = gsub("_.*", "", names(map_cluster_celltype))[match(AD_seurat$seurat_clusters, map_cluster_celltype)]

col_to_celltype = read.csv("output/col_to_celltype.csv")
AD_seurat$celltype_color = col_to_celltype$x[match(AD_seurat$celltype, col_to_celltype$X)]

png(file.path(output,paste0("UMAP_cell_type_All.png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(AD_seurat, reduction = "umap", pt.size = 0.75,
                cols = c(unique(AD_seurat$celltype_color[order(AD_seurat$celltype)])),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

for(gene in NK_signature){
  png(file.path(output,paste0("Violin_cell_type_",gene,".png")), width = 1800, height = 1200, res = 250)
  print(VlnPlot(AD_seurat, features = gene, group.by = "celltype",
                cols =c(unique(AD_seurat$celltype_color[order(AD_seurat$celltype)])) ))
  dev.off()
}


AD_seurat. = AD_seurat
AD_seurat.$celltype = factor(AD_seurat.$celltype, levels = unique(gsub("_.*", "", names(map_cluster_celltype))))
png(file.path(output, paste0("DotPlot_gold_markers_celltype.png")), width = 2400, height = 1800, res = 250)
Seurat::DotPlot(AD_seurat., features = gold_markers, group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

# Differential Analysis
library(IDclust)
list_res = list()

res = IDclust::differential_edgeR_pseudobulk_LRT(AD_seurat,
                                                 by = "celltype",
                                                 assay = "RNA",
                                                 biological_replicate_col = "sample_id",
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
  print(VlnPlot(AD_seurat, features = gene, group.by = "celltype",
                cols =c(unique(AD_seurat$celltype_color[order(AD_seurat$celltype)])) ))
  dev.off()
  
  png(file.path(output,paste0("UMAP_",gene,".png")), width = 1800, height = 1200, res = 250)
  print( FeaturePlot(AD_seurat,  cols = c("grey", "darkred"), features = gene)) 
  dev.off()
}

png(file.path(output, paste0("DotPlot_HLAs_celltype.png")), width = 2400, height = 1800, res = 250)
Seurat::DotPlot(AD_seurat., features = c("HLA-A", "HLA-B", "HLA-C", "HLA-E"), group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()


AD_seurat. = AD_seurat[, AD_seurat$celltype %in% c("Cytotoxic T cells", 
                                                   "Dendritic cells",
                                                   "Gamma Delta T cells",
                                                   "Macrophages and monocytes",
                                                   "Natural Killer cells",
                                                   "T helper cells",
                                                   "Gamma Delta T cells"
                                                   )]
png(file.path(output, "CellType x Sample x Condition Immune.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(AD_seurat., "celltype", group.by = "sample_id",
             color.panel = unique(AD_seurat.$celltype_color[order(AD_seurat.$celltype)]),
             main = "CellType x Sample x Condition")  + xlab("")
dev.off()


AD_seurat. = AD_seurat
AD_seurat.$is_NK = ifelse(AD_seurat.$celltype == "Natural Killer cells", "Yes", "No")
d =  dittoBarPlot(AD_seurat., "is_NK", group.by = "sample_id",
                  color.panel = sample_color$sample_color[match(levels(d$data$label),
                                                                sample_color$sample)],
                  main = "Natural Killer Cells",  scale = "count")  + xlab("")

png(file.path(output, "CellType x NK number.png"), width = 3000, height = 1600, res = 300)
df = d$data
df %>% filter(label == "Yes") %>%
  ggplot(aes(x = grouping, y = count, fill = grouping)) + geom_bar(stat = "identity") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = unique(
    AD_seurat.$sample_id_color)) + ggtitle("Natural Killer cells") + ylab("% Natural Killer cells") + xlab("")
dev.off()

png(file.path(output, "CellType x NK percent.png"), width = 3000, height = 1600, res = 300)
df = d$data %>% mutate(percent = 100 * percent) 
df %>% filter(label == "Yes") %>%
  ggplot(aes(x = grouping, y = percent, fill = grouping)) + geom_bar(stat = "identity") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = unique(
    AD_seurat.$sample_id_color)) + ggtitle("Natural Killer cells") + ylab("% Natural Killer cells") + xlab("")
dev.off()

d = dittoBarPlot(AD_seurat, "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  
png(file.path(output, "Sample x CellType x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(AD_seurat, "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

AD_seurat. = AD_seurat[,AD_seurat$celltype %in% unique(gsub("_.*","",names(map_cluster_celltype)))[c(3,4,5,6,7,8,9)]]
d = dittoBarPlot(AD_seurat., "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  

png(file.path(output, "Sample x CellType x Condition Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(AD_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "percent",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

png(file.path(output, "Sample x CellType x Condition Immune Compartment Count.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(AD_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

AD_seurat. = AD_seurat[,grep(group, AD_seurat$study)]

png(file.path(output,paste0("UMAP_celltype_", group,".png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(AD_seurat., reduction = "umap", pt.size = 0.75,
                cols =  unique(AD_seurat.$celltype_color[order(AD_seurat.$celltype)]),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
print(tplot)
dev.off()


png(file.path(output, "Condition x Clonality x Cell Type Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(AD_seurat., "is_clonal", group.by = "celltype", color.panel = c("grey", "darkred"),
             main = "Condition x Clonality x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()


png(file.path(output, "Sample x Clonality x Condition Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(AD_seurat., "is_clonal", group.by = "sample_id", color.panel = c("grey", "darkred"),
             main = "Condition x Cell Type Immune", scale = "percent",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()




pdf(file.path(seurat_dir_combined, paste0("Cell_type_Markers_expression_violin.pdf")), height = 10,
    width = 15, pointsize = 12)
for(i in cell_markers){
  print( VlnPlot(AD_seurat,  features = i,  pt.size = 1.5, group.by = "seurat_clusters"))
}
dev.off()



# Plot for each object the NK signature
for( i in 1:length(seurat_list)){
  seurat_obj = seurat_list[[i]]
  AD_seurat. = AD_seurat[,AD_seurat$sample_id == seurat_obj$sample_id[1]]
  seurat_obj$NK_Ucell = AD_seurat.$NK_UCell[match(colnames(seurat_obj),  gsub(paste0(".*",seurat_obj$sample_id[1], "_"),"",colnames(AD_seurat.)))]
  
  png(file.path(output,paste0("UMAP_NK_signature_", seurat_obj$sample_id[1],".png")), width = 1400, height = 1400, res = 150)
  print(FeaturePlot(seurat_obj, features = "NK_Ucell", label = F, reduction = "umap", cols = c("grey", "darkred")) + ggtitle(seurat_obj$sample_id[1]))
  dev.off()
  
}

# For AD study compare Tumor T cells vs Healthy T-cells

# For AD study compare Tumor T cells vs non-tumor T cells

# For AD study compare non-tumor T cells vs Healthy T-cells

# Save the Seurat objects
qs::qsave(AD_seurat, file.path(output, "AD_seurat.qs"))




# 
