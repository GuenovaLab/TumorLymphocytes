#################################################################################
# Analysis of cancerous T lymphocytes 
########################################

library(Seurat)
library(scater)

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

# Directories -------------------------------------------------------------
maindir = "/mnt/RECHERCHE/GUENOVA_LAB/Pacome/Results/scRNA_YunTsan_Revision/"
output = file.path(maindir, "output");  if(!dir.exists(output)){dir.create(output)}
QCdir = file.path(output, "QC");  if(!dir.exists(QCdir)){dir.create(QCdir)}

# Replace GEO path with the path to the GEO directory downloaded on your computer (scRNA MM468)
GEOdir_scRNA = file.path("/mnt/RECHERCHE/GUENOVA_LAB/Project_Yun-Tsan_ADCC_paper/dataset to analyze for ADCC paper/") 

#################################################################################
# Loading & QC
####################

# Load required libraries
library(Seurat)
library(harmony)

# Get a list of all subdirectories (studies)
studies <- list.files(GEOdir_scRNA, full.names = TRUE)

# Initialize a list to store Seurat objects
seurat_list <- list()

library(DoubletFinder)

# Loop through each study
pdf(file.path(output, "scRNA_sample_per_sample.pdf"))
for (study_dir in studies) {
  # Get the study name
  study_name <- basename(study_dir)
  print(study_name)
  
  samples = unique(gsub("_clonotypes.*|_filtered_.*|_barcodes.*|_features.*|_matrix.*", "", list.files(study_dir)))
  for(sample in samples){
    
    print(sample)
    # Read the matrix files
    barcode_file <- file.path(study_dir, paste0(sample, "_barcodes.tsv.gz"))
    features_file <- file.path(study_dir, paste0(sample, "_features.tsv.gz"))
    matrix_file <- file.path(study_dir, paste0(sample, "_matrix.mtx.gz"))
    
    # Read the TCR sequence file
    tcr_file <- file.path(study_dir, paste0(sample, "_filtered_contig_annotations.csv"))
    
    # Load the matrix
    matrix <- ReadMtx(mtx = matrix_file, cells = barcode_file, features = features_file)
    
    # Create a Seurat object
    seurat_obj <- CreateSeuratObject(counts = matrix)
    seurat_obj$study <- study_name
    seurat_obj$sample_id <- sample
    
    # Add TCR sequence using the 'barcode' column as ID
    if(file.exists(tcr_file)){
      tcr_data <- read.csv(tcr_file)
      seurat_obj$TCR <- tcr_data$contig_id[match(colnames(seurat_obj), tcr_data$barcode)]
      seurat_obj$TCR_chain <- tcr_data$chain[match(colnames(seurat_obj), tcr_data$barcode)]
      seurat_obj$TCR_v_gene <- tcr_data$v_gene[match(colnames(seurat_obj), tcr_data$barcode)]
      seurat_obj$TCR_d_gene <- tcr_data$d_gene[match(colnames(seurat_obj), tcr_data$barcode)]
      seurat_obj$TCR_j_gene <- tcr_data$j_gene[match(colnames(seurat_obj), tcr_data$barcode)]
      seurat_obj$TCR_c_gene <- tcr_data$c_gene[match(colnames(seurat_obj), tcr_data$barcode)]
      seurat_obj$TCR_raw_clonotype_id <- tcr_data$raw_clonotype_id[match(colnames(seurat_obj), tcr_data$barcode)]
      seurat_obj$TCR_raw_consensus_id <- tcr_data$raw_consensus_id[match(colnames(seurat_obj), tcr_data$barcode)]
    }
  
  
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")

    qc_thresholds <- c(min_genes = 400, max_genes = 8000, min_cells = 3, max_mito = 10)
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > qc_thresholds["min_genes"] &
                           nFeature_RNA < qc_thresholds["max_genes"] &
                           nCount_RNA > qc_thresholds["min_cells"] &
                           percent.mt < qc_thresholds["max_mito"])
    
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, npcs = 30)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
    
    print(DimPlot(seurat_obj) + ggtitle(sample))
    
    gc()
    sweep.res <- paramSweep_v3(seurat_obj, PCs = 1:30, sct = FALSE, num.cores = 10)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    
    bcmvn <- find.pK(sweep.stats)
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    homotypic.prop <- modelHomotypic(seurat_obj$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.075* nrow(seurat_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:10, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$MeanBC)])),
                                   nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE, )
    
    seurat_obj$doublet = seurat_obj@meta.data[,ncol(seurat_obj@meta.data)]
    
    print(DimPlot(seurat_obj, group.by = "doublet") + ggtitle(sample)) 
    
    seurat_obj = seurat_obj[, which(seurat_obj$doublet != "Doublet")]
    
    # Store the Seurat object in the list
    seurat_list[[sample]] <- seurat_obj
    gc()
}
}
dev.off()

qs::qsave(seurat_list, file.path(output, "seurat_list.qs"))


# Merge all Seurat objects
merged_seurat = Reduce(function(x,y) merge(x,y, add.cell.ids = c(x$sample_id[1],y$sample_id[1])) , seurat_list)
merged_seurat = merged_seurat[,-which(merged_seurat$sample_id == "MF303_Skin")]
gc()

# Run standard QC

# Run standard Seurat processing
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat, npcs = 30)
merged_seurat <- RunHarmony(merged_seurat, group.by.vars = "study")
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.6)


# Add metadata
merged_seurat$sample_id_color = sample_color$sample_color[match(merged_seurat$sample_id, sample_color$sample)]
merged_seurat$has_TCR = ifelse(is.na(merged_seurat$TCR_raw_consensus_id) | merged_seurat$TCR_raw_consensus_id == "None", FALSE, TRUE)
clones = table(merged_seurat$TCR_raw_consensus_id)
clones = clones[clones>100]
merged_seurat$is_clonal <- ifelse(merged_seurat$TCR_raw_consensus_id %in% names(clones), TRUE, FALSE)


# Perform additional downstream analysis as needed
pdf(file.path(output, paste0("scRNA.pdf")))
print(DimPlot(merged_seurat, reduction = "umap", group.by = "seurat_clusters"))
print(DimPlot(merged_seurat, reduction = "umap", group.by = "sample_id"))
print(DimPlot(merged_seurat, reduction = "umap", group.by = "study"))
print(FeaturePlot(merged_seurat,  features = "nCount_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
print(FeaturePlot(merged_seurat,  features = "nFeature_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
DimPlot(merged_seurat, reduction = "umap", pt.size = 0.2, cols = c( "grey",  "darkred"),
        group.by = "has_TCR") + ggtitle("")
DimPlot(merged_seurat, reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
        group.by = "is_clonal", ) + ggtitle("")
dev.off()

# Filter clusters with less than 200 cells (< 0.1%)
occ = table(merged_seurat$seurat_clusters)
merged_seurat = merged_seurat[,-which(merged_seurat$seurat_clusters %in% as.numeric(names(occ)[occ < 200]))]

png(file.path(output,paste0("UMAP_has_TCR.png")), width = 1400, height = 1200, res = 250)
DimPlot(merged_seurat, reduction = "umap", pt.size = 0.2, cols = c( "grey",  "darkred"),
                group.by = "has_TCR") + ggtitle("")
dev.off()

png(file.path(output,paste0("UMAP_is_clonal.png")), width = 1400, height = 1200, res = 250)
DimPlot(merged_seurat, reduction = "umap", pt.size = 0.2, cols = c(  "grey", "darkred"),
        group.by = "is_clonal", ) + ggtitle("")
dev.off()


# Save the Seurat objects
qs::qsave(merged_seurat, file.path(output, "merged_seurat.qs"))

cell_markers = c("PTPRC", "CD3E","CD4","CD8A", "FOXP3", "CD8B","TIGIT","IL2RA", 
                 "LAMP3","CD40",  "CD14", "CD16", "ITGAM", "HLA-DRA","HLA-DRQA1","HLA-DRQB1", "HLA-DRB1", "CD163", "CD1A", "CD207",
                 "MS4A2", "MS4A1",
                 "CDH5","PECAM1","CD34","KRT5","KRT10","KRT14","KRT1", "FABP5", "SPPR2G", "CDSN", 
                 "ITGA6", "ITGB1", "GRHL1", 
                 "CD200", "SOX9", "KRT19", "APOC1", "ACSL5", "ABCC3",
                 "LUM", "PDGFRB", "COL1A1",
                 "SOX10", "S100B", "DCT", "TYRP1", "MLANA")

NK_signature = c("NKG7",
                 "NKG2",
                 "KLRB1",
                 "KLRC1",
                 "KLRD1",
                 "KLRK1",
                 "CD7",
                 "TRA",
                 "GZMB",
                 "GNLY",
                 "PRF1",
                 "IL2RB",
                 "TBX21",
                 "NCAM1",
                 "CXCR4",
                 "TRAC",
                 "STYK1",
                 "TRBC1")

pdf(file.path(output, paste0("Markers_expression.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(cell_markers)){
  print(i)
  if(i %in% rownames(merged_seurat)) print( FeaturePlot(merged_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()

pdf(file.path(output, paste0("Markers_expression_NK.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(NK_signature)){
  print(i)
  if(i %in% rownames(merged_seurat)) print( FeaturePlot(merged_seurat,  cols = c("grey", "darkred"), features = i)) 
}
dev.off()


png(file.path(output,paste0("UMAP_sample_All.png")), width = 1400, height = 1200, res = 250)
tplot = DimPlot(merged_seurat, reduction = "umap", pt.size = 1.5, 
                cols = unique(merged_seurat$sample_id_color[order(merged_seurat$sample_id)]),
                group.by = "sample_id") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()


for(group in unique(merged_seurat$study)){
  merged_seurat. = merged_seurat[,grep(group, merged_seurat$study)]
  png(file.path(output, paste0("UMAP_sample_", group,".png")), width = 1400, height = 1200, res = 250)
  tplot = DimPlot(merged_seurat., reduction = "umap", pt.size = 0.2, cols = unique(merged_seurat.$sample_id_color[order(merged_seurat.$sample_id)]),
                  group.by = "sample_id") + ggtitle("")
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
  dev.off()
}

# Annotation
merged_seurat. = ScaleData(merged_seurat, features = cell_markers)

png(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.png")), width = 3200, height = 1600, res = 200)
Seurat::DoHeatmap(merged_seurat., features = unique(cell_markers), cells = sample(ncol(merged_seurat.), 50000),group.by = "seurat_clusters")
dev.off()


pdf(file.path(output,paste0("Heatmap_cell_markers_seurat_clusters.pdf")), width = 16, height = 8)
Seurat::DoHeatmap(merged_seurat., features =  unique(cell_markers), cells = 1:1000, group.by = "seurat_clusters")
dev.off()

png(file.path(output, paste0("DotPlot_cell_markers_seurat_clusters.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(merged_seurat, features = unique(cell_markers), group.by = "seurat_clusters", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

map_cluster_celltype = unlist(list(
  "Epidermis_" = c(1, 6, 7, 10,14, 24, 25) ,
  "Follicular epidermis_" = c(28),
  "B cells_" = c(15),
  "T helper cells_" = c(3),
  "Cytotoxic T cells_" = c(13),
  "Macrophages and monocytes_" = c(11),
  "Myeloid cells_" = c(5),
  "Dendritic cells_" = c(21),
  "Fibroblasts_" = c(0, 12, 16, 27),
  "Sebaceous glands_" = c(),
  "Melanocytes and Schwann cells_" = c(19),
  "Endothelial cells_" = c(4, 18, 26),
  "Smooth muscle cells_" = c(9),
  "Tumour cells_" = c(2,8,17,20,22,23)
))


merged_seurat$celltype = gsub("_.*", "", names(map_cluster_celltype))[match(merged_seurat$seurat_clusters, map_cluster_celltype)]

col_to_celltype = read.csv("output/col_to_celltype.csv")
merged_seurat$celltype_color = col_to_celltype$x[match(merged_seurat$celltype, col_to_celltype$X)]

png(file.path(output,paste0("UMAP_cell_type_All.png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(merged_seurat, reduction = "umap", pt.size = 0.75,
                cols = unique(merged_seurat$celltype_color[order(merged_seurat$celltype)]),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

for(group in unique(merged_seurat$study)){
  merged_seurat. = merged_seurat[,grep(group, merged_seurat$study)]
  png(file.path(output,paste0("UMAP_celltype_", group,".png")), width = 1800, height = 1200, res = 250)
  tplot = DimPlot(merged_seurat., reduction = "umap", pt.size = 0.75,
                  cols =  unique(merged_seurat.$celltype_color[order(merged_seurat.$celltype)]),
                  group.by = "celltype") + ggtitle("")
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
}

# Ucell signature NK
merged_seurat = UCell::AddModuleScore_UCell(merged_seurat, list("NK" = NK_signature))

Idents(merged_seurat) = merged_seurat$celltype
png(file.path(output,paste0("UMAP_NK_signature_", group,".png")), width = 1400, height = 1400, res = 150)
FeaturePlot(merged_seurat, features = "NK_UCell", label = T, reduction = "umap", cols = c("grey", "darkred"))
dev.off()





# Absolute number of NK cells per sample
library(dittoSeq)
merged_seurat$condition = gsub("[0-9]*","",gsub("_.*", "",merged_seurat$sample_id))

png(file.path(output, "Cluster x Clonality.png"), width = 3000, height = 1600, res = 300)
d = dittoBarPlot(merged_seurat, "is_clonal", group.by = "seurat_clusters",
                 main = "Condition x Cell Type", color.panel = c("grey", "darkred"))  
d
dev.off()

png(file.path(output, "Sample x Cluster.png"), width = 3000, height = 1600, res = 300)
d = dittoBarPlot(merged_seurat, "sample_id", group.by = "seurat_clusters",
                 main = "Sample x Cluster", color.panel = unique(merged_seurat$sample_id_color[order(merged_seurat$sample_id)]))  
d
dev.off()

d = dittoBarPlot(merged_seurat, "sample_id", group.by = "celltype",
                 color.panel =  unique(merged_seurat$sample_id_color[order(merged_seurat$sample_id)]),
                 main = "Condition x Cell Type")  
png(file.path(output, "CellType x Sample x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(merged_seurat, "sample_id", group.by = "celltype",
             color.panel = sample_color$sample_color[match(levels(d$data$label),
                                                           sample_color$sample)],
             main = "CellType x Sample x Condition")  + xlab("")
dev.off()



d = dittoBarPlot(merged_seurat, "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  
png(file.path(output, "Sample x CellType x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(merged_seurat, "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

merged_seurat. = merged_seurat[,merged_seurat$celltype %in% unique(gsub("_.*","",names(map_cluster_celltype)))[c(3,4,5,6,7,8,9)]]
d = dittoBarPlot(merged_seurat., "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  

png(file.path(output, "Sample x CellType x Condition Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(merged_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "percent",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

png(file.path(output, "Sample x CellType x Condition Immune Compartment Count.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(merged_seurat., "celltype", group.by = "sample_id", color.panel = col_to_celltype$x[match(levels(d$data$label), col_to_celltype$X)],
             main = "Condition x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

merged_seurat. = merged_seurat[,grep(group, merged_seurat$study)]

png(file.path(output,paste0("UMAP_celltype_", group,".png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(merged_seurat., reduction = "umap", pt.size = 0.75,
                cols =  unique(merged_seurat.$celltype_color[order(merged_seurat.$celltype)]),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
print(tplot)
dev.off()


png(file.path(output, "Condition x Clonality x Cell Type Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(merged_seurat., "is_clonal", group.by = "celltype", color.panel = c("grey", "darkred"),
             main = "Condition x Clonality x Cell Type Immune", scale = "count",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()


png(file.path(output, "Sample x Clonality x Condition Immune Compartment.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(merged_seurat., "is_clonal", group.by = "sample_id", color.panel = c("grey", "darkred"),
             main = "Condition x Cell Type Immune", scale = "percent",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()




pdf(file.path(seurat_dir_combined, paste0("Cell_type_Markers_expression_violin.pdf")), height = 10,
    width = 15, pointsize = 12)
for(i in cell_markers){
  print( VlnPlot(Seu_combined,  features = i,  pt.size = 1.5, group.by = "seurat_clusters"))
}
dev.off()



# Plot for each object the NK signature
for( i in 1:length(seurat_list)){
  seurat_obj = seurat_list[[i]]
  merged_seurat. = merged_seurat[,merged_seurat$sample_id == seurat_obj$sample_id[1]]
  seurat_obj$NK_Ucell = merged_seurat.$NK_UCell[match(colnames(seurat_obj),  gsub(paste0(".*",seurat_obj$sample_id[1], "_"),"",colnames(merged_seurat.)))]
  
  png(file.path(output,paste0("UMAP_NK_signature_", seurat_obj$sample_id[1],".png")), width = 1400, height = 1400, res = 150)
  print(FeaturePlot(seurat_obj, features = "NK_Ucell", label = F, reduction = "umap", cols = c("grey", "darkred")) + ggtitle(seurat_obj$sample_id[1]))
  dev.off()
  
}

# For MF study compare Tumor T cells vs Healthy T-cells

# For MF study compare Tumor T cells vs non-tumor T cells

# For MF study compare non-tumor T cells vs Healthy T-cells

# Save the Seurat objects
qs::qsave(merged_seurat, file.path(output, "merged_seurat.qs"))


# Clonal T cells vs Non Clonal T cells
seu_T = merged_seurat[,merged_seurat$celltype %in% c("T_cytotoxic","Natural Killer cells","T helper cells")]
seu_T = seu_T[,seu_T$is_clonal ]
library(ggpubr)
VlnPlot(seu_T, features = "IL32", group.by = "is_clonal", split.by = "condition")

gene = as.numeric(seu_T@assays$RNA@data["IL32",])
correlations <- apply(seu_T@assays$RNA@data, 1, function(x){cor(gene,as.numeric(x))})
hist(correlations, breaks = 150)
names(correlations)[which(correlations>0.5)]


# Boxplot stage
merged_seurat = qs::qread(file.path(output, "merged_seurat.qs"))
dim(merged_seurat)





