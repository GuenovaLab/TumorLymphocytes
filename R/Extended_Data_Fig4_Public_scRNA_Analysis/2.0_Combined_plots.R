#################################################################################
# Analysis of cancerous T lymphocytes 
########################################

library(Seurat)
library(scater)
library(dplyr)
library(dittoSeq)
library(cowplot)


maindir = "/home/localadmin/Documents/scRNA_YunTsan_Revision///"
setwd(maindir)
output = file.path(maindir, "output");  if(!dir.exists(output)){dir.create(output)}

MF_seurat = qs::qread("output/MF_seurat.qs")
AD_seurat = qs::qread("output/AD_seurat.qs")
BCC_seurat = qs::qread("output/BCC_seurat.qs")
CBCL_seurat = qs::qread("output/CBCL_seurat.qs")

table(MF_seurat$condition)

# UMAP Cell Type
plot0 <-  DimPlot(MF_seurat[,MF_seurat$condition == "HC"], reduction = "umap", pt.size = 0.75, 
                  cols = c(unique(MF_seurat$celltype_color[order(MF_seurat$celltype)])),
                  group.by = "celltype") + ggtitle("HC")
plot1 <-  DimPlot(MF_seurat[,MF_seurat$condition == "MF"], reduction = "umap", pt.size = 0.75, 
                  cols = c(unique(MF_seurat$celltype_color[order(MF_seurat$celltype)])),
                  group.by = "celltype") + ggtitle("MF")
plot2 <-  DimPlot(AD_seurat, reduction = "umap", pt.size = 0.75, 
                  cols = c(unique(AD_seurat$celltype_color[order(AD_seurat$celltype)])),
                  group.by = "celltype") + ggtitle("AD")
plot3 <-  DimPlot(BCC_seurat, reduction = "umap", pt.size = 0.75, 
                  cols = c(unique(BCC_seurat$celltype_color[order(BCC_seurat$celltype)])),
                  group.by = "celltype") + ggtitle("BCC")
plot4 <-  DimPlot(CBCL_seurat, reduction = "umap", pt.size = 0.75, 
                  cols = c(unique(CBCL_seurat$celltype_color[order(CBCL_seurat$celltype)])),
                  group.by = "celltype") + ggtitle("CBCL")

pdf(file.path(output, "UMAP_celltype_combined.pdf"), width = 15, height = 15)
combined_plot <- plot_grid(plot0, plot1, plot2, plot3, plot4, nrow = 3, ncol = 2)
print(combined_plot)
dev.off()

pdf(file.path(output, "UMAP_celltype_combined_wolegend.pdf"), width = 10, height = 15)
combined_plot <- plot_grid(plot0 + NoLegend(), plot1 + NoLegend(), plot2 + NoLegend(), plot3 + NoLegend(), plot4 + NoLegend(), nrow = 3, ncol = 2)
print(combined_plot)
dev.off()


# UMAP NK signature
plot0 <-  FeaturePlot(MF_seurat[,MF_seurat$condition == "HC"], reduction = "umap", pt.size = 0.75, 
                      cols = rev(viridis::magma(100)),
                      feature = "NK_UCell", min.cutoff = 0.1, max.cutoff = 0.9,
                      order =T) + ggtitle("HC")
plot1 <-  FeaturePlot(MF_seurat[,MF_seurat$condition == "MF"], reduction = "umap", pt.size = 0.75, 
                  cols = rev(viridis::magma(100)),
                  feature = "NK_UCell", min.cutoff = 0.1, max.cutoff = 0.9,
                  order =T) + ggtitle("MF")
plot2 <-  FeaturePlot(AD_seurat, reduction = "umap", pt.size = 0.75, 
                  cols = rev(viridis::magma(100)),
                  feature = "NK_UCell", min.cutoff = 0.1, max.cutoff = 0.9,
                  order =T) + ggtitle("AD")
plot3 <-  FeaturePlot(BCC_seurat, reduction = "umap", pt.size = 0.75, 
                  cols = rev(viridis::magma(100)),
                  feature = "NK_UCell", min.cutoff = 0.1, max.cutoff = 0.9,
                  order =T) + ggtitle("BCC")
plot4 <-  FeaturePlot(CBCL_seurat, reduction = "umap", pt.size = 0.75, 
                  cols = rev(viridis::magma(100)), 
                  feature = "NK_UCell", min.cutoff = 0.1, max.cutoff = 0.9,
                  order =T) + ggtitle("CBCL")

pdf(file.path(output, "UMAP_NK_signature_combined.pdf"), width = 3400, height = 15)
combined_plot <- plot_grid(plot0, plot1, plot2, plot3, plot4, nrow = 3, ncol = 2)
print(combined_plot)
dev.off()

pdf(file.path(output, "UMAP_NK_signature_combined_wolegend.pdf"), width = 10, height = 15)
combined_plot <- plot_grid(plot0 + NoLegend(), plot1 + NoLegend(), plot2 + NoLegend(), plot3 + NoLegend(), plot4 + NoLegend(), nrow = 3, ncol = 2)
print(combined_plot)
dev.off()

# Violin NK signature
plot1 <-  VlnPlot(MF_seurat, group.by = "celltype", pt.size = 0.05, 
                      cols =c(unique(MF_seurat$celltype_color[order(MF_seurat$celltype)])),
                      feature = "NK_UCell") + ggtitle("MF") + xlab("")
plot2 <-  VlnPlot(AD_seurat,  group.by = "celltype", pt.size = 0.05, 
                      cols = c(unique(AD_seurat$celltype_color[order(AD_seurat$celltype)])),
                      feature = "NK_UCell") + ggtitle("AD")  + xlab("")
plot3 <-  VlnPlot(BCC_seurat, group.by = "celltype",  pt.size = 0.05, 
                      cols =c(unique(BCC_seurat$celltype_color[order(BCC_seurat$celltype)])),
                      feature = "NK_UCell") + ggtitle("BCC")  + xlab("")
plot4 <-  VlnPlot(CBCL_seurat, group.by = "celltype", pt.size = 0.05, 
                      cols = c(unique(CBCL_seurat$celltype_color[order(CBCL_seurat$celltype)])), 
                      feature = "NK_UCell") + ggtitle("CBCL") + xlab("") 

pdf(file.path(output, "Violin_NK_signature_combined.pdf"), width = 20, height = 6)
combined_plot <- plot_grid(plot1, plot2, plot3, plot4, nrow = 1, ncol = 4)
print(combined_plot)
dev.off()

pdf(file.path(output, "Violin_NK_signature_combined_wolegend.pdf"),width = 15, height = 6, res = 250)
combined_plot <- plot_grid(plot1 + NoLegend(), plot2 + NoLegend(), plot3 + NoLegend(), plot4 + NoLegend(), nrow = 1, ncol = 4)
print(combined_plot)
dev.off()

# Markers

MF_res = readxl::read_xlsx("output/MF/One_vs_All_celltypes.xlsx")
AD_res = readxl::read_xlsx("output/AD/One_vs_All_celltypes.xlsx")
BCC_res = readxl::read_xlsx("output/BCC/One_vs_All_celltypes.xlsx")
CBCL_res = readxl::read_xlsx("output/CBCL/One_vs_All_celltypes.xlsx")

intersect(
intersect(
  MF_res$gene[MF_res$cluster == "Tumour cells"],
  BCC_res$gene[BCC_res$cluster == "Tumour cells"]
),
CBCL_res$gene[CBCL_res$cluster == "Tumour cells"]
)

markers = c(
  "PTPRC",
  "CD3E",
  "CD3D",
  "CD4",
  "NKG7",
  "CD8A",
  "MS4A1",
  "CD79A",
  "CD14",
  "CD163",
  "LAMP3",
  "MS4A2",
  "KRT6A",
  "KRT6B",
  "KRT19",
  "CD34",
  "PECAM1",
  "ELF5",
  "THBS1",
  "FLT4",
  "COL1A1",
  "COL1A2",
  "ACTA2",
  "DCT",
  "MLANA",
  "TRBV11-2",
  "IGLV4-69",
  "EPCAM",
  "STMN1",
  "IDH2",
  "HMGN1",
  "SNRPE",
  "ILF2",
  "ECH1" 
)

celltypes = unlist(list(
  "T helper cells",
  "Gamma Delta T cells",
  "Natural Killer cells",
  "Cytotoxic T cells",
  "B cells",
  "Macrophages and monocytes",
  "Dendritic cells",
  "Mast cells",
  "Epidermis" ,
  "Follicular epidermis",
  "Endothelial cells",
  "Lymphatic Endothelial cells",
  "Platelets",
  "Erythroid cells",
  "Fibroblasts",
  "Smooth muscle cells",
  "Melanocytes and Schwann cells",
  "Tumour cells",
  "Clonal T cells"
))

celltypes = forcats::as_factor(celltypes)
MF_seurat$celltype = factor(MF_seurat$celltype, levels = levels(celltypes))
AD_seurat$celltype = factor(AD_seurat$celltype, levels = levels(celltypes))
BCC_seurat$celltype = factor(BCC_seurat$celltype, levels = levels(celltypes))
CBCL_seurat$celltype = factor(CBCL_seurat$celltype, levels = levels(celltypes))

# DotPlot NK signature
plot1 <-  DotPlot(MF_seurat, group.by = "celltype", 
                  cols = c("grey", "gold"),
                  features = markers) + ggtitle("MF") + xlab("") + ylab("") + theme(axis.text.x =  element_text(angle = 90))
plot2 <-  DotPlot(AD_seurat,  group.by = "celltype",
                  cols = c("grey", "gold"),
                  features = markers) + ggtitle("AD")  + xlab("") + ylab("") + theme(axis.text.x =  element_text(angle = 90))
plot3 <-  DotPlot(BCC_seurat, group.by = "celltype",
                  cols = c("grey", "gold"),
                  features = markers) + ggtitle("BCC")  + xlab("") + ylab("") + theme(axis.text.x =  element_text(angle = 90))
plot4 <-  DotPlot(CBCL_seurat, group.by = "celltype",
                  cols = c("grey", "gold"),
                  features = markers) + ggtitle("CBCL") + xlab("") + ylab("") + theme(axis.text.x =  element_text(angle = 90)) 

pdf(file.path(output, "DotPlot_celltypes_combined.pdf"), width = 15, height = 3000)
combined_plot <- plot_grid(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
print(combined_plot)
dev.off()

pdf(file.path(output, "DotPlot_celltypes_combined_no_legend.pdf"),width = 15, height = 3000)
combined_plot <- plot_grid(plot1 + NoLegend(), plot2 + NoLegend(), plot3 + NoLegend(), plot4 + NoLegend(), nrow = 2, ncol = 2)
print(combined_plot)
dev.off()


# Composition Analysis
MF_seurat$celltype = as.character(MF_seurat$celltype)
AD_seurat$celltype = as.character(AD_seurat$celltype)
BCC_seurat$celltype = as.character(BCC_seurat$celltype)
CBCL_seurat$celltype = as.character(CBCL_seurat$celltype)

plot1 <- dittoBarPlot(MF_seurat, "celltype", group.by = "sample_id",
                      color.panel = unique(MF_seurat$celltype_color[order(MF_seurat$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("MF")
plot2 <- dittoBarPlot(AD_seurat, "celltype", group.by = "sample_id",
                      color.panel = unique(AD_seurat$celltype_color[order(AD_seurat$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("AD")
plot3 <- dittoBarPlot(BCC_seurat, "celltype", group.by = "sample_id",
                      color.panel = unique(BCC_seurat$celltype_color[order(BCC_seurat$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("BCC")
plot4 <- dittoBarPlot(CBCL_seurat, "celltype", group.by = "sample_id",
                      color.panel = unique(CBCL_seurat$celltype_color[order(CBCL_seurat$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("CBCL")

rel_widths <- c(15,6,4,2)

pdf(file.path(output, "Composition_celltypes_combined.pdf"), width = 20, height = 6)
combined_plot <- plot_grid(plot1, plot2, plot3, plot4, nrow = 1, ncol = 4, rel_widths = rel_widths)
print(combined_plot)
dev.off()

pdf(file.path(output, "Composition_celltypes_combined_no_legend.pdf"),width = 15, height = 6)
combined_plot <- plot_grid(plot1 + NoLegend(), plot2 + NoLegend(), plot3 + NoLegend(), plot4 + NoLegend(),
                           nrow = 1, ncol = 4, rel_widths = rel_widths)
print(combined_plot)
dev.off()

# Immune
immune_celltypes = c("Cytotoxic T cells", 
                     "Dendritic cells",
                     "Gamma Delta T cells",
                     "Macrophages and monocytes",
                     "Natural Killer cells",
                     "T helper cells",
                     "B cells",
                     "Gamma Delta T cells",
                     "Mast cells"
)

MF_seurat_healthy = MF_seurat[, (MF_seurat$celltype %in% immune_celltypes) & (MF_seurat$condition == "HC")]
MF_seurat_healthy$celltype[MF_seurat_healthy$celltype == "Tumour cells"] = "T helper cells"
MF_seurat_non_healthy = MF_seurat[, (MF_seurat$celltype %in% immune_celltypes) & (MF_seurat$condition != "HC")]
AD_seurat. = AD_seurat[, AD_seurat$celltype %in% immune_celltypes]
BCC_seurat. = BCC_seurat[, BCC_seurat$celltype %in% immune_celltypes]
CBCL_seurat. = CBCL_seurat[, CBCL_seurat$celltype %in% immune_celltypes]

plot1 <- dittoBarPlot(MF_seurat_healthy, "celltype", group.by = "sample_id",
                      color.panel = unique(MF_seurat_healthy$celltype_color[order(MF_seurat_healthy$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("Healthy")
plot1bis <- dittoBarPlot(MF_seurat_non_healthy, "celltype", group.by = "sample_id",
                      color.panel = unique(MF_seurat_non_healthy$celltype_color[order(MF_seurat_non_healthy$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("MF")

plot2 <- dittoBarPlot(AD_seurat., "celltype", group.by = "sample_id",
                      color.panel = unique(AD_seurat.$celltype_color[order(AD_seurat.$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("AD")
plot3 <- dittoBarPlot(BCC_seurat., "celltype", group.by = "sample_id",
                      color.panel = unique(BCC_seurat.$celltype_color[order(BCC_seurat.$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("BCC")
plot4 <- dittoBarPlot(CBCL_seurat., "celltype", group.by = "sample_id",
                      color.panel = unique(CBCL_seurat.$celltype_color[order(CBCL_seurat.$celltype)]),
                      main = "CellType x Sample x Condition")  + xlab("") + ylab("") + ggtitle("CBCL")

rel_widths <- c(4,11,6,4,2)

pdf(file.path(output, "Composition_immune_celltypes_combined.pdf"), width = 20, height = 6)
combined_plot <- plot_grid(plot1, plot1bis, plot2, plot3, plot4, nrow = 1, ncol = 4, rel_widths = rel_widths)
print(combined_plot)
dev.off()

pdf(file.path(output, "Composition_immune_celltypes_combined_no_legend.pdf"), width = 20, height = 7)
combined_plot <- plot_grid(plot1 + NoLegend(), plot1bis + NoLegend(), plot2 + NoLegend(), plot3 + NoLegend(), plot4 + NoLegend(),
                           nrow = 1, ncol = 5, rel_widths = rel_widths)
print(combined_plot)
dev.off()


# Real number
plot1 <- dittoBarPlot(MF_seurat_healthy, "celltype", group.by = "sample_id",
                      color.panel = unique(MF_seurat_healthy$celltype_color[order(MF_seurat_healthy$celltype)]),
                      main = "CellType x Sample x Condition", scale = "count", min = 0, max = 8500) + xlab("") + ylab("") + ggtitle("Healthy")
plot1bis <- dittoBarPlot(MF_seurat_non_healthy, "celltype", group.by = "sample_id",
                         color.panel = unique(MF_seurat_non_healthy$celltype_color[order(MF_seurat_non_healthy$celltype)]),
                         main = "CellType x Sample x Condition", scale = "count", min = 0, max = 8500)  + xlab("") + ylab("") + ggtitle("MF")

plot2 <- dittoBarPlot(AD_seurat., "celltype", group.by = "sample_id",
                      color.panel = unique(AD_seurat.$celltype_color[order(AD_seurat.$celltype)]),
                      main = "CellType x Sample x Condition", scale = "count", min = 0, max = 8500)  + xlab("") + ylab("") + ggtitle("AD")
plot3 <- dittoBarPlot(BCC_seurat., "celltype", group.by = "sample_id",
                      color.panel = unique(BCC_seurat.$celltype_color[order(BCC_seurat.$celltype)]),
                      main = "CellType x Sample x Condition", scale = "count", min = 0, max = 8500)  + xlab("") + ylab("") + ggtitle("BCC")
plot4 <- dittoBarPlot(CBCL_seurat., "celltype", group.by = "sample_id",
                      color.panel = unique(CBCL_seurat.$celltype_color[order(CBCL_seurat.$celltype)]),
                      main = "CellType x Sample x Condition", scale = "count", min = 0, max = 8500)  + xlab("") + ylab("") + ggtitle("CBCL")

rel_widths <- c(4,11,6,4,2)

pdf(file.path(output, "Composition_immune_celltypes_number_combined.pdf"), width = 20, height = 7)
combined_plot <- plot_grid(plot1, plot1bis, plot2, plot3, plot4, nrow = 1, ncol = 5, rel_widths = rel_widths)
print(combined_plot)
dev.off()

pdf(file.path(output, "Composition_immune_celltypes_number_combined_no_legend.pdf"), width = 20, height = 7)
combined_plot <- plot_grid(plot1 + NoLegend(), plot1bis + NoLegend(), plot2 + NoLegend(), plot3 + NoLegend(), plot4 + NoLegend(),
                           nrow = 1, ncol = 5, rel_widths = rel_widths)
print(combined_plot)
dev.off()


