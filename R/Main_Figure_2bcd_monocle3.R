if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("cachem")
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')




library(monocle3)
library(ggplot2)
library(paletteer)



setwd("M:/DER/RECHERCHE/GUENOVA_LAB/Yun-Tsan Chang/CTCL_research/project_scRNA_seq/GeneVia/raw_data/ALL/four_definitions/skin/monocle3_test/")

matrix <- read.csv("skin_clone_raw_counts_TCR_without_Tc_matrix.csv", header = TRUE, row.names = 1)
cell <- read.csv("skin_clone_raw_counts_TCR_without_Tc_meta.csv", header=TRUE, row.names = 1)
gene <- read.csv("skin_raw_counts_TCR_without_Tc_GeneList.csv", header = TRUE, row.names = 1)

matrix <- as.matrix(matrix)



# Make the CDS object
cds <- new_cell_data_set(matrix,
                         cell_metadata = cell,
                         gene_metadata = gene)


# pre-process the data
cds <- preprocess_cds(cds, num_dim = 161)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "group")


# Reduve dimensionality and visualize the cells
cds <- reduce_dimension(cds)

plot_cells(cds, color_cells_by="StageNow", cell_size = 1, label_cell_groups = F)




cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

plot_cells(cds, reduction_method="UMAP",color_cells_by="StageNow", label_cell_groups=F,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5, cell_size = 1) 


cds <- order_cells(cds)
  

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1)



plot_cells(cds,
           color_cells_by = "StageNow",
           label_cell_groups=FALSE,
           label_leaves=T,
           label_branch_points=F,
           graph_label_size=1.5, cell_size = 1) +
  scale_color_manual (values = c("#0080FF", "#CCFF99", "#FF8000", "#FFCCFF", "#FF3333"))


AFD_genes <- c("BEST1", "RNF213", "EIF4A2", "HLA-B")


AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$definition %in% c("MainClone")]


plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=300, ncol = 2)







?plot_genes_in_pseudotime









cds = cluster_cells(cds, reduction_method = "UMAP")


plot_cells(cds,
           color_cells_by = "StageNow",
           label_groups_by_cluster=FALSE,
           label_leaves=T,
           label_branch_points=T, graph_label_size=1.5, cell_size = 1)


plot_cells(cds,
           color_cells_by = "StageNow",
           label_cell_groups=FALSE,
           label_leaves=T,
           label_branch_points=F,
           graph_label_size=1.5, cell_size = 1) +
  scale_color_manual (values = c("#0080FF", "#CCFF99", "#FF8000", "#FFCCFF", "#FF3333"))

cds <- order_cells(cds)


plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1)



plot_cells(cds, show_trajectory_graph=FALSE)




AFD_genes <- c("BEST1", "RNF213", "EIF4A2", "HLA-A", "IFITM1", "IFITM2", "RGS1", "SERF2", "SQSTM1", "TOMM20")

AFD_genes <- c("HLA-A", "HLA-B", "HLA-C", "IL32")


AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$definition %in% c("MainClone")]


plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=300, ncol = 1)



?plot_genes_in_pseudotime

View(cds)





cds_subset <- choose_cells(cds, reduction_method="UMAP")


cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)
plot_cells(cds,reduction_method = "PCA", cell_size = 1)




  















