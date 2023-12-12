setwd("/mnt/RECHERCHE/GUENOVA_LAB/Pacome/Results/scRNA_YunTsan/")
library(Seurat)

Seu = qs::qread("../scRNA_YunTsan/output/Seurat/Seu.qs")

table(Seu$tissue, Seu$TCR)

Seu = Seu[,which(!is.na(Seu$TCR))]
Seu$tumor = "NonTumor"
Seu$tumor[Seu$TCR %in% c("MainClone", "RelateMainClone")] = "Tumor"

res = FindMarkers(Seu, ident.1 = "Tumor", ident.2 = "NonTumor", group.by = "tumor")

library(dplyr)
View(res %>% arrange(desc(avg_log2FC)))


table(Seu$tumor)

AI_genes = c("RGS1", "SERF2", "EIF4A2", "SQSTM1", "BEST1", "RNF213",
             "IFITM2", "IFITM1", "HLA-A", "PGK1", "TOMM20", "RBM5",
             "PLIN2", "FYCO1", "HLA-E", "ZNF500", "S100A11", "ANXA2",
             "ACTG1", "ENSG00000262633", "KLF6", "ESYT2", "DNAJB6", 
             "HNRNPR", "MT-RNR1", "CCT6A", "SAT1", "MRNIP", "C4orf51",
             "MT-RNR2", "VDAC3", "TSPAN13", "CD96", "EEF1A1", "SLC2A3",
             "FTH1", "MT-CO1", "TRBC2", "DDX5", "UBC", "B2M", "TMSB4X", 
             "LGALS8", "HNRNPH1", "GLG1", "NCOR1", "MYL6", "ABCC6", "HSPA8", "LCP1")

res. = res %>% filter(p_val_adj < 0.01)

res["KLF6",]

Seu["PKG1",]

length(intersect(rownames(res.), AI_genes))  / length(AI_genes)
length(intersect(genes_fig, AI_genes))  / length(AI_genes)

WriteXLS::WriteXLS(res, ExcelFileName = "../scRNA_YunTsan/output/Tumor_vs_Non_Tumor_TCR_cells.xlsx")
