library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(dplyr)
library(readxl)

# Defines a colorblind-friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Imports dataset called "combined" that has the columns "replicate," "variable," and "value"
mat_filtered = read.csv("../../../Project_Yun-Tsan_ADCC_paper/raw data and R code/Figure_6/6A_ridge plot/all_healthy_raw_counts_TCR_without_Tc.csv")
metadata = read.csv("../../../Project_Yun-Tsan_ADCC_paper/raw data and R code/Figure_6/6A_ridge plot/all_healthy_raw_counts_TCR_without_Tc_group.csv")

metadata$clone_info = gsub("_.*","",metadata$labels)
metadata = metadata %>% dplyr::filter(group != "Healthy" & compartment == "skin" & clone_info %in% c("MainClone","RelateMainClone"))
metadata$sample_id = gsub("_.*","",gsub(".*Clone_","",metadata$labels))
mat_filtered. = mat_filtered[,match(metadata$labels, colnames(mat_filtered))]
mat_filtered. = as.matrix(mat_filtered.)
rownames(mat_filtered.) = mat_filtered$Gene
rownames(metadata) = metadata$labels

Seu = CreateSeuratObject(mat_filtered., meta.data = metadata)
Seu = Seurat:::SCTransform(Seu)

meta = Seu@meta.data

meta$StageNow_Simplified = meta$StageNow
meta$StageNow_Simplified[meta$StageNow_Simplified %in% c("IB")] = "I"
meta$StageNow_Simplified[meta$StageNow_Simplified %in% c("IIB")] = "II"
meta$StageNow_Simplified[meta$StageNow_Simplified %in% c("IVA1", "IVA2", "IVB")] = "IV"


for(gene in c("HLA.A","HLA.B","HLA.C","BEST1", "EIF4A2", "IFITM1", "IFITM2", "RGS1", "RNF213", "SERF2", "SQSTM1", "TOMM20")){
  meta. = meta
  
  if(gene %in% c("HLA.A","HLA.B","HLA.C")){
    meta.$value = colSums(Seu@assays$SCT@data[grep(paste0("^", gene),rownames(Seu)),])
    
  } else{
    meta.$value = Seu@assays$SCT@data[gene,]
    
  }
  
  
  # Orders the variables on x-axis
  meta.$StageNow_Simplified <- factor(meta.$StageNow_Simplified, levels = c("I", "II", "IV"))
  rownames(meta.) = NULL
  # Calculates averages of each replicate
  ReplicateAverages <- meta. %>% select(sample_id, StageNow_Simplified, value) %>%
    group_by(StageNow_Simplified, sample_id) %>%
    summarise_each(list(mean))
  ReplicateAverages
  
  # Calculates total averages
  TotalAverages <- ReplicateAverages %>% summarise_each(list(mean))
  TotalAverages
  
  
  # Gives the p-value for the t-Test of variable 1 and 2
  ttest_IV_I <- t.test(x = meta.$value[meta.$StageNow_Simplified == "IV"],
                       y = meta.$value[meta.$StageNow_Simplified == "I"],
                       alternative="two.sided", var.equal = TRUE)[["p.value"]]
  cat("P.Value ",gene, " T-Test IV vs I:", ttest_IV_I, "\n")
  
  ttest_II_I <- t.test(x = meta.$value[meta.$StageNow_Simplified == "II"],
                       y = meta.$value[meta.$StageNow_Simplified == "I"],
                       alternative="two.sided", var.equal = TRUE)[["p.value"]]
  cat("P.Value ",gene, " T-Test II vs I:", ttest_II_I, "\n")
  
  
  # Plots Superplot based on biological replicate averages
  p = meta. %>% ggplot(aes(x=StageNow_Simplified, y= value, color=factor(sample_id))) +
    
    # Adds individu al data points
    geom_beeswarm(cex=2) +
    
    # Adds mean values as bars
    stat_summary(data = TotalAverages, fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = 0.25, color = "black") +
    
    # Adds error bars
    stat_summary(data = ReplicateAverages, fun.data = mean_se,
                 geom = "errorbar", width = 0.1, color = "black", size= 1) +
      
    #Adds Replicative averages as points (argument "cex" can be used to spread the data points if the averages are close together)
    geom_beeswarm(aes(fill = sample_id), color = "black", data=ReplicateAverages, size=5, pch = 23) +
    
    #Cosmetics and labeling
    theme_bw() + theme(axis.line = element_line(size = 1, colour = "black"),
                       legend.position = "none",
                       axis.ticks = element_line(size = 1, color = "black"), 
                       axis.ticks.length = unit(2, "mm"),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(), 
                       panel.border = element_blank()) +
    xlab("")
  
  pdf(file.path(output, paste0("SuperPlot_",gene,".pdf")))
  print(p)
  dev.off()
}



