# loading packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggpubr")

library("ggpubr")


setwd("M:/DER/RECHERCHE/GUENOVA_LAB/Yun-Tsan Chang/CTCL_research/project_scRNA_seq/GeneVia/raw_data/ALL/four_definitions/skin and blood/")

# loading the data
data <- read.csv("all_healthy_raw_counts_TCR_without_Tc_group.csv", header = TRUE, row.names = NULL) 
names <- make.unique(as.character(data$labels))
rownames(data) <- names
data <- data[, -1] # get rid of old names


head(data, 6)

ggscatter(data, x = "IL32", y = "HLA.ABC",
          add = "reg.line",         # Add regression line
          conf.int = TRUE,          # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
         ) +
  stat_cor(method = "pearson")


# Change the point shape
ggscatter(data, x = "IL32", y = "HLA.A",
          shape = 18)

# show shapes
show_point_shapes()



ggscatter(data, x = "IL32", y = "HLA.ABC",
          add = "reg.line",                         # Add regression line
          conf.int = FALSE,                          # Add confidence interval
          color = "definitions", palette = "jco"         # Color by groups "tumor"
) +
stat_cor(aes(color = definitions))           # Add correlation coefficient




# Association between IL32 and HLA_A
# Color points by dataset
# Add correlation coefficient by dataset
ggscatter(data, x = "IL32", y = "HLA.ABC", size = 1.5, 
          rug = TRUE,                                # Add marginal rug
          color = "definitions", palette = "jco") +
  stat_cor(aes(color = definitions), method = "spearman")



# Facet/split by data set, add regression line and confidence interval:

data1 <- data                                                 # Replicate original data
data1$definitions <- factor(data1$definitions,                                    # Change ordering manually
                  levels = c("Healthy", "SingleBystanders", "BystanderGroups", "RelateMainClone", "MainClone"))


P1 <- ggscatter(data1, x = "IL32", y = "HLA.ABC", size = 1,
          color = "definitions", palette = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"),
          facet.by = "definitions", #scales = "free_x",
          add = "reg.line", conf.int = F, ncol = 5) +
  stat_cor(aes(color = definitions), method = "spearman", label.y = 4.5, size = 4.5)


P2 <- ggscatter(data1, x = "IL2RG", y = "HLA.ABC", size = 1,
                color = "definitions", palette = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"),
                facet.by = "definitions", #scales = "free_x",
                add = "reg.line", conf.int = F, ncol = 5) +
  stat_cor(aes(color = definitions), method = "spearman", label.y = 4.5, size = 4.5)

P3 <- ggscatter(data1, x = "IL7R", y = "HLA.ABC", size = 1,
                color = "definitions", palette = c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"),
                facet.by = "definitions", #scales = "free_x",
                add = "reg.line", conf.int = F, ncol = 5) +
  stat_cor(aes(color = definitions), method = "spearman", label.y = 4.5, size = 4.5)

CombinePlots(plots = list(P1, P2, P3), ncol = 1)

