# loading the data
setwd("/Users/yun-tsan.chang/Desktop/")
data <- read.csv("upset_sample_stage.csv", header = TRUE, row.names = NULL) 
names <- make.unique(as.character(data$X))
rownames(data) <- names
data <- data[, -1] # get rid of old names


install.packages("UpSetR")
library(UpSetR)


upset(data, 
      nsets = 13, sets = c("stage_IIB", "stage_IIA", "stage_IB", "stage_IA", "B1", "B0", "M0", "N2", "N1", "N0", "T3", "T2", "T1"), 
      keep.order = T, 
      matrix.color = "black", 
      main.bar.color = "black", 
      sets.bar.color = "black", 
      mainbar.y.label = "# patients meeting combined criteria",
      sets.x.label = "# patients meeting criterion",
      order.by = "degree", 
      decreasing = T,
      text.scale = c(1.5, 1.2, 1.2, 1.2, 1.2, 0.75),
      mb.ratio = c(0.65, 0.35))
