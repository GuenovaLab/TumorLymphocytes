install.packages("ggmosaic")


library(Seurat)
library(ggplot2)
library(ggmosaic)


# load data
setwd("/Users/yun-tsan.chang/Desktop/")
data <- read.csv("clonality_mosaic.csv", header = TRUE, row.names = NULL) 

count <- table(data$patient_1, data$TRBV)
mosaicplot(count)


mosaicplot(count, main = "patients and TRBVs",
           xlab = "patients",
           ylab = "TRBVs",
           las = 1)















