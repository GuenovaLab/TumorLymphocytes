

setwd("/Users/yun-tsan.chang/Desktop/")

# loading the data
data <- read.csv("01_WaG_SE_TCR.csv", header = TRUE, row.names = NULL) 



library(tidyr)


data_1 <- pivot_wider(data, names_from = "to", values_from = "value")

# loading the data
names <- make.unique(as.character(data_1$from))
rownames(data_1) <- names

data_1

data_2 <- as.matrix(data_1)


write.csv(data_2, file="01_WaG_SE_TCR_chord.csv", row.names = TRUE, col.names = TRUE)





