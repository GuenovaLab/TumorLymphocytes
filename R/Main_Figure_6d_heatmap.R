install.packages("pheatmap")
install.packages("heatmaply")
install.packages("gplot")

setwd ("M:/DER/RECHERCHE/GUENOVA_LAB/Yun-Tsan Chang/CTCL_research/nanostring/nanoString_FFPEvsCRYO/")

# load data
my_data <- read.csv("20201210_NormalizedData_CRYO_percentage_1.csv", header = TRUE, row.names = 1)

my_data <- as.matrix(my_data)

heatmap(my_data, scale = "none")
        
?heatmap()

# install.packages("gplots")
library("gplots")
heatmap.2(my_data, scale = "none", col = greenred(50),
          trace = "none", density.info = "none")

library("pheatmap")
pheatmap(my_data, cutree_cols = 2, annotation_names_row = T)






library(heatmaply)
mtcars_2 <- percentize(mtcars)
heatmaply(my_data, k_row = 4, k_col = 2)


?pheatmap()
