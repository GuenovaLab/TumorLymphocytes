install.packages("viridis")

# library
library(ggplot2)
library(viridis)
library(hrbrthemes)



setwd("M:/DER/RECHERCHE/GUENOVA_LAB/Project_Yun-Tsan_ADCC_paper/raw data and R code/Figure_6/6A_ridge plot/")

# loading the data
data <- read.csv("barplot_cluster_stage.csv", header = TRUE, row.names = NULL) 



# grouped boxplot
ggplot(data, aes(x=cluster, y=percentage, fill=StageNow)) +
  scale_fill_viridis(discrete = T, option = "D") +
  geom_boxplot() +
  ggtitle("Stage at end of follow-up") +
  theme_ipsum()+
  xlab("")


