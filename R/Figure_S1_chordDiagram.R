install.packages("circlize")
library(circlize)

setwd("M:/DER/RECHERCHE/GUENOVA_LAB/Project_Yun-Tsan_ADCC_paper/result/Chord gram/")

# loading the data
data <- read.csv("Healthy_T_cells_chord.csv", header = TRUE, row.names = NULL) 
names <- make.unique(as.character(data$X))
rownames(data) <- names
data <- data[, -1] # get rid of old names

data <- as.matrix(data)
data


nm = unique(unlist(dimnames(data)))
nm


# library(magrittr)
nm_1 <- gsub("_", "", nm)
nm_1

group = structure(gsub("[0-9]*", "", nm_1), names = nm)
group

circos.clear()

circos.par(start.degree = 80)

# grid.col <- setNames(rainbow(length(unlist(dimnames(data)))), union(rownames(data), colnames(data)))
# grid.col <-  c(TRBV20.1.TRBJ2.1 = "red", TRBV26_1_TRBJ53 = "red")


chordDiagram(t(data), group = group, annotationTrack = c("grid", "axis"), grid.col = c("gray40"))

# chordDiagram(t(data), annotationTrack = c("grid", "axis"), grid.col = grid.col  ,row.col = c("#FF000080", "#00FF0010", "#0000FF10"))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1]+3, CELL_META$sector.index,
              facing = "clockwise", cex = 0.6 ,niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)

title("healthy", cex = 0.6)



circos.clear()



?circos.text


?chordDiagram()


