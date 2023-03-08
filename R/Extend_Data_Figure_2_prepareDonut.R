# Prepare a donut plot
library(ggplot2)
library(dplyr)

### Functions
# Function rounds values such that they still add up to 100%
round_preserve_sum <- function(x, digits = 0) {
  up = 10 ^ digits
  x = x * up
  y = floor(x)
  indices = tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] = y[indices] + 1
  y / up
}

# Parse command line arguments 
args = commandArgs(trailingOnly = TRUE)

# Grouping data 
group.data <- read.table(args[1], sep = "\t", header = T)

# Add the unidentified cells to as alpha beta 
metadata = read.csv(args[2])

sample_id = paste(unlist(strsplit(basename(args[1]), split = "_"))[1:3], collapse = "_")

# Find the row in metadata corresponding to input file name
row.match <- rep(F, nrow(metadata))
for (i in 1:nrow(metadata)) {
  row.match[i] = sample_id == metadata[i,1]
}


if (args[3] == "all") {
  
  # Reduce the observed counts from the total counts
  unidentified.count = metadata[row.match, 2] - nrow(group.data)
  
  counts <- as.data.frame(table(group.data$Group))
  
  if ("unidentified" %in% counts$Var1) {
    counts$Freq[counts$Var1 == "unidentified"] = counts$Freq[counts$Var1 == "unidentified"]  + unidentified.count
  } else {
    df = data.frame("Var1" = "unidentified" , "Freq" = unidentified.count)
    counts = rbind(counts, df)	
  }
  
  counts$Var1 = factor(counts$Var1, levels = c("main clone", "related to main clone", "bystander groups", "single bystanders", "unidentified"))  
  
  # Arrange according to factor levels 
  counts <- counts[order(counts$Var1),]
  
  counts$fraction = counts$Freq / sum(counts$Freq)
  counts$ymax = cumsum(counts$fraction)
  
  # Compute the bottom of each rectangle
  counts$ymin = c(0, head(counts$ymax, n=-1))
  counts$labelPosition <- (counts$ymax + counts$ymin)/2
  
  
  counts$perc <- round_preserve_sum(100 * counts$fraction, 1)
  counts$lab = paste0(counts$Freq, " (", counts$perc ,"%)")
  
  counts = counts[counts$Freq > 0,]
  
  gg <- ggplot(counts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
    geom_rect() +
    geom_text(x = 3.5, aes(y=labelPosition, label=lab), size=8) +
    coord_polar(theta="y") +
    xlim(c(2, 4))  +
    theme_void() +
    scale_fill_manual(values = c("main clone" = "#EC4C4C", 
                                 "related to main clone" = "#FDC37E",
                                 "bystander groups" = "#D29FF8",
                                 "single bystanders" = "#4F4FC3",
                                 "unidentified" = "grey")) + 
    theme(legend.position = "none") 
  
  print( gsub("_group.tsv", "_donut.pdf", args[1]))
  ggsave(plot = gg, filename = gsub("_group.tsv", "_all_donut.pdf", args[1]), dpi = 800)
  
} else {
  
  # Remove cells with unidentified TCR
  group.data <- group.data %>% dplyr::filter(Group != "unidentified") 
  group.data$Group = factor(group.data$Group, levels = c("main clone","related to main clone","bystander groups","single bystanders"))
  
  counts <- as.data.frame(table(group.data$Group))
  
  # Arrange according to factor levels 
  counts <- counts[order(counts$Var1),]
  
  counts$fraction = counts$Freq / sum(counts$Freq)
  counts$ymax = cumsum(counts$fraction)
  
  # Compute the bottom of each rectangle
  counts$ymin = c(0, head(counts$ymax, n=-1))
  counts$labelPosition <- (counts$ymax + counts$ymin)/2
  
  counts$perc <- round_preserve_sum(100 * counts$fraction, 1)
  counts$lab = paste0(counts$Freq, " (", counts$perc ,"%)")
  
  saveRDS(counts, "test.rds")
  
  counts = counts[counts$Freq > 0,]
  
  gg <- ggplot(counts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
    geom_rect() +
    geom_text(x = 3.5, aes(y=labelPosition, label=lab), size=8) +
    coord_polar(theta="y") +
    xlim(c(2, 4))  +
    theme_void() +
    scale_fill_manual(values = c("main clone" = "#EC4C4C", 
                                 "related to main clone" = "#FDC37E",
                                 "bystander groups" = "#D29FF8",
                                 "single bystanders" = "#4F4FC3")) + 
    theme(legend.position = "none") 
  
  print( gsub("_group.tsv", "_donut.pdf", args[1]))
  ggsave(plot = gg, filename = gsub("_group.tsv", "_identified_donut.pdf", args[1]), dpi = 800)
}


