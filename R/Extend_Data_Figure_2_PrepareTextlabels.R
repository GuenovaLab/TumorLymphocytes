# Prepare a donut plot
library(ggplot2)
library(dplyr)
library(stringr)
library(Cairo)

# Parse command line arguments 
args = commandArgs(trailingOnly = TRUE)

# Read in meta data and sample name
metadata = read.csv(args[1]) 
sample_id = args[2]

# Find the row in metadata corresponding to input file name
row.match <- rep(F, nrow(metadata))
for (i in 1:nrow(metadata)) {
  row.match[i] = sample_id == metadata[i,1]
}

# Find the malignant chains 
malignant.chain.alpha = as.character(metadata[row.match, 4])
malignant.chain.beta = as.character(metadata[row.match, 5])

label.title = "TCR associated to main clone"
label.a = paste0("\u03b1", "-chain ", ": ", malignant.chain.alpha)
label.b = paste0("\u03b2", "-chain ", ": ", malignant.chain.beta)

#gg <- ggplot() + 
#   annotate(geom = "text", x = 260, y = 0.9 , size = 4, label = label.title, colour = "black") +
#   annotate(geom = "text", x = 260, y = 0.7 , size = 4, label = label.a, colour = "#EC4C4C") +
#   annotate(geom = "text", x = 430, y = 0.5,  size = 4, label = label.b, colour = "#EC4C4C") +
#   xlim(c(0, 1000)) + ylim(c(0, 2)) +
#   theme_void()


df = data.frame(x = c(1, 1, 1), y = c(0.9, 0.7, 0.5), label = c(label.title, label.a, label.b), colour = c("title","chain","chain"))

gg <- ggplot() + geom_text(data = df , mapping = aes(x = x, y = y, label = label, colour = colour), size = 4, hjust = 0) +
   xlim(c(0, 1000)) + ylim(c(0, 2)) + scale_color_manual(values = c("title" = "black", "chain" = "#EC4C4C")) +
   theme_void() +
   theme(legend.position = "none")

ggsave(plot = gg, filename = file.path(args[3], paste0(sample_id, "_labels.pdf")), device = cairo_pdf, dpi = 800)

