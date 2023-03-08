# Prepare a donut plot
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)
library(purrr)

# Parse command line arguments 
args = commandArgs(trailingOnly = TRUE)

# Grouping data 
group.data <- read.table(args[1], sep = "\t", header = T)
metadata = read.csv(args[2]) 

sample_id = paste(unlist(strsplit(basename(args[1]), split = "_"))[1:3], collapse = "_")

# Find the row in metadata corresponding to input file name
row.match <- rep(F, nrow(metadata))
for (i in 1:nrow(metadata)) {
  row.match[i] = sample_id == metadata[i,1]
}

# Find the malignant chains 
malignant.chain.alpha = as.character(metadata[row.match, 4])
malignant.chain.beta = as.character(metadata[row.match, 5])

#### Alpha

#  Counts 
tcra.count <- as.data.frame(table(group.data$TCRA)) %>% dplyr::filter(Var1 != "beta") 
tcrb.count <- as.data.frame(table(group.data$TCRB)) %>% dplyr::filter(Var1 != "alpha") 

# Annotate based on the malignant clone 
chain.type.tcra <- ifelse(tcra.count$Var1 == malignant.chain.alpha, "malignant", "non-malignant")

# Add the chain type
tcra.count$chain.type <- factor(as.character(chain.type.tcra), levels = c("non-malignant", "malignant")) 

# Shorten the chain names
alpha_chains = unlist(map(strsplit(as.character(tcra.count$Var1), split = "_"), function(x){paste(x[c(1,length(x))], collapse = "...")}))
tcra.count$Short <- as.factor(paste0(alpha_chains, "." , seq(1, length(alpha_chains))))

#### Beta

# Annotate based on the malignant clone 
chain.type.tcrb <- ifelse(tcrb.count$Var1 == malignant.chain.beta, "malignant", "non-malignant")

# Add the chain type
tcrb.count$chain.type <- factor(as.character(chain.type.tcrb), levels = c("non-malignant","malignant"))

# Shorten the chain names
beta_chains = unlist(map(strsplit(as.character(tcrb.count$Var1), split = "_"), function(x){paste(x[c(1,length(x))], collapse = "...")}))
tcrb.count$Short <- as.factor(paste0(beta_chains, ".", seq(1, length(beta_chains))))

# Scale the barplots according to maximum height 
max.height = max(c(tcra.count$Freq, tcrb.count$Freq))

# We want to make the bars equally wide. We will add bars to 38
if ((37 - length(tcra.count$Var1) > 0)) {
   fill.a.df = data.frame(Var1 = as.character(paste0("C", seq(0, 37 - length(tcra.count$Var1)))) , 
                          Freq = rep(0, 38 - length(tcra.count$Var1)), chain.type = "non-malignant", 
                          Short = as.factor(seq(0, 37 - length(tcra.count$Var1))))
   fill.a.df$chain.type <- factor(fill.a.df$chain.type, levels = c("non-malignant","malignant"))
   tcra.count = rbind(tcra.count, fill.a.df)
}

if ((37 - length(tcrb.count$Var1) > 0)){
   fill.b.df = data.frame(Var1 = as.character(paste0("C", seq(0, 37 - length(tcrb.count$Var1)))) , 
                          Freq = rep(0, 38 - length(tcrb.count$Var1)), 
                          chain.type = "non-malignant",
                          Short = as.factor(seq(0, 37 - length(tcrb.count$Var1))))
   fill.b.df$chain.type <- factor(fill.b.df$chain.type, levels = c("non-malignant","malignant"))
   tcrb.count = rbind(tcrb.count, fill.b.df)
} 

# Arrange by frequency
tcra.count$Short <- factor(tcra.count$Short, levels = tcra.count$Short[order(tcra.count$Freq, decreasing = T)])

# Arrange by frequency
tcrb.count$Short <- factor(tcrb.count$Short, levels = tcrb.count$Short[order(tcrb.count$Freq, decreasing = T)])

# Function for plotting the barplot 
plotBar <- function(counts, highlighted.counts = NULL, chain, breaks, malignant.chain){
  
  # Create horizontal lines representing the percentage for different frequencies
  idx = which(!duplicated(counts$perc))
  uniq.counts <- counts[idx,]
  hline.annot = data.frame(y = uniq.counts$Freq, lab = paste(uniq.counts$perc, "%"))
  
  # If no highlighted chains
  if (is.null(highlighted.counts)) {
    
    # If chain is alpha
    if (chain == "alpha") {
      
      gg <- ggplot() + 
        geom_col(data = counts, mapping = aes(x = Short, y = Freq, fill = chain.type)) + 
        ggtitle(paste0("\u03b1", "-chains")) +
        theme_classic() +
        theme(title = element_text(size=20),
              plot.title = element_text(hjust = 0.5),
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_blank(),
              axis.text.y = element_text(size=20, colour = "black"),
              legend.position = "none") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max.height + max.height/8 + 1), breaks = breaks) +
        scale_fill_manual(values = c("malignant" = "#EC4C4C",
                                     "non-malignant" = "black")) +
        geom_hline(data = hline.annot, mapping = aes(yintercept = y), linetype = "dashed") +
        geom_text(data = hline.annot, mapping = aes(y = y, x = 30, label = lab))
      
    # Chain is beta  
    } else {
      
      gg <- ggplot() + 
        geom_col(data = counts, mapping = aes(x = Short, y = Freq, fill = chain.type)) + 
        ggtitle(paste0("\u03b2", "-chains")) +
        theme_classic() + 
        theme(title = element_text(size=20),
              plot.title = element_text(hjust = 0.5),
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_blank(),
              axis.text.y = element_text(size=20, colour = "black"),
              legend.position = "none") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max.height + max.height/8 + 1), breaks = breaks) +
        scale_fill_manual(values = c("malignant" = "#EC4C4C",
                                     "non-malignant" = "black")) +
        geom_hline(data = hline.annot, mapping = aes(yintercept = y), linetype = "dashed") +
        geom_text(data = hline.annot, mapping = aes(y = y, x = 30, label = lab))
      
    }
    
  # Highlighted chains exist
  } else {
    
    if (chain == "alpha") {
    
      gg <- ggplot() + 
        geom_col(data = counts, mapping = aes(x = Short, y = Freq, fill = chain.type, colour = chain.type)) + 
        geom_segment(data = highlighted.counts, mapping = aes(x = seq(1, nrow(highlighted.counts)), 
                                                xend = seq(1, nrow(highlighted.counts)) + 10, 
                                                y = Freq, 
                                                yend = Freq + max.height/8,colour = chain.type), 
                     size = 1) +
        geom_text(data = highlighted.counts, mapping = aes(x = seq(1, nrow(highlighted.counts)) + 10, y = Freq + max.height/8, label = Short, colour = chain.type)) +
        ggtitle(paste0("\u03b1", "-chains")) +
        theme_classic() +
        theme(title = element_text(size=20),
            plot.title = element_text(hjust = 0.5),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line = element_blank(),
            axis.text.y = element_text(size=20, colour = "black"),
            legend.position = "none") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max.height + max.height/8 + 1), breaks = breaks) +
        scale_fill_manual(values = c("malignant" = "#EC4C4C",
                                   "non-malignant" = "black")) +
        scale_colour_manual(values = c("malignant" = "#EC4C4C",
                                     "non-malignant" = "black")) +
        geom_hline(data = hline.annot, mapping = aes(yintercept = y), linetype = "dashed") +
        geom_text(data = hline.annot, mapping = aes(y = y, x = 30, label = lab))
      
      # Chain is beta
      } else {
        
        gg <- ggplot() + 
          geom_col(data = counts, aes(x = Short, y = Freq, fill = chain.type, colour = chain.type)) + 
          geom_segment(data = highlighted.counts, mapping = aes(x = seq(1, nrow(highlighted.counts)), 
                                                  xend = seq(1, nrow(highlighted.counts)) + 10, 
                                                  y = Freq, 
                                                  yend = Freq + max.height/8, colour = chain.type), 
                       size = 1) +
          geom_text(data = highlighted.counts, mapping = aes(x = seq(1, nrow(highlighted.counts)) + 10, y = Freq + max.height/8, label = Short, colour = chain.type)) +
          ggtitle(paste0("\u03b2", "-chains")) +
          theme_classic() + 
          theme(title = element_text(size=20),
                plot.title = element_text(hjust = 0.5),
                axis.title = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line = element_blank(),
                axis.text.y = element_text(size=20, colour = "black"),
                legend.position = "none") +
          scale_y_continuous(expand = c(0, 0), limits = c(0, max.height + max.height/8 + 1), breaks = breaks) +
          scale_fill_manual(values = c("malignant" = "#EC4C4C",
                                       "non-malignant" = "black")) +
          scale_colour_manual(values = c("malignant" = "#EC4C4C",
                                         "non-malignant" = "black")) +
          geom_hline(data = hline.annot, mapping = aes(yintercept = y), linetype = "dashed") +
          geom_text(data = hline.annot, mapping = aes(y = y, x = 30, label = lab))
      }
      
  }
  
  malignant.annot = NULL
  
  if (!(malignant.chain %in% highlighted.counts$Var1)) {
    malignant.annot = counts[counts$Var1 == malignant.chain,]
    malignant.short = as.character(malignant.annot$Short)
    malignant.y = malignant.annot$Freq
  } else {
    malignant.annot = NULL
  }
  
  if (!is.null(malignant.annot)) {
    gg <- gg + 
         annotate(geom= "segment", x = which(malignant.short == levels(counts$Short)), 
                                   xend = which(malignant.short == levels(counts$Short)) + 10, 
                                   y = malignant.y, 
                                   yend =  malignant.y + max.height/8, size = 1, colour = "#EC4C4C") +
         annotate(geom = "text", x = which(malignant.short == levels(counts$Short)) + 10, y = malignant.y + max.height/8, label = malignant.short, colour = "#EC4C4C")
  }
  return(gg)
}


########## Prepare plot for alpha chains 

# Determine the breaks
breaks.alpha = sort(c(seq(0, max.height, by = 5), max(tcra.count$Freq)))

# Calculate the percentages corresponding to breaks 
tcra.count$perc <- round(100 * tcra.count$Freq/sum(tcra.count$Freq), 1)

# Generate annotation
generateAnnotation <- function(counts) {
  
  # Select only chains which correpond over 10 percent of the total
  # identified chains
  annot = counts[counts$perc > 10,]
  
  # Test if there are any chains meeting the condition
  if (nrow(annot) > 0) {
    
    # Drop levels
    annot$Short <- droplevels(annot$Short) 
    
    # Order based on levels
    annot <- annot[order(annot$Short),]
    
    # Remove the suffix (. + number)
    suff.rm.a = gsub("[.]\\d+","", as.character(annot$Short))
    
    # If duplicates exist we need to change those to full names
    if (any(duplicated(suff.rm.a))) {
      
      
      # Duplicates
      dups = as.vector(suff.rm.a[duplicated(suff.rm.a)])
      
      print(annot$Short)
      print(gsub("[.]\\d+","", as.character(annot$Short)))
      
      
      # Assign full names
      findFullName = function(x, counts){
        Var1 = as.character(counts$Var1)
        Short = as.character(counts$Short)
        return(Var1[x == Short])
      }
      
      suff.rm.a <- unlist(map_if(as.character(annot$Short), gsub("[.]\\d+","", as.character(annot$Short)) %in% dups, findFullName, counts = counts))
      suff.rm.a = gsub("[.]\\d+","",suff.rm.a)
      annot$Short = suff.rm.a
      
      annot$Short = factor(suff.rm.a, levels = suff.rm.a)
      
    } else {
      
      # No duplicates 
      annot$Short = factor(suff.rm.a, levels = suff.rm.a)
    }
    
    
  } else {
    # We need to add the malignant chain
    annot = NULL
  }
  return(annot)
}

annotation.alpha = generateAnnotation(tcra.count)

saveRDS(tcra.count, file = gsub("_group.tsv", "_TCRA_barplot.rds", args[1]))

if (!is.null(annotation.alpha)) {
  
  gg.1 <- plotBar(counts = tcra.count, 
                  highlighted.counts = annotation.alpha,
                  chain = "alpha",
                  breaks = breaks.alpha,
                  malignant.chain = malignant.chain.alpha)
  
} else {
  gg.1 <- plotBar(counts = tcra.count, 
                  highlighted.counts = NULL,
                  chain = "alpha",
                  breaks = breaks.alpha,
                  malignant.chain = malignant.chain.alpha) 
}

ggsave(plot = gg.1, filename = gsub("_group.tsv", "_TCRA_barplot.pdf", args[1]), device = cairo_pdf, dpi = 800)

########## Prepare plot for beta chains 

# Determine the breaks
breaks.beta  = sort(c(seq(0, max.height, by = 5), max(tcrb.count$Freq)))

# Calculate the percentages
tcrb.count$perc <- round(100 * tcrb.count$Freq/sum(tcrb.count$Freq), 1)

# Generate annotation
annotation.beta = generateAnnotation(tcrb.count)

saveRDS(tcrb.count, file = gsub("_group.tsv", "_TCRB_barplot.rds", args[1]))

if (!is.null(annotation.beta))  {
  
  gg.2 <- plotBar(counts = tcrb.count, 
                  highlighted.counts = annotation.beta,
                  chain = "beta",
                  breaks = breaks.beta,
                  malignant.chain = malignant.chain.beta)
  
} else {
  gg.2 <- plotBar(counts = tcrb.count, 
                  highlighted.counts = NULL,
                  chain = "beta",
                  breaks = breaks.beta,
                  malignant.chain = malignant.chain.beta) 
}

ggsave(plot = gg.2, filename = gsub("_group.tsv", "_TCRB_barplot.pdf", args[1]), device = cairo_pdf, dpi = 800)
