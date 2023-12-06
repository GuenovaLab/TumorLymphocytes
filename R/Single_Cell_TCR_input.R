# Loading scRNA TCR-containing T cells
input_dir = "/home/localadmin/Téléchargements/GSE224449_processed_data/"

# Read in matrices
files = list.files(input_dir, pattern = "MF|Healthy", full.names = T)
names(files) = gsub(".csv","",basename(files))
list = lapply(files, read.csv)
names(list) = NULL
all_matrices = do.call("cbind", list)
colnames(all_matrices) = gsub("^X", "", colnames(all_matrices))
genes = all_matrices[,1]

# Read in Cell information
Cell_information = read.csv(file.path(input_dir, "Cell information.csv"))
Cell_information$labels = gsub("-", ".", Cell_information$labels)

# Select cells with a TCR detected or CD8+ T cells
Cell_information_TCR_containing = Cell_information[which(Cell_information$definitions != "undefined" & Cell_information$definitions != "CD8 T"),]

# Average over different transcripts when it is the case
all_matrices = all_matrices[,grep("Gene.HGNC", colnames(all_matrices), invert = TRUE)]
all_matrices$gene = genes

all_matrices = all_matrices %>% group_by(gene) %>% summarise(across(everything(), mean))
genes = all_matrices$gene
rownames(all_matrices) = genes
all_matrices = all_matrices[,setdiff(colnames(all_matrices),"gene")]
colnames(all_matrices) = gsub(".X","_",colnames(all_matrices))


# Select cells
all_matrices_TCR = all_matrices[,match(Cell_information_TCR_containing$labels, colnames(all_matrices))]

# Save metadata
write.csv(Cell_information_TCR_containing, "all_healthy_raw_counts_TCR_without_Tc_group.csv", header = TRUE, row.names = NULL)

# Save matrix
write.csv(all_matrices_TCR, "all_healthy_raw_counts_TCR_without_Tc.csv", header = TRUE, row.names = NULL) 

