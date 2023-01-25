### Prepares clonograms 
### Takes a count matrix (adjacency matrix as input) and the TCRs related to the malignant clone 

# Features to add 
# Add edge from zero alpha to zero beta

library("network")
library("ggnetwork")
library("sna")
library("ggplot2")
library("dplyr")
library("ggraph")
library("tibble")
library("Cairo")

###### Functions 

# Function for the project 

# Define the group for each each. If an edge invoves a TCR related to the 
# malignant clone set as malignant else set as non-clonal bystander 
groupEdge <- function(x, malignant.TCR){
  # Vertices 
  v.1 <- x[1]
  v.2 <- x[2]
  
  # Find the chains corresponding to these vertices
  vertice.names <- attr(edge.list.matrix, "vnames") 
  v.1.name = vertice.names[v.1]
  v.2.name = vertice.names[v.2]
  
  # Finally test if either of the vertices correspond to 
  # malignant clone 
  if ((v.1.name == malignant.TCR["TCRA"]) | (v.1.name == malignant.TCR["TCRB"])) {
    return("malignant")
  } else if ((v.2.name == malignant.TCR["TCRA"]) | (v.2.name == malignant.TCR["TCRB"])){
    return("malignant")
  } else {
    return("non-malignant bystander")
  }
  
}

# Assign the end vertices corresponding each edge
assignEndVertices <- function(network.obj, edge.list.matrix){
  
  # Select only rows describing the true edges between vertices (no loops)
  idx.no.loop <- which(!(network.obj$x == network.obj$xend & network.obj$y == network.obj$yend))
  network.subset <- network.obj[idx.no.loop,]
  
  # For each edge ID find the corresponding row in edge.list.matrix and retrieve the end vertex
  getEndVertex <- function(row, edge.list.matrix){
    id <- as.numeric(as.character(row["edge.ID"]))
    end.vertex.id <- edge.list.matrix[id, 2]
    vertice.names <- attr(edge.list.matrix, "vnames")
    return(vertice.names[end.vertex.id])
    
  }
  
  end.vertices <- apply(network.subset, 1, getEndVertex, edge.list.matrix = edge.list.matrix)
  df.out = data.frame(network.subset, end.vertices = end.vertices)
  
  # Assign the self-loops end.vertices to to vertex.names
  idx.loop <- which(network.obj$x == network.obj$xend & network.obj$y == network.obj$yend)
  network.subset.loop <- network.obj[idx.loop,]
  df.loop.out <- data.frame(network.subset.loop, end.vertices = network.subset.loop$vertex.names)
  
  # Merge
  df.out = rbind(df.out, df.loop.out)
  
  return(df.out)
}


# Organise in four columns
clonogram_layout <- function(network.obj, start.x, start.y, row.dist, col.dist, group.1, group.2) {
  
  # Find the largest group. The layout will be made based on this group
  if (length(group.1) > length(group.2)) {
    
    # Group 1 is larger
    start.x.coord.group.1 <- rep(start.x + col.dist, length(group.1))
    start.y.coord.group.1 <- seq(start.y, by = row.dist, length.out = length(group.1))
    
    # Adjust Group 2 coordinates based on Group 1 
    # the nodes are aligned such that the nodes are equally spaced and 
    # the farthest nodes are equally distant from the group 1 farthest nodes in y-direction
    start.x.coord.group.2 <- rep(start.x + 2 * col.dist, length(group.2)) 
    
    adjusted.row.dist = (start.y.coord.group.1[length(start.y.coord.group.1)] - start.y.coord.group.1[1])/(length(group.2) + 1)
    start.y.coord.group.2 <- seq(adjusted.row.dist + start.y.coord.group.1[1], by = adjusted.row.dist, length.out = length(group.2)) 
    
    # Adjust the alpha zero coordinates
    alpha.x = start.x
    alpha.y = (start.y.coord.group.1[length(start.y.coord.group.1)] - start.y.coord.group.1[1])/2 +  start.y.coord.group.1[1]
    beta.x = start.x + 3 * col.dist    
    beta.y = (start.y.coord.group.1[length(start.y.coord.group.1)] - start.y.coord.group.1[1])/2 + start.y.coord.group.1[1]
    
  } else if ((length(group.2) > length(group.1))){
    # Group 2 is larger
    start.x.coord.group.2 <- rep(start.x + 2 * col.dist, length(group.2))
    start.y.coord.group.2 <- seq(start.y, by = row.dist, length.out = length(group.2))
    
    # Adjust Group 2 coordinates based on Group 1 
    # the nodes are aligned such that the nodes are equally spaced and 
    # the farthest nodes are equally distant from the group 1 farthest nodes in y-direction
    start.x.coord.group.1 <- rep(start.x + col.dist, length(group.1)) 
    
    adjusted.row.dist = (start.y.coord.group.2[length(start.y.coord.group.2)] - start.y.coord.group.2[1])/(length(group.1) + 1)
    start.y.coord.group.1 <- seq(adjusted.row.dist + start.y.coord.group.2[1], by = adjusted.row.dist, length.out = length(group.1)) 
    
    # Adjust the alpha zero coordinates
    alpha.x = start.x
    alpha.y = (start.y.coord.group.2[length(start.y.coord.group.2)] - start.y.coord.group.2[1])/2 + start.y.coord.group.2[1]
    beta.x = start.x + 3 * col.dist    
    beta.y = (start.y.coord.group.2[length(start.y.coord.group.2)] - start.y.coord.group.2[1])/2 + + start.y.coord.group.2[1]
  }
  else {
    # Groups are equal in size
    start.x.coord.group.1 <- rep(start.x + col.dist, length(group.1))
    start.x.coord.group.2 <- rep(start.x + 2 * col.dist, length(group.2))  
    
    start.y.coord.group.1 <- seq(start.y, by = row.dist, length.out = length(group.1))
    start.y.coord.group.2 <- seq(start.y, by = row.dist, length.out = length(group.2))
    
    # Adjust the alpha zero coordinates
    alpha.x = start.x
    alpha.y = start.y.coord.group.1[length(start.y.coord.group.1)]/2 + start.y.coord.group.1[1]
    beta.x = start.x + 3 * col.dist    
    beta.y = start.y.coord.group.1[length(start.y.coord.group.1)]/2 + start.y.coord.group.1[1]
  }
  
  # Store the new coordinates into a dataframe
  group.1.new.coords <- data.frame(x = start.x.coord.group.1, 
                                   y = start.y.coord.group.1)
  rownames(group.1.new.coords) <- group.1
  
  group.2.new.coords <- data.frame(x = start.x.coord.group.2, 
                                   y = start.y.coord.group.2)
  rownames(group.2.new.coords) <- group.2
  
  # Create vectors to store the coordinates
  alpha.coordinates = c(alpha.x, alpha.y)
  beta.coordinates = c(beta.x, beta.y)
  
  # Update coordinates 
  updateCoord <- function(row,  group.1.new.coords, group.2.new.coords, alpha.coordinates, beta.coordinates){
    vertex.name = row["vertex.names"]
    end.vertex.name = row["end.vertices"]
    weight = row["weight"]
    group = row["group"]
    edge.ID = row["edge.ID"]
    group = row["group"]
    vertex.group = row["vertex.group"]
    
    # 
    if (vertex.name %in% rownames(group.1.new.coords)) {
      # Vertex is in group 1
      x.coord.new = group.1.new.coords$x[vertex.name == rownames(group.1.new.coords)]
      y.coord.new = group.1.new.coords$y[vertex.name == rownames(group.1.new.coords)]
    } else if (vertex.name %in% rownames(group.2.new.coords)){
      # Vertex is in group 2 
      x.coord.new = group.2.new.coords$x[vertex.name == rownames(group.2.new.coords)]
      y.coord.new = group.2.new.coords$y[vertex.name == rownames(group.2.new.coords)]      
    } else if (vertex.name == "alpha") {
      x.coord.new = alpha.coordinates[1]
      y.coord.new = alpha.coordinates[2]
    } else {
      # vertex name is beta
      x.coord.new = beta.coordinates[1]
      y.coord.new = beta.coordinates[2]
    }
    
    if (end.vertex.name %in% rownames(group.1.new.coords)) {
      # Vertex is in group 1
      x.coord.end.new = group.1.new.coords$x[end.vertex.name == rownames(group.1.new.coords)]
      y.coord.end.new = group.1.new.coords$y[end.vertex.name == rownames(group.1.new.coords)]
    } else if (end.vertex.name %in% rownames(group.2.new.coords)) {
      # Vertex is in group 2 
      x.coord.end.new = group.2.new.coords$x[end.vertex.name == rownames(group.2.new.coords)]
      y.coord.end.new = group.2.new.coords$y[end.vertex.name == rownames(group.2.new.coords)]      
    } else if (end.vertex.name == "alpha")  {
      x.coord.end.new = alpha.coordinates[1]
      y.coord.end.new = alpha.coordinates[2]
    } else {
      x.coord.end.new = beta.coordinates[1]
      y.coord.end.new = beta.coordinates[2]        
    }
    
    # Collect to output line 
    out.row <- c(x.coord.new, y.coord.new, as.character(vertex.name), 
                 x.coord.end.new, y.coord.end.new, edge.ID, group,
                 weight, as.character(end.vertex.name), vertex.group)
    
    return(out.row)
  }
  
  # Adjust the coordinates 
  clono.network <- t(apply(network.obj, 1, updateCoord, group.1.new.coords, group.2.new.coords, alpha.coordinates, beta.coordinates))
  
  # Set the column names
  colnames(clono.network) <- colnames(network.obj)
  
  # Convert data as as numeric
  clono.network = as.data.frame(clono.network)
  clono.network$x <- as.numeric(as.character(clono.network$x))
  clono.network$y <- as.numeric(as.character(clono.network$y))
  clono.network$weight <- as.numeric(as.character(clono.network$weight))
  clono.network$xend <- as.numeric(as.character(clono.network$xend))
  clono.network$yend <- as.numeric(as.character(clono.network$yend))
  
  return(clono.network)
  
}

assignEdgesToGroups <- function(edge.ls.matrix, group.data) {
  
  # Get vertice names
  vertice.names <- attr(edge.list.matrix, "vnames") 
  
  assignVerticeName <- function(row) {
    TCRB <- vertice.names[row[1]]
    TCRA <- vertice.names[row[2]]
    
    df <-data.frame(TCRA, TCRB)
    return(df)
  }
  
  # Convert edge list matrix into something that can be merged with group data
  temp <- apply(edge.ls.matrix, 1, assignVerticeName)
  temp <- do.call("rbind", temp)
  
  temp$ID <- 1:nrow(temp)
  
  # Harmonize the columns 
  # All dashes and underscores are converted into .
  temp[,1] <- gsub("_",".", gsub("-",".", as.character(temp[,1])))
  temp[,2] <- gsub("_",".", gsub("-",".", as.character(temp[,2])))
  
  group.data[,2] <- gsub("_",".", gsub("-",".", group.data[,2]))
  group.data[,3] <- gsub("_",".", gsub("-",".", group.data[,3]))
  
  # If alpha, beta edge (= unidentified chains) is not present in group.data
  # add it to ensure that this edge is correctly assigned to a group
  if (!(any((group.data$TCRA == "beta") & (group.data$TCRB == "alpha")))) {
    # We add 
    vect <- c(97,"beta", "alpha", "unidentified")
    names(vect) <- c("Row", "TCRA", "TCRB", "Group")
    df.add = as.data.frame(t(vect))
    group.data <- rbind(group.data, df.add) 
    group.data$Row <- as.vector(group.data$Row)
    group.data$TCRA <- as.vector(group.data$TCRA)
    group.data$TCRB <- as.vector(group.data$TCRB)
    group.data$Group <- as.vector(group.data$Group)
  }
  
  # Merge the group data
  df.merged <- left_join(temp, group.data, by = c("TCRA","TCRB"))
  
  df.selected <- df.merged %>% dplyr::select(ID, Group) 
  
  # Remove duplicate rows
  df.unique = unique(df.selected)
  
  #return(df.unique$Group)
  return(df.unique$Group)
}

######


# Parse command line arguments 
args = commandArgs(trailingOnly = TRUE)

### Read in data and preprocess
adjacency_matrix <- read.csv(args[1])
adjacency_matrix <- tibble::column_to_rownames(adjacency_matrix , "X")

# Grouping data 
group.data <- read.table(args[2], sep = "\t", header = T)

# Metadata 
metadata = read.csv(args[3]) 

sample_id = paste(unlist(strsplit(basename(args[1]), split = "_"))[1:3], collapse = "_")

print(sample_id)

# Find the row in metadata corresponding to input file name
row.match <- rep(F, nrow(metadata))
for (i in 1:nrow(metadata)) {
  #print(as.character(metadata[i,1]))
  #sample_id = paste(unlist(strsplit(args[1], split = "_"))[1:3], collapse = "_")
  row.match[i] = sample_id == metadata[i,1]
}

# Find the main chains
main_alpha = gsub("_",".", gsub("-",".", as.character(metadata[row.match, 4])))
main_beta = gsub("_",".", gsub("-",".",as.character(metadata[row.match, 5])))

### Generate the network
tcr.network <- network(adjacency_matrix, 
                    matrix.type = "bipartite",
                    ignore.eval = F,
                    names.eval = 'weight')

# Extract and edgelist matrix for adding End vertice
edge.list.matrix <- as.matrix(tcr.network, matrix.type="edgelist")

group <- assignEdgesToGroups(edge.list.matrix, group.data)

# Add edge ID which correponds to edgelist row 
set.edge.attribute(x = tcr.network, 
                   "edge.ID", 1:(network.edgecount(tcr.network)))

set.edge.attribute(x = tcr.network, 
                   "group", as.character(group))

# Generate the ggnetwork object 
tcr.network.gg <- ggnetwork(tcr.network, layout = "fruchtermanreingold")

# Find the end vertices for each edge
tcr.network.gg <- assignEndVertices(tcr.network.gg, edge.list.matrix)


# Add the vertex group 
vertex.names.mod = gsub("_",".", gsub("-",".", tcr.network.gg$vertex.names))
tcr.network.gg$vertex.group = ifelse(((vertex.names.mod == main_alpha) | (vertex.names.mod  == main_beta)), "main clone", "other")

print(head(tcr.network.gg))

### Change the layout of the graph to a clonogram 

# Define some constants 
start.clono.x = 0.5 # The x coordinate where to place the first node from the left
start.clono.y = 1   # The y coordinate where to place the first node from the bottom 
row.dist.clono = 1  # Row distance 
col.dist.clono = 2  # Col distance 

# Change the layout
tcr.network.gg <- clonogram_layout(tcr.network.gg, 
                       start.x = start.clono.x,
                       start.y = start.clono.y, 
                       row.dist = row.dist.clono, 
                       col.dist = col.dist.clono, 
                       group.1 = rownames(adjacency_matrix)[1:(nrow(adjacency_matrix) - 1)],
                       group.2 = colnames(adjacency_matrix)[1:(ncol(adjacency_matrix) - 1)])

### Prepare the figure 

# Define the label positions based on the nodes 
no.beta.label.x <- start.clono.x - 0.1
no.beta.label.y <- unique(tcr.network.gg$y[tcr.network.gg$vertex.names == "alpha"]) + 1
no.alpha.label.x <- start.clono.x + 3 * col.dist.clono + 0.1
no.alpha.label.y <- unique(tcr.network.gg$y[tcr.network.gg$vertex.names == "beta"]) + 1
alpha.chain.label.x = start.clono.x + col.dist.clono 
alpha.chain.label.y = start.clono.y - 1
beta.chain.label.x =  start.clono.x + 2 * col.dist.clono 
beta.chain.label.y = start.clono.y - 1

# Save object for checking 
saveRDS(tcr.network.gg, file = gsub("_adjmatrix.csv", "_clonogram.rds", args[1]))  

# Change the factor levels 
tcr.network.gg$group <- factor(tcr.network.gg$group, levels = c("main clone","related to main clone",
                                                                "bystander groups","single bystanders",
                                                                "unidentified"))

# Prepare the plot 
g <- ggplot(tcr.network.gg, aes(x = x, y = y, xend = xend, yend = yend)) +
     geom_edges(aes(size = weight, colour = group)) + 
     geom_nodes(aes(colour = vertex.group), size = 3) +
     theme_blank() + 
     theme(legend.position = "none", 
           legend.text = element_text(size=10),
           legend.box.background = element_rect(colour = "black")) +
     labs(color = "Clonality group") +
     annotate(geom = "text", size = 8, x = no.beta.label.x, y = no.beta.label.y, label = paste0("no ", "\u03b2")) +
     annotate(geom = "text", size = 8, x = no.alpha.label.x, y = no.alpha.label.y, label = paste0("no ", "\u03b1")) + 
     annotate(geom = "text", size = 8, x = alpha.chain.label.x, y = alpha.chain.label.y, label = paste0("\u03b1","-chains")) +
     annotate(geom = "text", size = 8, x = beta.chain.label.x, y = beta.chain.label.y, label = paste0("\u03b2", "-chains")) +
     scale_colour_manual(values = c("main clone" = "#EC4C4C", 
                                    "related to main clone" = "#FDC37E",
                                    "bystander groups" = "#D29FF8",
                                    "single bystanders" = "#4F4FC3",
                                    "unidentified" = "grey",
				    "other" = "black"))


# Save the plot 
ggsave(plot = g, filename = gsub("_adjmatrix.csv", "_clonogram.pdf", args[1]), device = cairo_pdf, dpi = 800)



