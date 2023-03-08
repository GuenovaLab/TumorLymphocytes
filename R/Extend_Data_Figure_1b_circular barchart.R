# library
library(tidyverse)

Data <- read.csv("clone_FACS_vs_scRNAseq/01_WaG_SE_clone_percentage.csv", header = TRUE, row.names = NULL)





# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(Data$group), ncol(Data)) )
colnames(to_add) <- colnames(Data)
to_add$group <- rep(levels(Data$group), each=empty_bar)
data <- rbind(Data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id)) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# Make the plot
ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+5, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = -10, xend = end, yend = -10), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,0), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)


p
