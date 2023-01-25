library("ggplot2")

coordinates = data.frame(x.min = rep(1, 5),
			 x.max = rep(1, 5) + 0.5,
			 y.min = seq(1, 5),
			 y.max = seq(1, 5) + 0.5,
			 clonal.group = factor(c("main clone",
                                   "related to main clone", 
                                   "bystander groups", 
                                   "single bystanders", 
                                   "unidentified"), levels = c("main clone", "related to main clone", "bystander groups", "single bystanders", "unidentified")),
			 fill = c("#EC4C4C", "#FDC37E", "#D29FF8", "#4F4FC3", "grey"))

gg = ggplot(coordinates, aes(xmin = x.min, xmax = x.max, ymin = y.min, ymax = y.max, fill = clonal.group)) + geom_rect() + 
     scale_fill_manual(values = c("main clone" = "#EC4C4C", 
                                    "related to main clone" = "#FDC37E",
                                    "bystander groups" = "#D29FF8",
                                    "single bystanders" = "#4F4FC3",
                                    "unidentified" = "grey")) +
     guides(size = FALSE, guide_legend(override.aes = list(size=12, stroke = 12))) +
     theme_void()

ggsave(plot = gg, filename = "legend.pdf", width = 4, height = 4)
