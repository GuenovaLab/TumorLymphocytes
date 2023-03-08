install.packages("gprofiler2")
library(gprofiler2)

gostres <- gost(query = c("ABCC6", "ACTG1",
                          "ANXA2",
                          "B2M",
                          "BEST1",
                          "C4orf51",
                          "CCT6A",
                          "CD96",
                          "DDX5",
                          "DNAJB6",
                          "EEF1A1",
                          "EIF4A2",
                          "ESYT2",
                          "FTH1",
                          "FYCO1",
                          "HLA-A",
                          "HLA-E",
                          "HNRNPH1",
                          "HNRNPR",
                          "HSPA8",
                          "IFITM1",
                          "IFITM2",
                          "KLF6",
                          "LCP1",
                          "LGALS8",
                          "MRNIP",
                          "MT-CO1",
                          "MT-RNR1",
                          "MT-RNR2",
                          "MYL6",
                          "PGK1",
                          "PLIN2",
                          "RBM5",
                          "RGS1",
                          "RNF213",
                          "S100A11",
                          "SAT1",
                          "SERF2",
                          "SLC2A3",
                          "SQSTM1",
                          "TMSB4X",
                          "TOMM20",
                          "TRBC2",
                          "TSPAN13",
                          "UBC",
                          "VDAC3"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

names(gostres)

head(gostres$result)


gostplot(gostres, capped = TRUE, interactive = TRUE)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p


pp <- publish_gostplot(p, highlight_terms = c("GO:0042612", "GO:0002480", "GO:0002716", "KEGG:04612"), 
                       width = NA, height = NA, filename = NULL )
pp


publish_gosttable(gostres, highlight_terms = gostres$result[c(1:2,10,100:102,120,124,125),],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)

