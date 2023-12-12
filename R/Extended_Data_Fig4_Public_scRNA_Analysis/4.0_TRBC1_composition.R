library(Seurat)
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Pacome/Results/scRNA_YunTsan_Revision/output/")
Seu =  qs::qread("merged_seurat.qs")

DimPlot(Seu, group.by = "celltype", cells = sample(colnames(Seu), 100000))

Seu. = Seu[,Seu$study == "GSE173205__GSE165623_MF_TCR"]
table(Seu.$study, Seu.$sample_id)
table(Seu.$sample_id, Seu.$has_TCR)
table(Seu.$sample_id, Seu.$is_clonal)

library(ggpubr)
png("TRBC1_MF_vs_Healthy.png")
VlnPlot(Seu., features = "TRBC1", group.by = "condition") + 
  stat_compare_means(aes( label = after_stat(p.signif)),
                     method = "t.test", ref.group = "MF")
dev.off()

pdf("TRBC1_Clonal_vs_Non_Clonal.pdf")
meta = Seu.@meta.data
# Seu. = SCTransform(Seu.)
meta$TRBC1 = (Seu.@assays$RNA@counts["TRBC1",] > 0)
meta %>% 
  dplyr::filter(celltype != "T_helper") %>%
  dplyr::filter(condition != "HD") %>%
   group_by(sample_id, is_clonal) %>% 
  dplyr::summarise(percent_TRBC1 = 100 *sum(TRBC1) / n()) %>% 
  ggplot(aes(x = sample_id, y = percent_TRBC1, fill = is_clonal)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#d1d1d1ff", "#681fc2")) +
  theme_classic() + xlab("") + 
  theme(axis.text.x = element_text(angle = 90)) 
dev.off()