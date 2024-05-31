
library(tidyverse)

scenicz = read.csv("./data/albeck_pyscenic_analysis_v2_condition_specific_regulons_zscores.csv")

table(scenicz$Conditions)

# removing oscar's typo in the name of the condition

scenicz = 
  scenicz %>% 
  mutate(Conditions = case_when(Conditions == "base_i_maging_medium_im" ~ "base_imaging_medium", 
                                TRUE ~ Conditions))

scenicz$Conditions = factor(scenicz$Conditions)

scenicz %>% 
  ggplot(aes(x = regulon,
             fill = Zscores,
             y = Conditions))+
  theme_bw()+
  geom_raster()+
  scale_fill_gradient2(low = "blue",
                      mid = "white",
                      high = "red",)+
  theme(axis.text.x = element_text(angle = 80, hjust=1, size = 8))+
  scale_y_discrete(limits = rev(levels(scenicz$Conditions))) + # to get alphabetical order
  labs(title = "All significantly enriched regulons across conditions",
        subtitle = "scRNAseq data, Oscar Davalos pySCENIC v2")


ggsave(plot = last_plot(), "./figures/TF_all_regulons_zscore.svg", width = 35, height = 14,
       units = "cm",dpi = 300, scale = 1.3)


# let's try to plot it as a heatmap so I can cluster it



# Do the scaling to get z-scores
# z-score don't require the log tranformation to scale your conditions
# in a way where you can compare them

library(pheatmap)
library(RColorBrewer)

mats = 
  scenicz %>% 
  select(-X) %>% 
  pivot_wider(names_from = regulon,
                   values_from = Zscores) %>% 
  select(-Conditions) 

rownames(mats) = unique(scenicz$Conditions)

# annot = as.data.frame(annot$ResponseType)
# rownames(annot) = rownames(mats)
# colnames(annot) = "ResponseType"


svg(file="./figures/heatmap_TF_all_regulons.svg", width = 18, height = 8, pointsize = 7)

pheatmap(mats, main="All significantly enriched regulons across conditions", 
         subtitle = "scRNAseq data, Oscar Davalos pySCENIC v2",
         color = rev(brewer.pal(9,"RdBu")),
         cluster_cols=T,
         fontsize_row=10, border_color=NA, 
         legend = T)
dev.off()

# Remove the uninteresting regulons

roi = mats[,abs(colMeans(mats))>.001]
rownames(roi) = rownames(mats)

svg(file="./figures/heatmap_TF_regulons_top107.svg", width = 16, height = 8, pointsize = 7)

pheatmap(roi, main="Top 107 significantly enriched regulons across conditions", 
         subtitle = "scRNAseq data, Oscar Davalos pySCENIC v2",
         color = rev(brewer.pal(9,"RdBu")),
         cluster_cols=T,
         fontsize_row=10, border_color=NA, 
         legend = T)
dev.off()

# more stringent clustering

roi = mats[,abs(colMeans(mats))>.005]
rownames(roi) = rownames(mats)

svg(file="./figures/heatmap_TF_regulons_top64.svg", width = 14, height = 8, pointsize = 7)

pheatmap(roi, main="Top 64 significantly enriched regulons across conditions", 
         subtitle = "scRNAseq data, Oscar Davalos pySCENIC v2",
         color = rev(brewer.pal(9,"RdBu")),
         cluster_cols=T,
         fontsize_row=10, border_color=NA, 
         legend = T)
dev.off()


















