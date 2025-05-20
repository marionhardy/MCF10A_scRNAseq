
library(tidyverse)

tfheat = read.csv("./data/cellxgene.all_data_scaled_seurat.heatmapData_zscore.csv")

tfheat_long =
  pivot_longer(tfheat,
               names_to = "TF",
               values_to = "Z_score",
               c(2:14))


# Plot heatmap

tfheat_long %>% 
  ggplot(aes(x = TF,
             fill = Z_score,
             y = treatment))+
  theme_bw()+
  geom_raster()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",)+
  theme(axis.text.x = element_text(angle = 80, hjust=1, size = 8))+
  scale_y_discrete(limits = rev(levels(tfheat_long$treatment))) + # to get alphabetical order
  labs(title = "Transcription factors z-scores across conditions",
       subtitle = "scRNAseq 10A data")

# As a heatmap in order to cluster

library(pheatmap)
library(RColorBrewer)

# Compute means in long format
tfheat_long = 
  tfheat_long %>% 
  group_by(treatment, TF) %>% 
  mutate(mean_zscore = mean(Z_score)) 

tf_mean = distinct(tfheat_long, treatment, mean_zscore)

# Convert to wide for heatmap
# Make matrix

mats = 
  tf_mean %>% 
  pivot_wider(names_from = TF,
              values_from = mean_zscore)

mats = mats[,-1]
rownames(mats) = unique(tf_mean$treatment)

# annot = as.data.frame(annot$ResponseType)
# rownames(annot) = rownames(mats)
# colnames(annot) = "ResponseType"


svg(file="./figures/heatmap_TF_JGA.svg", width = 18, height = 8, pointsize = 7)

pheatmap(mats, main="Transcription factors mean(z-scores) across conditions", 
         color = rev(brewer.pal(9,"RdBu")),
         cluster_cols=T,
         fontsize_row=10, border_color=NA, 
         legend = T)
dev.off()









