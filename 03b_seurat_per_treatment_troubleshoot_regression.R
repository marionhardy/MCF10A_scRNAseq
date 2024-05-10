library(Seurat)
library(tidyverse)
library(Matrix)
library(SeuratDisk)


seurat_sct <- readRDS("./data_output/per_treatment/IM_Oligomycin_seurat_scaled.rds")

Idents(seurat_sct) = seurat_sct@meta.data$replicate

DimPlot(seurat_sct, reduction = "umap")+
  labs(title = "n neighbors = 10")

seurat_sct@meta.data <-
  seurat_sct@meta.data %>%
  mutate(batch =
           ifelse(replicate_name == paste0(condition), 1, 2))


for (i in 1:length(samples)){
  condition = samples[i]
  subseurat = subset(seurat, subset = treatment == condition)
  
  # scaling per condition
  
  seurat_sct <- SCTransform(subseurat, vars.to.regress = "replicate")
  write_rds(subseurat, paste0("./data_output/per_treatment/",condition,"_seurat_raw.rds"))
  
  print(paste0("Scaling done for condition ",condition," which is ",i,"/32"))
  
  # QC plots
  Idents(seurat_sct) = seurat_sct@meta.data$treatment
  p1 = VlnPlot(subseurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), 
          ncol = 3, pt.size = 0.0001)
  ggsave(plot = p1, paste0("./figures/per_treatment/",condition,"_violin_features.png"))
  p2 = FeatureScatter(subseurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    theme(legend.position="none")
  ggsave(plot =p2, paste0("./figures/per_treatment/",condition,"_scatter_features.png"))

  # pca 
  
  Idents(seurat_sct) = seurat_sct@meta.data$treatment
  seurat_sct <- RunPCA(seurat_sct, features = VariableFeatures(object = seurat_sct))
  
  p3 = DimPlot(seurat_sct, reduction = "pca")+
    labs(title = "All data filtered and scaled")
  ggsave(plot =p3, paste0("./figures/per_treatment/",condition,"_PCA.png"))
  
  # umap
  
  seurat_sct <- FindNeighbors(seurat_sct, dims = 1:40)
  seurat_sct <- FindClusters(seurat_sct, resolution = 0.5)
  
  seurat_sct <- RunUMAP(seurat_sct, dims = 1:40, seed.use = 054057, n.neighbors = 5)
  p4 = DimPlot(seurat_sct, reduction = "umap")+
    labs(title = "n neighbors = 5")
  ggsave(plot =p4, paste0("./figures/per_treatment/",condition,"_UMAP_nn5.png"))
  
  seurat_sct <- RunUMAP(seurat_sct, dims = 1:40, seed.use = 054057, n.neighbors = 10)
  p5 = DimPlot(seurat_sct, reduction = "umap")+
    labs(title = "n neighbors = 10")
  ggsave(plot =p5, paste0("./figures/per_treatment/",condition,"_UMAP_nn10.png"))
  
  print(paste0("UMAP done for condition ",condition," which is ",i,"/32"))
  
  write_rds(seurat_sct, paste0("./data_output/per_treatment/",condition,"_seurat_scaled.rds"))
  
  print(paste0("Saving final rds object of ",condition," which is ",i,"/32"))
  
  seurat_sct[["RNA"]] <- as(object = seurat_sct[["RNA"]], Class = "Assay")
  seurat_sct[["SCT"]] <- as(object = seurat_sct[["SCT"]], Class = "Assay")
  
  SaveH5Seurat(seurat_sct, filename = paste0("./data_output/per_treatment/",
                                             condition,"_seurat_sct.h5Seurat"), 
               overwrite = T)
  Convert(paste0("./data_output/per_treatment/",
                 condition,"_seurat_sct.h5Seurat"), dest = "h5ad", overwrite = TRUE)
  
  rm(seurat_sct)
  
}











