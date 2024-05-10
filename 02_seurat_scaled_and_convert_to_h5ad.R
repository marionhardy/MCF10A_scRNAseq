library(Seurat)
library(tidyverse)
library(Matrix)
library(SeuratDisk)

seurat = readRDS("./data_output/seurat_mfiltered_scaled.rds")
rowdata = read.csv("./data/all_genes_filled.csv")
coldata = read.csv("./data/cell_metadata.csv")


# explore ----------------------------------------------------------------------

seurat@meta.data$replicate_name = coldata$sample
seurat@meta.data$treatment = coldata$group

# rename from id to gene symbols -----------------------------------------------
# in theory seurat v5 allows to do that with a function-------------------------

temp = data.frame(orig_gene = rownames(seurat))
temp2 = data.frame(updated = rowdata$gene_name,
                   orig_gene = rowdata$gene_id) %>%
  filter(!duplicated(updated))

# if temp2$updated is NA, fill with temp

temp3 = left_join(temp, temp2) # now temp3 is ordered like the genes in seurat

temp3 =
  temp3 %>%
  mutate(updated = coalesce(updated,orig_gene)) #works

seurat@assays$RNA@layers$counts@Dimnames[[1]] <- temp3$updated
# seurat@assays$RNA@layers$meta.features = rowdata
rownames(seurat@assays$RNA@features@.Data) = temp3$updated
rownames(seurat@assays$RNA@features) = temp3$updated
rownames(seurat)

# seurat@assays$SCT@counts@Dimnames[[1]] <- temp3$updated # has not been created yet
# seurat@assays$SCT@data@Dimnames[[1]] <- temp3$updated
# seurat@assays$SCT@meta.features = rowdata

# get the replicates into rowdata

treatmentlist = c("IM","IM_IL6","IM_Galloflavin", "IM_PD","IM_NoCT","IM_NoHC",
                  "IM_NoINS","IM_MK8722","IM_NoGluc","IM_UK5099","IM_Ipasertib",
                  "IM_Glut","IM_EGF","IM_Rapamycin","Full_GM","IM_Oligomycin")

seurat@meta.data <-
  seurat@meta.data %>%
  mutate(replicate=
           ifelse(replicate_name %in% treatmentlist, 1,2))

strrep = sub(pattern = "\\_2","", seurat@meta.data$replicate_name)
seurat@meta.data$treatment = strrep

head(seurat@meta.data)

# QC and filtering -------------------------------------------------------------

# How I would usually calculate mitochondrial percentage
## get the features that have mt in them

seurat$mitoPercent = PercentageFeatureSet(seurat, pattern='^MT-')
head(seurat@meta.data) # we now have the mitopercent added to the metadata

write_rds(seurat, "./data_output/seurat_raw_annotated.rds")

# Just like Cell Ranger output, feature in the following results represents gene. 
# nFeature_ is the number of genes detected in each cell. nCount_ is the 
# total number of molecules detected within a cell.

Idents(seurat) = seurat@meta.data$treatment
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), 
        ncol = 3, pt.size = 0.0001)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none") #  show correlations between features and gene counts

## Rerunning

seurat = subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & mitoPercent<20)

write_rds(seurat, "./data_output/seurat_raw_annotated_filtered.rds")

seurat_sct <- SCTransform(seurat, vars.to.regress = "replicate", 
                          verbose = TRUE)
rm(seurat)

# pca 

Idents(seurat_sct) = seurat@meta.data$treatment
seurat <- RunPCA(seurat_sct, features = VariableFeatures(object = seurat_sct))

DimPlot(seurat_sct, reduction = "pca")+
  labs(title = "All data filtered and scaled")

# umap

seurat_sct <- RunUMAP(seurat_sct, dims = 1:40, seed.use = 054057, n.neighbors = 5)
p1 = DimPlot(seurat_sct, reduction = "umap")+
  labs(title = "n neighbors = 5")

seurat_sct <- RunUMAP(seurat_sct, dims = 1:40, seed.use = 054057, n.neighbors = 10)
p2 = DimPlot(seurat_sct, reduction = "umap")+
  labs(title = "n neighbors = 10")

p1
p2

# Annotate

treatmentlist = c("IM","IM_IL6","IM_Galloflavin", "IM_PD","IM_NoCT","IM_NoHC",
                  "IM_NoINS","IM_MK8722","IM_NoGluc","IM_UK5099","IM_Ipasertib",
                  "IM_Glut","IM_EGF","IM_Rapamycin","Full_GM","IM_Oligomycin")

seurat_sct@meta.data <-
  seurat_sct@meta.data %>%
  mutate(replicate=
           ifelse(replicate_name %in% treatmentlist, 1,2))

strrep = sub(pattern = "\\_2","", seurat_sct@meta.data$replicate_name)
seurat_sct@meta.data$treatment = strrep

head(seurat_sct@meta.data)

# save 

write_rds(seurat_sct, "./data_output/seurat_mfiltered_scaled.rds")

# Convert to h5ad for cellxgene
# SeuratDisk

seurat_sct[["RNA"]] <- as(object = seurat_sct[["RNA"]], Class = "Assay")
seurat_sct[["SCT"]] <- as(object = seurat_sct[["SCT"]], Class = "Assay")

SaveH5Seurat(seurat_sct, filename = "./data_output/seurat_sct.h5Seurat", overwrite = T)
Convert("./data_output/seurat_sct.h5Seurat", dest = "h5ad", overwrite = TRUE)






