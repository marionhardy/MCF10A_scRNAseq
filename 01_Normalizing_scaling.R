library(Seurat)
# library(ggplot2)
library(tidyverse)
# BiocManager::install('glmGamPoi')
# remotes::install_github("NMikolajewicz/scMiko")
# BiocManager::install(c("Seurat","gridExtra","viridis","Polychrome","circlize"))

# The top three lines of this file are header lines. The third line contains 
# the total number of rows in all the three files in this folder 
# (genes.tsv, barcodes.tsv, matrix.mtx).
# The next lines (line number 4 onwards) have three columns:
# - The first column refers to the "gene id" index.
# - The second column refers to "cell id" index.
# - The third column represents the total UMI count per cell and gene combination.
# The 'gene id' and 'cell id' indices correspond to the entries in the 
# barcodes.tsv and genes.tsv files. The index in the MEX file is 1-based.

temp = read.csv("./data/all_genes.csv", na.strings=c("","NA"))
write.csv(temp, "./data/all_genes_filled.csv")

# One way to read it
counts =
  Seurat::ReadMtx("./data/DGE.mtx", 
                  "./data/cell_metadata.csv",
                  "./data/all_genes_filled.csv", 
                  cell.column = 1, 
                  feature.column = 2,
  cell.sep = ",", feature.sep = ",", skip.cell = 1, skip.feature = 1,
  mtx.transpose = TRUE,
  unique.features = TRUE,
  strip.suffix = FALSE)


seurat <- CreateSeuratObject(counts = counts, min.cells=5) # min cells = 1 we are conservative with filtering
rm(counts)
gc()

# seurat_f = subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
# rm(seurat)
# 
# seurat_sct <- SCTransform(seurat_f,verbose = TRUE)
# 
# ScaledCounts = seurat@assays$SCT@data

write_rds(seurat, "./data_output/seurat_raw.rds")
# write_rds(seurat_f, "./data_output/seurat_filtered.rds")
# write_rds(seurat_sct, "./data_output/seurat_filtered_scaled.rds")
# write.table(ScaledCounts, "./data_output/ScaledCounts.mtx")













