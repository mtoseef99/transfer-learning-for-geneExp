# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(patchwork)


# get the data from the following link and download all samples, place them under four differenet folders, such as:
# GSM3371684_Control
# GSM3371685_aPD1
# GSM3371686_aCTLA4
# GSM3371687_aPD1_aCTLA4

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119352

# get data location
dirs <- list.dirs(path = 'data/GSE119352/', recursive = F, full.names = F)

# reading all .gz files into cts (countmatrix) and then creating SeuratObject for each data sample
for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('data/GSE119352/',x,'/matrix.mtx.gz'),
                 features = paste0('data/GSE119352/',x,'/genes.tsv.gz'),
                 cells = paste0('data/GSE119352/',x,'/barcodes.tsv.gz'))
  ReadMtx
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}


# merge datasets
merged_seurat <- merge(GSM3371684_Control, y = c(GSM3371685_aPD1, GSM3371686_aCTLA4, GSM3371687_aPD1_aCTLA4),
                       add.cell.ids = ls()[32:35], # this is to give index number using ls() 
                       project = 'Dataset12')


# QC & filtering -----------------------
View(merged_seurat@meta.data)
#create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# replace the columns with aPD1_aCTLA4 to aPD1-aCTLA4
merged_seurat$sample <- gsub("aPD1_aCTLA4", "aPD1-aCTLA4", merged_seurat$sample)


# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Therapy', 'Barcode'), 
                                    sep = '_')


# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

######################################################
# saveRDS(merged_seurat, file = "GSE119352.seurat.rds")
######################################################




# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 10)

merged_seurat_filtered
merged_seurat




# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
# to select dimensions for PCA 
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:10)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:10)



# plot
p31 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient') +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
p32 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Therapy')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
p33 <- DoHeatmap(merged_seurat_filtered,
                 features = VariableFeatures(merged_seurat_filtered)[1:50],
                 group.by= "ident",
                 label = FALSE,
                 group.bar.height = 0.02,
                 angle = 4) + NoLegend()

# 
# grid.arrange(p31, p32, p33, ncol = 3, nrow = 1)
# 
# cowplot_dat3 <- cowplot::plot_grid(p31, p32, p33, NULL, rel_widths = c(1, 1, 1, 0), ncol = 2, nrow = 2, labels = c('(A)', '(B)', '(C)')) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.3))
# 
# cowplot_dat3


patchwork.plots <- p33 | (p31 / p32)
patchwork.plots + plot_annotation(tag_levels = c('A', '1'), tag_prefix = '(', tag_suffix = ')') +
  plot_layout(widths = c(1.2, 1.8))












