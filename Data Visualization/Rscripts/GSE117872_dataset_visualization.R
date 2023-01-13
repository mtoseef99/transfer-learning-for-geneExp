library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(ComplexHeatmap)


## Load expression matrix and metadata
# =========================================
# for windows
ma.dat2 <- read.delim("D:/RProjects/Figures_and_graphs/data/GSE117872/GSE117872_good_Data_TPM.txt.gz", header = T, sep = "\t")

# metadata
ma.metadata <- read.delim("D:/RProjects/Figures_and_graphs/data/GSE117872/GSE117872_good_Data_cellinfo.txt.gz", header = T, sep = "\t")

# creating a Seurat object from for this data
# =========================================
ma.seurat.obj <- CreateSeuratObject(counts = ma.dat2, project = "Dataset2", min.cells = 10)

status <- factor(ma.metadata$drug_status)
celltype <- factor(ma.metadata$cluster)
cellcolor <- factor(ma.metadata$cell_color)
origin <- factor(ma.metadata$origin)

ma.seurat.obj@meta.data$'Drug Status' <- status
ma.seurat.obj@meta.data$'Cell Type' <- celltype
ma.seurat.obj@meta.data$'Original Sample' <- cellcolor
ma.seurat.obj@meta.data$'Cell Origin' <- origin

# saving seurat without any data pre-processing
# saveRDS(ma.seurat.obj, file = "GSE117872.seurat.rds")

seurat.integrated <- ma.seurat.obj

seurat.integrated <- NormalizeData(seurat.integrated)
seurat.integrated <- FindVariableFeatures(seurat.integrated, selection.method = "vst", nfeatures=2000)



# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)

ElbowPlot(seurat.integrated, ndims=50)

seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)

# umap_tx = seurat.integrated@reductions$umap@cell.embeddings %>% 
#   as.data.frame() %>% cbind(tx = seurat.integrated@meta.data$`Drug Status`)
# 
# ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, color=tx)) + geom_point()





p21 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = "Original Sample")
p22 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = "Drug Status") # it is cancer type and not drug status
p23 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = "Cell Type",
              cols = c('red','green','blue'))

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
p24 <- DoHeatmap(seurat.integrated,
                features = VariableFeatures(seurat.integrated)[1:500],
                group.by= "ident",
                label = FALSE,
                group.bar.height = 0.0002,
                angle = 4) + NoLegend()


grid.arrange(p21, p22, p23, p24, ncol = 2, nrow = 2)

cowplot_ma.dat <- cowplot::plot_grid(p21, p22, p23, p24, rel_widths = c(1, 1, 1, 1), rel_heights = c(0.6, 0.6, 0.6), ncol = 2, nrow = 2, labels = c('(A)', '(B)', '(C)', '(D)')) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.3))

cowplot_ma.dat













# mad2 <- SCTransform(seurat.integrated, verbose = FALSE)
# anchors <- FindTransferAnchors(
#   reference = reference,
#   query = seurat.integrated,
#   normalization.method = "SCT",
#   reference.reduction = "spca",
#   dims = 1:30
# )
# 
# mad2 <- MapQuery(
#   anchorset = anchors,
#   query = seurat.integrated,
#   reference = reference,
#   refdata = list(
#     celltype.l1 = "celltype.l1",
#     celltype.l2 = "celltype.l2",
#     predicted_ADT = "ADT"
#   ),
#   reference.reduction = "spca",
#   reduction.model = "wnn.umap"
# )
# 
# 
# p11 = DimPlot(mad2, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE)
# p12 = DimPlot(mad2, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE)
# p11






















