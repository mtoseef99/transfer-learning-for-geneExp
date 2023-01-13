library(Seurat)
library(SeuratDisk)
library(ProjecTILs)
library(ggplot2)
library(gridExtra)
library(reshape)
library(cowplot)
library(tidyverse)
library(RColorBrewer)
library(viridis)


# from projectR_ICI study in the TL review paper
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575


################################################################################################
# Steps in 
# 1- load the seurat object
# 2- load data matrix file to get sample column for meta data processing
# 3- load metadata and remove matrix from step 2
# 4- create a new seurat object to be sure
# 5- run standard processes such as normalization (not here because it is TPM), PCA, UMAP
# 6- Create and run TIL anchors and predictions
################################################################################################

# loading the save Seurat object
sf.seurat.obj <- readRDS("SadeFeldman.seurat.rds")




## Load expression matrix and metadata
# ==========================================================================================
# for windows
sf.dat <- read.delim("D:/RProjects/Figures_and_graphs/data/GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz", header = F, sep = "\t")

# for macOS
# sf.dat <- read.delim(sprintf("/Users/mtoseef/Data/TLRpaper/GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"), header = F, sep = "\t")


## data manipulation to save in desired count matrixformat
# ==========================================================================================

# genes <- sf.dat[c(-1, -2), 1]
# cells <- as.vector(t(sf.dat[1, 2:16292]))
# samples <- as.factor(t(sf.dat[2, 2:16292]))
# #
# sf.dat.modified <- sf.dat[c(-1, -2), 2:16292]
# colnames(sf.dat.modified) <- cells
# rownames(sf.dat.modified) <- genes


## ======= Create Seurat object
# sf.seurat.obj1 <- CreateSeuratObject(counts = sf.dat.modified, project = "SFDataset1", min.cells = 10)
# rm(sf.dat)
# rm(sf.dat.modified)


## ======= Adding meta data
sf.metadata <- read.delim('D:/RProjects/Figures_and_graphs/data/GSE120575/GSE120575_patient_ID_single_cells.txt.gz',
                       header = T, sep = "\t", skip = 19, nrows = 16291)

sf.metadata <- sf.metadata[, 1:7]


treat <- factor(ifelse(grepl("Post", samples), "Post", "Pre"))
response <- factor(sf.metadata$characteristics..response)
therapy <- factor(sf.metadata$characteristics..therapy)


sf.seurat.obj@meta.data$Sample <- samples
sf.seurat.obj@meta.data$Time <- treat
sf.seurat.obj@meta.data$Response <- response
sf.seurat.obj@meta.data$Therapy <- therapy

# saveRDS(sf.seurat.obj, file = "SadeFeldman.seurat.rds")
sf.seurat.backup <- sf.seurat.obj

###################################################
#
# Run standard steps

# seurat.integrated <- NormalizeData(seurat.integrated) # not here
seurat.integrated <- FindVariableFeatures(sf.seurat.backup, selection.method = "vst", nfeatures=2000)
# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)

ElbowPlot(seurat.integrated, ndims=50)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:18)

View(seurat.integrated@meta.data)


## ======= start - Mapping human peripheral blood cells
# reference_TIL <- LoadH5Seurat("./pbmc_multimodal.h5seurat")
# # 
# # # sctransform based normalization
# sf.TIL <- SCTransform(sf.seurat.obj, verbose = FALSE)
# # 
# # 
# anchors <- FindTransferAnchors(
#   reference = reference_TIL,
#   query = sf.TIL,
#   normalization.method = "SCT",
#   reference.reduction = "spca",
#   dims = 1:20
# )
# # 
# sf.TIL <- MapQuery(
#   anchorset = anchors,
#   query = sf.TIL,
#   reference = reference_TIL,
#   refdata = list(
#     celltype.l1 = "celltype.l1",
#     celltype.l2 = "celltype.l2",
#     predicted_ADT = "ADT"
#   ),
#   reference.reduction = "spca",
#   reduction.model = "wnn.umap"
# )


###################################################################
# define color palette
# Save as variable to global environment

# my_cols = brewer.pal(8, "Dark2")
p11_cols = c("purple", "orange")
p12_cols = c("orangered3", "yellow2", "mediumvioletred")
p13_cols = brewer.pal(8, "Dark2")

################################################################################################

p11 <- DimPlot(sf.TIL, reduction = 'ref.umap', group.by = "Response", cols = alpha(p13_cols, 0.4), pt.size = 1) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
        legend.position = "top") + 
  scale_x_continuous(name="UMAP 1") +
  scale_y_continuous(name="UMAP 2") 

p12 <- DimPlot(sf.TIL, reduction = 'ref.umap', group.by = "Therapy", cols = alpha(p13_cols, 0.5), pt.size = 1) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
        legend.position = "top") +
  scale_x_continuous(name="UMAP 1") +
  scale_y_continuous(name="UMAP 2")

p11 | p12

p13 <- DimPlot(sf.TIL, reduction = 'ref.umap', group.by = 'predicted.celltype.l1' ,  cols = alpha(p13_cols, 0.4), pt.size = 1)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
        legend.position = "top")+ 
  scale_x_continuous(name="UMAP 1") +
  scale_y_continuous(name="UMAP 2") +
  ggtitle("Predicted Cell Type based \n on ProjecTILs") 


# seurat.integrated.markers <- FindAllMarkers(seurat.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# seurat.integrated.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10

mapal <- colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(256)
p14 <- DoHeatmap(seurat.integrated,
                features = VariableFeatures(seurat.integrated)[1:20],
                # features = top10$gene,
                group.by= "ident",
                label = FALSE,
                group.bar.height = 0.0002,
                angle = 4,
                group.colors = c("red", "blue")) + NoLegend() 

cowplot_sf.dat <- cowplot::plot_grid(p11, p12, p13, p14, rel_widths = c(1, 1, 1, 1), ncol = 2, nrow = 2, labels = c('(A)', '(B)', '(C)', '(D)')) # + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.3))

cowplot_sf.dat


################################################################################################
################################################################################################


# # Idents(sf.TIL) <- 'predicted.celltype.l2'
# # VlnPlot(sf.TIL, features = c("CD4", "LILRA4"), sort = TRUE) + NoLegend()
# 
# treg_markers <- FindMarkers(sf.TIL, ident.1 = "Treg", only.pos = TRUE, logfc.threshold = 0.1)
# print(head(treg_markers))
# 
# 
# clusters <- DimPlot(sf.TIL, reduction = 'ref.umap', group.by = 'predicted.celltype.l1', label = TRUE) +
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
# condition.Response <- DimPlot(sf.TIL, reduction = 'ref.umap', group.by = 'Response') +
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
# condition.Therapy <- DimPlot(sf.TIL, reduction = 'ref.umap', group.by = 'Therapy') +
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
# p142 <- DoHeatmap(sf.TIL,
#                  features = VariableFeatures(sf.TIL)[1:100],
#                  group.by= "ident",
#                  label = FALSE,
#                  group.bar.height = 0.0002,
#                  angle = 4) + NoLegend()

# condition.Therapy|condition.Response|clusters|p14

# cowplot_sf.dat.til <- cowplot::plot_grid(clusters, condition.Response, condition.Therapy, p14, rel_widths = c(1, 1, 1, 1), ncol = 2, nrow = 2, labels = c('(A)', '(B)', '(C)', '(D)')) 
# 
# cowplot_sf.dat.til



## ======= end - Mapping human peripheral blood cells
# Diffuse intrinsic pontine glioma (DIPG)4 --> 1
# ANCA-associated vasculitis 4 --> 1
# Retina Development 2 --> 1

################################################################################################



# query.projected <- make.projection(sf.seurat.obj1, ref = reference_TIL)
# which.types <- table(query.projected$functional.cluster) > 20
# 
# stateColors_func <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc",
#                       "#FF0000", "#87f6a5", "#e812dd")
# states_all <- levels(ref$functional.cluster)
# names(stateColors_func) <- states_all
# cols_use <- stateColors_func[names(which.types)][which.types]
# 
# # Responder vs non Responder
# query.projected$functional.cluster <- factor(query.projected$functional.cluster,
#                                              levels = states_all)
# query.list <- SplitObject(query.projected, split.by = "Response")
# 
# norm.c <- table(query.list[["Non-responder"]]$functional.cluster)/sum(table(query.list[["Non-responder"]]$functional.cluster))
# norm.q <- table(query.list[["Responder"]]$functional.cluster)/sum(table(query.list[["Responder"]]$functional.cluster))
# 
# foldchange <- norm.q[which.types]/norm.c[which.types]
# foldchange <- sort(foldchange, decreasing = T)
# 
# tb.m <- melt(foldchange)
# colnames(tb.m) <- c("Cell_state", "Fold_change")
# pll <- list()
# p14 <- ggplot(tb.m, aes(x = Cell_state, y = Fold_change, fill = Cell_state)) + geom_bar(stat = "identity") +
#   scale_fill_manual(values = cols_use) + geom_hline(yintercept = 1) + scale_y_continuous(trans = "log2") +
#   theme(axis.text.x = element_blank(), legend.position = "left") + ggtitle("Responder vs. Non-responder")

################################################################################################
################################################################################################
