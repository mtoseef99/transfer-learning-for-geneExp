library(affy)
library(GEOquery)
library(tidyverse)
library(ggplot2)

# TransComp-R paper by Brubaker
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16879
# reading in .cel files
raw.data.GSE16879 <- ReadAffy(celfile.path = "data/GSE16879/")


# ======================================================================================================

# visualizing some initial plots
hist(raw.data.GSE16879)
boxplot(raw.data.GSE16879)

raw.exp <- exprs(raw.data.GSE16879)
ma.plot.bf.nm <- ma.plot(rowMeans(log2(raw.exp)), log2(raw.exp[,2]) - log2(raw.exp[,3]), cex = 1)
ma.plot.af.nm <- ma.plot(rowMeans(log2(raw.data.GSE16879)), log2(raw.data.GSE16879[,2]) - log2(raw.data.GSE16879[,3]), cex = 1)



# ======================================================================================================
# getting series matrix file and creating metadata from index 1 of GSE6008
gse16879 <- getGEO("GSE16879")
meta.gse16879 <- pData(gse16879[[1]])

#rename cols names 
meta.GSE16879.mod <- meta.gse16879 %>%
  dplyr::rename(tissue = characteristics_ch1) %>%
  dplyr::rename(disease = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(disease = gsub("disease: ", "", disease))

statuscol_gse16879 <- as.numeric(factor(meta.GSE16879.mod$disease)) + 1

View(meta.GSE16879.mod)

# ======================================================================================================



# Plot the chart.
boxplot(raw.data.GSE16879,
        main = "GSE16879 - Before Normalization",
        xlab = "Samples",
        xaxt = "n",
        ylab = "Expression Values",
        col = statuscol,
        border = "black"
)
legend("topright", legend = c("Clear-cell ovarian carcinoma", "Endometrioid carcinoma", "Mucinous carcinoma", "Serous carcinoma", "Normal") , 
       col = c(rgb(223, 83, 107, maxColorValue = 255) , rgb(97, 208, 79, maxColorValue = 255), rgb(34, 151, 230, maxColorValue = 255), rgb(205, 11, 188, maxColorValue = 255), rgb(40, 226, 229, maxColorValue = 255)),
       # col = statuscol,
       pch=20 , pt.cex = 1.5, cex = 0.5, horiz = FALSE, 
       inset = c(0.00003, 0.000001)
)


# ======================================================================================================
# performing rma normalization
rma.exp.GSE16879 <- rma(raw.data.GSE16879)
## These steps will be performed
# Background correcting
# Normalizing
# Calculating Expression


# get the expression matrix for all samples and genes after Background correcting ,Normalizing, and Calculating Expression
rma.GSE16879 <- exprs(rma.exp.GSE16879)


#writing and saving expression matrix to a .txt file
write.table(rma.GSE16879, file = "data/GSE16879/rma.GSE16879.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# ======================================================================================================
# Plot the chart.
boxplot(rma.GSE16879,
        main = "GSE16879 - After Normalization",
        xlab = "Samples",
        xaxt = "n",
        ylab = "Expression Values",
        col = statuscol,
        border = "black"
)




# ======================================================================================================
GSE16879.mas5 = mas5(raw.data.GSE16879)
## background correction: mas 
# PM/MM correction : mas 
# expression values: mas 
# background correcting...done.
# 54675 ids to be processed

exprSet.nologs = exprs(GSE16879.mas5)
colnames(exprSet.nologs)

colnames(raw.exp)
heatmap(exprSet.nologs)


