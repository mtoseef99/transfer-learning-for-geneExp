library(affy)
library(GEOquery)
library(tidyverse)
library(ggplot2)

# https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE6008
# reading in .cel files
raw.data.GSE6008 <- ReadAffy(celfile.path = "data/GSE6008/")


# ======================================================================================================

# visualizing some initial plots
hist(raw.data.GSE6008)
boxplot(raw.data.GSE6008)

raw.exp <- exprs(raw.data.GSE6008)

# MA plot before normalization
ma.plot(rowMeans(log2(raw.exp)), log2(raw.exp[,2]) - log2(raw.exp[,3]), cex = 1)

# MA plot after normalization
ma.plot(rowMeans(log2(rma.GSE6008)), log2(rma.GSE6008[,2]) - log2(rma.GSE6008[,3]), cex = 1)



# ======================================================================================================

# getting series matrix file and creating metadata from index 1 of GSE6008
gseGSE6008 <- getGEO("GSE6008")
meta.GSE6008 <- pData(gseGSE6008[[1]])


#rename cols names 
meta.GSE6008.mod <- meta.GSE6008 %>%
  rename(cell_status = source_name_ch1) %>%
  rename(tumor = characteristics_ch1) %>%
  mutate(tumor = gsub("Tumor_Type: ", "", tumor))

statuscol <- as.numeric(factor(meta.GSE6008.mod$tumor)) + 1

View(meta.GSE6008.mod)

# ======================================================================================================


# creating boxp;lot before normalization
# graphics.off()
par(mar=c(9,3,2,0.2))
boxplot(raw.data.GSE6008, xaxt = "n", col = statuscol, main = "GSE6008 - Before Normalization color with Cancer types")
mtext("Samples", side = 1, line = 1, cex = 1)
mtext("Expression Values", side = 2, line = 2, cex = 1)


# Plot the chart.
boxplot(raw.data.GSE6008,
        main = "GSE6008 - Before Normalization",
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
rma.exp.GSE6008 <- rma(raw.data.GSE6008)
## These steps will be performed
# Background correcting
# Normalizing
# Calculating Expression



# get the expression matrix for all samples and genes after Background correcting ,Normalizing, and Calculating Expression
rma.GSE6008 <- exprs(rma.exp.GSE6008)


#writing and saving expression matrix to a .txt file
write.table(rma.GSE6008, file = "data/GSE6008/rma_GSE6008.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# ======================================================================================================


# creating boxplot after visualization
# graphics.off()
# par(mar=c(9,3,2,0.2))
# boxplot(rma.GSE6008, xaxt = "n", col = statuscol, main = "GSE6008 - After Normalization")
# mtext("Samples", side = 1, line = 1, cex = 1)
# mtext("Expression Values", side = 2, line = 2, cex = 1)


# Plot the chart.
boxplot(rma.GSE6008,
        main = "GSE6008 - After Normalization",
        xlab = "Samples",
        xaxt = "n",
        ylab = "Expression Values",
        col = statuscol,
        border = "black"
)
# Add a legend
legend("bottomright", legend = c("Clear-cell ovarian carcinoma", "Endometrioid carcinoma", "Mucinous carcinoma", "Serous carcinoma", "Normal") , 
       col = c(rgb(223, 83, 107, maxColorValue = 255) , rgb(97, 208, 79, maxColorValue = 255), rgb(34, 151, 230, maxColorValue = 255), rgb(205, 11, 188, maxColorValue = 255), rgb(40, 226, 229, maxColorValue = 255)),
       # col = statuscol,
       pch=20 , pt.cex = 1.5, cex = 0.5, horiz = FALSE, 
       inset = c(0.00003, 0.000001)
       )




# ======================================================================================================
# Before batch correction
pcDAT <- prcomp(t(rma.GSE6008))
fviz_pca_ind(pcDAT, geom.ind = "point")

















