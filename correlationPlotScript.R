options(bitmapType='cairo-png') # this line of code needed BEFORE knitr code in chunk below
suppressPackageStartupMessages({
  library(Seurat)
  library(glmGamPoi)
  library(magrittr)
  library(RColorBrewer)
  library(Matrix)
  library(patchwork)
  library(DOSE)
  library(heatmap3)
  library(reshape2)
  library(Matrix)
  library(qlcMatrix)
  library(tidyverse)
})
# Both heatmaps have been saved to either outputs50K/DAOY-CBO/figures or outputs50K/ONS76-CBO/figures. Example path used: well/ludwig/users/ikb229/carp_project/analysis50K/outputs50K/DAOY-CBO/figures
# correlation plot for DAOY-CBO

# load the datasets; cbo first
cbo.annotate <- 
  readRDS("outputs50K/CBO/cbo.filtRiboMito.CellAnnot.SCT.CCdiffRegressed.clustered.clusterMap.rds")

# load the phase.daoyCbo annotated object
dyCbo.annotate <- 
  readRDS("outputs50K/DAOY-CBO/DaoyCbo.CellAnnot.filtRiboMito.SCT.CCdiffRegressed.clustered.clusterMap.rds")

av.exp.cbo <- sapply(sort(unique(cbo.annotate$celltype)), function(ct) rowMeans(cbo.annotate$RNA@data[,which(cbo.annotate$celltype == ct)] ))

avg.exp.dyCbo <- sapply(levels(dyCbo.annotate@active.ident), function(ct) rowMeans(dyCbo.annotate$RNA@data[,which(dyCbo.annotate@active.ident == ct)]))


genes2cor <- intersect(VariableFeatures(cbo.annotate), rownames(dyCbo.annotate))
corr2cbo.cl <- cor(avg.exp.dyCbo[genes2cor,], av.exp.cbo[genes2cor,], method = "spearman")
corr2cbo.cl

graphics.off()
par("mar")
par(mar=c(1,1,1,1))
dyCbo.cbo.correl.heatmap <- heatmap3(corr2cbo.cl, scale="none", margins=c(15,17),
         labRow = colnames(avg.exp.dyCbo), labCol = colnames(av.exp.cbo), cexRow=0.8, cexCol=0.8,
         col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))
png("./outputs50K/DAOY-CBO/figures/correlnPlot.png",
    width = 480, height = 480, units = "px", bg = "white", res = 300) # this is not saving, not sure why!
dev.off()


# correlation plot for ONS76-CBO

# load the stored annotated ons76-cbo seurat object
o76cbo.annotate <- readRDS("outputs50K/ONS76-CBO/ons76cbo.CellAnnot.filtRiboMito.SCT.CCdiffRegressed.clustered.clusterMap.rds")


av.exp.cbo <- sapply(sort(unique(cbo.annotate$celltype)), function(ct) rowMeans(cbo.annotate$RNA@data[,which(cbo.annotate$celltype == ct)] ))

avg.exp.ons76Cbo <- sapply(levels(o76cbo.annotate@active.ident), function(ct) rowMeans(o76cbo.annotate$RNA@data[,which(o76cbo.annotate@active.ident == ct)]))

genes2cor <- intersect(VariableFeatures(cbo.annotate), 
                       rownames(o76cbo.annotate))
corr2ons76.cbo.cl <- cor(avg.exp.ons76Cbo[genes2cor,], av.exp.cbo[genes2cor,], method = "spearman")

corr2ons76.cbo.cl

png("./outputs50K/ONS76-CBO/figures/correlnPlot.png",
    width = 480, height = 480, units = "px", bg = "white", res = 300) # this is not saving, not sure why!
graphics.off()
par("mar")
par(mar=c(1,1,1,1))
heatmap3(corr2ons76.cbo.cl, scale="none", margins=c(15,17),
         labRow = colnames(avg.exp.ons76Cbo), labCol = colnames(av.exp.cbo), cexRow=0.8, cexCol=0.8,
         col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))
dev.off()
