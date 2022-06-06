library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(cluster)
library(clustree)
library(mclust)
library(cowplot)
library(gridExtra)



obj.integrated <- readRDS(file = "obj.umap_clusters.Rds")

Idents(obj.integrated) <- "seurat_clusters"
markers.all <- FindAllMarkers(object = obj.integrated,
                              min.pct = 0.1,
                              logfc.threshold = 0.1,
                              test.use = "LR",
                              only.pos = TRUE)

if(length(markers.all) > 0){
    markers.all <- markers.all %>% group_by(cluster) %>%
            arrange(desc(avg_logFC), .by_group = TRUE)

    markers <- split(markers.all, markers.all$cluster)
    WriteXLS::WriteXLS(lapply(markers, subset,
                              subset = p_val_adj < 0.05),
                           ExcelFileName = "markers.xlsx",
                           SheetNames = names(markers))
}

