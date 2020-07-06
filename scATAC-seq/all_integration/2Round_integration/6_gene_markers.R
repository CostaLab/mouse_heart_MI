mute <- suppressPackageStartupMessages
mute(library(ggplot2))
mute(library(Seurat))
mute(library(stringr))
mute(library(magrittr))
mute(library(readr))
mute(library(Matrix))
mute(library(tidyr))
mute(library(dplyr))
mute(library(plotly))
mute(library(Signac))
mute(library(Seurat))
mute(library(GenomeInfoDb))
mute(library(EnsDb.Mmusculus.v79))
mute(library(cluster))
mute(library(clustree))
mute(library(mclust))
mute(library(cowplot))
mute(library(gridExtra))
mute(library(future.apply))

plan("multiprocess", workers = 100)
options(future.globals.maxSize = 100000* 1024^2)




obj.integrated <- readRDS(file = "obj.integrated.Rds")

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

