library(Seurat)
atac <- readRDS(file="save_res/trans_add_tdtomato_annotation_res1.0.Rds")

Idents(atac) <- "integrated_snn_res.0.1"
keep = setdiff(unique(atac$`integrated_snn_res.0.1`), c("4", "5", "6"))
atac <- subset(atac, idents=keep)

saveRDS(atac, file="save_res/trans_add_tdtomato_annotation_res1.0_remove_res0.1_456.Rds")





