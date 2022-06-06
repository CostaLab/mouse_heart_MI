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
library(mclust)
library(cowplot)
library(gridExtra)

dir.in <- "../../Fabian_scRNA/current_final_RDS_19_04_2021/"

ref_rna <- readRDS(file=file.path(dir.in, "VSMC&Pericytes_new_added_sub_cell_type.rds"))
Idents(ref_rna) <- "orig.geno"
ref_rna <- subset(ref_rna, idents="Myh11")
ref_rna <- FindVariableFeatures(ref_rna)


atac <- readRDS(file="save_res/obj.integrated_res.Rds")


transfer.anchors <- FindTransferAnchors(reference = ref_rna, query = atac,
                     features = VariableFeatures(object = ref_rna),
                     reference.assay = "RNA", query.assay = "gene.activity",
                     reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = ref_rna$sub_cell_type,
                        weight.reduction = atac[["umap_pca_integrated"]], dims=1:2)

celltype.predictions2 <- TransferData(anchorset = transfer.anchors, refdata = ref_rna$sub_cell_type,
                        weight.reduction = atac[["pca_integrated"]], dims=1:30)


#atac <- AddMetaData(atac, metadata = celltype.predictions)

saveRDS(list(celltype.predictions, celltype.predictions), "save_res/label_transfer.predict.Rds")

atac$predict_umap_celltype <- celltype.predictions[colnames(atac), "predicted.id"]
atac$predict_pca_celltype <- celltype.predictions2[colnames(atac), "predicted.id"]
saveRDS(atac, "save_res/trans.Rds")


