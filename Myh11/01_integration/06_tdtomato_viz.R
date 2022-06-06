library(ggplot2)
library(dplyr)
library(ggsci)
library(Seurat)


myh11 <- readRDS(file="save_res/trans_add_tdtomato_annotation_res1.0_remove_res0.1_456.Rds")

pdf("viz/myh11_totomato_without_imputation.pdf")
DefaultAssay(myh11) <- "gene.tdtomato"

FeaturePlot(myh11, features= "tdtomato", cols=c("lightgrey", "red"),
            order=T,
            reduction="umap_pca_integrated")

DimPlot(myh11,  group.by="integrated_snn_res.1",label=T) + ggsci::scale_color_simpsons()
VlnPlot(myh11, features="tdtomato", group.by="integrated_snn_res.1")

DimPlot(myh11,  group.by="predict_pca_celltype",label=T) + ggsci::scale_color_simpsons()
VlnPlot(myh11, features="tdtomato", group.by="predict_pca_celltype")

DimPlot(myh11,  group.by="predict_umap_celltype",label=T) + ggsci::scale_color_simpsons()
VlnPlot(myh11, features="tdtomato", group.by="predict_umap_celltype")


dev.off()



pdf("viz/myh11_totomato_with_imputation.pdf")
DefaultAssay(myh11) <- "MAGIC_gene.tdtomato"
FeaturePlot(myh11, features= "tdtomato", cols=c("lightgrey", "red"),
            order=T,
            reduction="umap_pca_integrated")

DimPlot(myh11,  group.by="integrated_snn_res.1",label=T) + ggsci::scale_color_simpsons()
VlnPlot(myh11, features="tdtomato", group.by="integrated_snn_res.1")

DimPlot(myh11,  group.by="predict_pca_celltype",label=T) + ggsci::scale_color_simpsons()
VlnPlot(myh11, features="tdtomato", group.by="predict_pca_celltype")

DimPlot(myh11,  group.by="predict_umap_celltype",label=T) + ggsci::scale_color_simpsons()
VlnPlot(myh11, features="tdtomato", group.by="predict_umap_celltype")
dev.off()

