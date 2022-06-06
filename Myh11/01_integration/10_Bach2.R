library(Seurat)
library(ggplot2)
library(Rmagic)

obj <- readRDS(file="save_res/trans_add_tdtomato_annotation_res1.0_remove_res0.1_456.Rds")
DefaultAssay(obj) <- "gene.activity"

obj$seurat_clusters <- obj$integrated_snn_res.0.1

umap_colors <- c("#fed439", "#709ae1", "#709ae1", "#d2af81")


pdf("viz/UMAPs.pdf", height=4, width=4)
DimPlot(obj, group.by="integrated_snn_res.0.1", label=F, cols=ggsci::pal_simpsons()(18))
DimPlot(obj, group.by="integrated_snn_res.0.1", label=F, cols=ggsci::pal_simpsons()(18)) + theme(legend.position='none')
DimPlot(obj, group.by="integrated_snn_res.0.1", label=F, cols=ggsci::pal_simpsons()(18)) + theme_void() + theme(legend.position='none')
DimPlot(obj, group.by="integrated_snn_res.0.1", label=F, cols=umap_colors) + theme_void() + theme(legend.position='none')
DimPlot(obj, group.by="time.ident", label=F, cols=ggsci::pal_simpsons()(18))

obj$seurat_clusters <- obj$integrated_snn_res.1
DimPlot(obj, group.by="integrated_snn_res.1", label=F, cols=ggsci::pal_simpsons()(18))

dev.off()



library(reticulate)
use_python("/home/sz753404/miniconda3/envs/schema/bin/python")

obj <- Rmagic::magic(obj, genes=rownames(obj))
DefaultAssay(obj) <- "MAGIC_gene.activity"



pdf("viz/Bach2.pdf")
Idents(obj) <- "annotation"
obj$time.ident <- factor(obj$time.ident, levels=c("Sham", "Day3", "Day10"))

##
a_sub <- subset(obj, idents=c("Fibroblasts", "Myofibroblasts"))
a_sub$annotation <- as.character(a_sub$annotation)
VlnPlot(a_sub, features="Bach2", group.by="time.ident", split.by="annotation", split.plot=T) + ggsci::scale_color_simpsons()

##
a_sub <- subset(obj, idents=c("Pericytes", "Myofibroblasts"))
a_sub$annotation <- as.character(a_sub$annotation)
VlnPlot(a_sub, features="Bach2", group.by="time.ident", split.by="annotation", split.plot=T) + ggsci::scale_color_simpsons()


##
a_sub <- subset(obj, idents=c("Pericytes", "Fibroblasts"))
a_sub$annotation <- as.character(a_sub$annotation)
VlnPlot(a_sub, features="Bach2", group.by="time.ident", split.by="annotation", split.plot=T) + ggsci::scale_color_simpsons()


dev.off()

