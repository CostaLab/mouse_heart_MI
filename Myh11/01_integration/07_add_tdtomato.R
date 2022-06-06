library(Seurat)
library(Rmagic)



myh11 <- readRDS(file="save_res/trans.Rds")
cnms <- colnames(myh11)
tdtomato <- rep(0, length(cnms))
names(tdtomato) <- cnms


td_df <- readRDS("../../data/cellranger_tdTomato_combined/src/tdtomato_counts.Rds")
dict_v <- td_df$Freq
names(dict_v) <- td_df$Var1
notin_idx <- which(!(td_df$Var1 %in% colnames(myh11)))
dict_v <- dict_v[-notin_idx]


tdtomato[names(dict_v)] <- dict_v
gene_access <- myh11@assays$gene.activity@counts
gene_access <- rbind(gene_access, tdtomato)
myh11[["gene.tdtomato"]] <- CreateAssayObject(counts=gene_access)

DefaultAssay(myh11) <- "gene.tdtomato"
all_genes <- rownames(myh11)
myh11 <- magic(myh11, genes=all_genes)


saveRDS(myh11, file="save_res/trans_add_tdtomato.Rds")


