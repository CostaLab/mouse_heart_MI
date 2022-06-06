library(optparse)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(reshape2)
library(rlist)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dplyr)
library(glue)
library(tibble)



AllOptions <- function(){
    parser <- OptionParser()
    parser <- add_option(parser, c("-t", "--timepoint"), type="character", default="Day10",
                    help="timepoint [default %default]",
                    metavar="character")
    return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)

tm   = pa$timepoint


celltypes = c('Cardiomyocytes',
              'Fibroblasts',
              'Myofibroblasts',
              'Endothelial cells',
              'Macrophages',
              'Anti-inflammatory macrophages',
              'inflammatory macrophages',
              'Pericytes',
              'Pericytes.vSMC') 

all_tfs <- c()
for (celltype in celltypes){
    Ivs_dir <- glue("../../../Diff1VS_myo/{tm}/{celltype}/{tm}_statistics.txt")
    if(!file.exists(Ivs_dir)){
        next 
    }

    df <- read.csv(Ivs_dir,sep="\t", header=T)
    rownames(df) <- df$Motif
    a_sel_tf <- rownames(subset(df, P_values < 0.05 & Num > 100))
    all_tfs <- c(all_tfs, a_sel_tf)
}


sel_tf <- unique(all_tfs)
sel_tf <- sel_tf[grep("var", sel_tf, invert=TRUE)]


tbl_txt <- sprintf("../Diff/%s_statistics.txt", tm) 

df <- read.csv(tbl_txt, sep="\t", header = TRUE)

all_celltypes <- sub("TC_", "", names(df)[grepl("TC_", names(df))])


keep_vec <- c() 
for (celltype in all_celltypes){
    keep_vec <- c(keep_vec, sprintf("%s_%s", tm, celltype))
    df[, sprintf("%s_%s", tm, celltype)] <- df[, sprintf("Protection_Score_%s", celltype)] +
                                df[, sprintf("TC_%s", celltype)]
}


df <- subset(df, select = c("Motif", "Num", keep_vec))



df_filter <- df %>% filter(Motif %in% sel_tf) 

rownames(df_filter) <- df_filter$Motif
df_filter$Motif <- NULL
df_filter$Num <- NULL


write.csv(df_filter, file = sprintf("pval_%s_TF_filter.csv", tm), quote = FALSE)

df_filter_scale <- t(apply(df_filter, 1, scale))
colnames(df_filter_scale) <- colnames(df_filter)
rownames(df_filter_scale) <- rownames(df_filter)

df_filter <- as.data.frame(df_filter_scale)


#mtx <- as.matrix(df_filter)
p <- Heatmap(as.matrix(df_filter),
             name = "Z-score",
             cluster_columns = TRUE,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 15),
             clustering_method_rows = "centroid",
             clustering_method_columns = "centroid",
             clustering_distance_rows = "pearson",
             clustering_distance_columns = "pearson",
             col = rev(brewer.pal(n = 11, name = "RdBu")))

pdf(sprintf("pval_%s_heatmap_with_name.pdf", tm), width = 8, height = 40)
print(p)
dev.off()

p <- Heatmap(as.matrix(df_filter),
             name = "Z-score",
             cluster_columns = TRUE,
             show_row_names = FALSE,
             show_column_dend = TRUE,
             show_row_dend = FALSE,
             row_names_gp = gpar(fontsize = 4),
             column_names_gp = gpar(fontsize = 10),
             clustering_method_rows = "centroid",
             clustering_method_columns = "centroid",
             clustering_distance_rows = "pearson",
             clustering_distance_columns = "pearson",
             col = rev(brewer.pal(n = 11, name = "RdBu")))

pdf(sprintf("pval_%s_heatmap_without_name.pdf", tm), width = 4, height = 6)
print(p)
dev.off()



celltypes <- names(df_filter)


all_names <- c()
df_list = list()
for (celltype in colnames(df_filter)){

   adf <- df_filter %>% rownames_to_column("rowname") %>% 
                        top_n(10, df_filter[celltype])
   adf <- adf  %>% filter(!(rownames(adf) %in% all_names))

   if (nrow(adf) == 0){
        next
   }
   adf <- adf %>% column_to_rownames("rowname")
   all_names <- c(all_names, rownames(adf))
   df_list[[celltype]] <- adf
}

names(df_list) <- NULL
df <- do.call(rbind, df_list)




#f1 <-  circlize::colorRamp2(c(-4, 0, +4), c("purple", "black", "yellow"))
f1 <-  circlize::colorRamp2(c(-2, 0, +2), c("purple", "black", "yellow"))

p <- Heatmap(as.matrix(df),
             name = "TF Activity",
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             row_names_gp = gpar(fontsize = 8), 
             column_names_gp = gpar(fontsize = 15),
             clustering_method_rows = "single",
             clustering_distance_rows = "pearson",
             col = f1) 

pdf(glue("pval_{tm}_heatmap_with_name_top10.pdf"), width = 6, height = 8)
print(p)
dev.off()


p <- Heatmap(as.matrix(df),
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             row_names_gp = gpar(fontsize = 8), 
             column_names_gp = gpar(fontsize = 15),
             clustering_method_rows = "single",
             clustering_distance_rows = "pearson",
             col = f1, 
             show_column_names = FALSE,
             show_row_names = FALSE,
             show_heatmap_legend = FALSE)

pdf(glue("pval_{tm}_heatmap_without_name_top10.pdf"), width = 6, height = 8)
print(p)
dev.off()

write.table(df, file = glue("pval_{tm}_TF_top_10.txt"), sep = "\t", quote = FALSE)
write.table(unique(all_names), file = glue("pval_{tm}_top10name.txt"), sep = "\t", quote = FALSE)





all_names <- c()
df_list = list()
for (celltype in colnames(df_filter)){

   adf <- df_filter %>% rownames_to_column("rowname") %>% 
                        top_n(5, df_filter[celltype])
   adf <- adf  %>% filter(!(rownames(adf) %in% all_names))

   if (nrow(adf) == 0){
        next
   }
   adf <- adf %>% column_to_rownames("rowname")
   all_names <- c(all_names, rownames(adf))
   df_list[[celltype]] <- adf
}

names(df_list) <- NULL
df <- do.call(rbind, df_list)




#f1 <-  circlize::colorRamp2(c(-4, 0, +4), c("purple", "black", "yellow"))
f1 <-  circlize::colorRamp2(c(-2, 0, +2), c("purple", "black", "yellow"))

p <- Heatmap(as.matrix(df),
             name = "TF Activity",
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             row_names_gp = gpar(fontsize = 8), 
             column_names_gp = gpar(fontsize = 15),
             clustering_method_rows = "single",
             clustering_distance_rows = "pearson",
             col = f1) 

pdf(glue("pval_{tm}_heatmap_with_name_top5.pdf"), width = 6, height = 8)
print(p)
dev.off()


p <- Heatmap(as.matrix(df),
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             row_names_gp = gpar(fontsize = 8), 
             column_names_gp = gpar(fontsize = 15),
             clustering_method_rows = "single",
             clustering_distance_rows = "pearson",
             col = f1, 
             show_column_names = FALSE,
             show_row_names = FALSE,
             show_heatmap_legend = FALSE)

pdf(glue("pval_{tm}_heatmap_without_name_top5.pdf"), width = 6, height = 8)
print(p)
dev.off()

write.table(df, file = glue("pval_{tm}_TF_top_5.txt"), sep = "\t", quote = FALSE)
write.table(unique(all_names), file = glue("pval_{tm}_top5name.txt"), sep = "\t", quote = FALSE)


