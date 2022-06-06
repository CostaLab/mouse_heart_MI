library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(optparse)
library(dplyr)
library(glue)
library(tibble)


celltp <- "Macrophages"


AllOptions <- function(){
    parser <- OptionParser()
    #parser <- add_option(parser, c("-c", "--celltype"), type="character", default="Macrophages",
    parser <- add_option(parser, c("-c", "--celltype"), type="character", default="Cardiomyocytes",
                    help="timepoint [default %default]",
                    metavar="character")
    return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)

celltp   = pa$celltype


tbl_txt <- glue("../../Diff/{celltp}/{celltp}_statistics.txt")


df <- read.csv(tbl_txt, sep="\t", header = TRUE)

all_TimePoints <- sub("TC_", "", names(df)[grepl("TC_", names(df))])


keep_vec <- c()
for (TimePoint in all_TimePoints){
    keep_vec <- c(keep_vec, sprintf("%s_%s", celltp, TimePoint))
    df[, sprintf("%s_%s", celltp, TimePoint)] <- df[, sprintf("Protection_Score_%s", TimePoint)] +
                                df[, sprintf("TC_%s", TimePoint)]
}


df <- subset(df, select = c("Motif", "Num", keep_vec))

col_order <-c("Motif", "Num", paste0(celltp, c("_Healthy", "_Day3", "_Day10")))
df <- df[, col_order]




df$Mean <- apply(df[3:ncol(df)], 1, mean)
df$Var <- apply(df[3:ncol(df)], 1, sd)
rownames(df) <- df$Motif

sel_tfs=grep("var",rownames(df),invert=TRUE)

df=df[sel_tfs,]


pdf(sprintf("%s_Num_Mean_Var.pdf", celltp))
hist(df$Num)
hist(df$Mean)
hist(df$Var)
dev.off()

write.table(df, file = sprintf("%s_TF.txt", celltp),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)


var<- subset(df$Var, df$Num > 1000)
Num=df$Num
df$Motif <- NULL
df$Mean <- NULL
df$Var <- NULL

df_plot <- subset(df, Num > 1000)
df_plot$Num <- NULL
df$Num <- NULL

df_plot_scale <- t(apply(df_plot, 1, scale))
colnames(df_plot_scale) <- colnames(df_plot)
rownames(df_plot_scale) <- rownames(df_plot)

df_filter = as.data.frame(df_plot_scale[var>0.02, ])

write.csv(df_filter, file = sprintf("%s_TF_filter_0.02.csv", celltp), quote = FALSE)


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

pdf(sprintf("%s_heatmap_with_name_0.02.pdf", celltp), width = 8, height = 40)
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

pdf(sprintf("%s_heatmap_without_name_0.02.pdf", celltp), width = 4, height = 6)
print(p)
dev.off()



TimePoints <- names(df_filter)


all_names <- c()
df_list = list()
for (TimePoint in TimePoints){

   adf <- df_filter %>% rownames_to_column("rowname") %>%
                        top_n(10, df_filter[TimePoint])
   adf <- adf  %>% filter(!(rownames(adf) %in% all_names))

   if (nrow(adf) == 0){
        next
   }
   adf <- adf %>% column_to_rownames("rowname")
   all_names <- c(all_names, rownames(adf))
   df_list[[TimePoint]] <- adf
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

pdf(glue("{celltp}_heatmap_with_name_top10.pdf"), width = 6, height = 8)
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

pdf(glue("{celltp}_heatmap_without_name_top10.pdf"), width = 6, height = 8)
print(p)
dev.off()

write.table(df, file = glue("{celltp}_TF_top_10.txt"), sep = "\t", quote = FALSE)
write.table(unique(all_names), file = glue("{celltp}_top10name.txt"), sep = "\t", quote = FALSE)





all_names <- c()
df_list = list()
for (TimePoint in TimePoints){

   adf <- df_filter %>% rownames_to_column("rowname") %>%
                        top_n(5, df_filter[TimePoint])
   adf <- adf  %>% filter(!(rownames(adf) %in% all_names))

   if (nrow(adf) == 0){
        next
   }
   adf <- adf %>% column_to_rownames("rowname")
   all_names <- c(all_names, rownames(adf))
   df_list[[TimePoint]] <- adf
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

pdf(glue("{celltp}_heatmap_with_name_top5.pdf"), width = 6, height = 8)
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

pdf(glue("{celltp}_heatmap_without_name_top5.pdf"), width = 6, height = 8)
print(p)
dev.off()

write.table(df, file = glue("{celltp}_TF_top_5.txt"), sep = "\t", quote = FALSE)
write.table(unique(all_names), file = glue("{celltp}_top10name.txt"), sep = "\t", quote = FALSE)



