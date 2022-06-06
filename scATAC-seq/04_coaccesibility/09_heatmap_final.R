suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr)) ## for str_remove
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(fclust))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(tiff))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(magick))

plan("multiprocess", workers = 20)
options(future.globals.maxSize = 80000 * 1024^2)
#

source("config_heatmap_orders.R")


#----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

parser <- OptionParser()


parser <- add_option(parser, c("-t", "--TimePoint"), type="character", default="Day3",
                    help="sample name for ATAC integration [default %default]",
                    metavar="character")

parser <- add_option(parser, c("-c", "--color"), type="character", default="A",
                     help="colors option [default %default]",
                     metavar="CHAR")
parser <- add_option(parser, c("-p", "--plots"), type="character", default="0.6_0.05",
                    help="plots for different pval_option [default %default]",
                    metavar="character")

cols = c('Cardiomyocytes' = '#800000',
         'Fibroblasts' = '#911eb4',
         'Myofibroblasts' = '#f032e6',
         'Endothelial cells' = '#000075',
         'Macrophages' = '#e6194B',
         'Anti-inflammatory macrophages' = '#568198',
         'Anti−inflammatory macrophages' = '#568198',
         'inflammatory macrophages' = '#469990',
         'Pericytes/vSMC' = '#aaffc3')



pa = parse_args(parser)

color.option <- pa$color 
TimePoint = pa$TimePoint

rna_name = TimePoint

plots_opt <-  ifelse(pa$plots == "",  "",  paste0("_", pa$plots))


print(plots_opt)
dir.create(sprintf("plots%s", plots_opt))
print(sprintf("plots%s", plots_opt))



atac_mtx <- readRDS(file=glue("save/{TimePoint}/scale_atac_mtx.Rds"))

atac_rownames <- rownames(atac_mtx)
atac_rownames <- gsub("-", "_", atac_rownames)
rownames(atac_mtx) <- atac_rownames

rna_mtx <- readRDS(file=glue("save/{TimePoint}/scale_rna_mtx.Rds"))
fdf <- read.csv(file=glue("save/{TimePoint}/corr_pval_final{plots_opt}.tsv"), sep="\t", stringsAsFactors = F)
order_atac_mtx <- as.matrix(atac_mtx[fdf$peak, ])#[1:1000, 1:1000]


cluster_anno <- get_cluster_anno(plots_opt, TimePoint) 

cluster_orders <- as.integer(names(cluster_anno)) 


number.clusters = length(cluster_orders)

row_order_list <- readRDS(glue("plots{plots_opt}/row_order_celltype_{number.clusters}_{TimePoint}.Rds"))



rorder <- as.character(cluster_orders)
nlist <- row_order_list[rorder]
row_order_ <- unlist(nlist)
row_order_anno <-lapply(seq_along(nlist),
                function(y, n, i) { rep(n[[i]] ,length(y[[i]])) },
                y=nlist, n=names(nlist))

number.clusters <- length(row_order_list)



col_order_list <- readRDS(glue("plots{plots_opt}/col_order_celltype_{number.clusters}_{TimePoint}.Rds"))

col_order_ <- unlist(col_order_list)


#colors <- c("#f7fcf0",  "#e0f3db",  "#ccebc5",  "#a8ddb5",  "#7bccc4",  "#4eb3d3",  "#2b8cbe",  "#0868ac",  "#084081")
#col_fun = colorRamp2(c(-4:4), colors)


cn <-unlist(col_order_list)

atac_dir <- glue("save/{TimePoint}/atac_with_promoter.Rds")

atac <- readRDS(atac_dir)



meta_data <- atac@meta.data
meta_data <- meta_data[which(rownames(meta_data) %in% colnames(order_atac_mtx)), ]
meta_data$celltype <- factor(meta_data$celltype, levels=names(cols)) 
meta_data$celltype <- droplevels(meta_data$celltype)
meta_data <- meta_data[order(meta_data$celltype), ]
meta_data$celltype <- as.character(meta_data$celltype) 

celltypes <- as.character(meta_data[cn,]$celltype)

table(celltypes)


#        colors <- c("#9A6324","#800000")
#        names(colors) <- c(cells_1_name, cells_2_name)
#        #colors <- factor(colors, levels=c("#800000", "#9A6324"))
#        aux <- ColorMapping(name="celltype", colors=colors, levels=names(colors))
#        ta = HeatmapAnnotation(celltype=col_category, 
#              annotation_name_gp = gpar(fontsize=0), 
#              col=list(celltype=colors),
#              annotation_legend_param = list(object = aux, plot=F))
#

aux <- ColorMapping(name = "celltype", colors=cols[unique(celltypes)], levels = names(cols[unique(celltypes)]))

ha = HeatmapAnnotation(celltype = celltypes,
                         annotation_name_gp = gpar(fontsize = 0),
                         col = list(celltype = cols),
                         annotation_legend_param = list(object = aux, plot=F))

lens <- sapply(row_order_anno, length)
x = 1
vec = c(x)
for (i in 1:(length(lens)-1)){
    x = x + lens[i]
    vec = c(vec, x)
}


annotated_cluster <- cluster_anno[unlist(row_order_anno)]

ca = rowAnnotation(clusters = annotated_cluster,
                   annotation_name_gp = gpar(fontsize = 0),
                   col = list(clusters = cols[unique(annotated_cluster)]),
                   show_legend = FALSE
                   )



all_peaks <- rownames(order_atac_mtx)[row_order_]


col_fun <-  circlize::colorRamp2(c(-2, 0, +2), c("purple", "black", "yellow"))

pdf(glue("plots{plots_opt}/reorder_{color.option}_{TimePoint}_heatmap_{number.clusters}.pdf"), width = 8, height = 8)

h1 <- Heatmap(order_atac_mtx[all_peaks, col_order_],
          cluster_rows = F, show_row_dend = F,
          #split = cluster_split,
          #row_order=r_o,
          #row_title = str_c(rev(rorder), collapse=" <- "),
          show_column_dend = F, cluster_columns = F,
          show_row_names = F,
          show_column_names = F,
          column_title = "Enhancer",
          name = sprintf("Enhancer %s", TimePoint),
          top_annotation = ha,
          col = col_fun,
          left_annotation = ca,
          use_raster=T,
          raster_device="CairoTIFF"
        )



#draw(h1)

#all_rna <- atac
#DefaultAssay(all_rna) <- "promoter"
#Idents(all_rna) <- "condition"
#rna <- subset(all_rna, ident=TimePoint)
#???? rna



col_order_i <- match(col_order_, colnames(order_atac_mtx))
rna_metadata <- atac@meta.data

cn_rna <- colnames(rna_mtx)
cn_rna <- str_remove(cn_rna, "\\.\\d+$")
rna_celltypes  <- as.character(rna_metadata[cn_rna,]$celltype)
col_fun <-  circlize::colorRamp2(c(-2, 0, +2), c("purple", "black", "yellow"))
ha = HeatmapAnnotation(rna_celltype = rna_celltypes[col_order_i],
                         show_legend = F,
                         annotation_name_gp = gpar(fontsize = 0),
                         col = list(rna_celltype = cols))

h2 <- Heatmap(rna_mtx[fdf$gene[row_order_], col_order_i],
          cluster_rows = F, show_row_dend = F,
          cluster_columns = F, show_column_dend = F,
          show_row_names = F,
          column_title = "Promoter",
          show_column_names = F,
          name = sprintf("Promoter %s", TimePoint),
          top_annotation=ha,
          col = col_fun,
          use_raster=T,
          raster_device="CairoTIFF"
        )

draw(h1 + h2)
dev.off()


peaks_clusters <- lapply(row_order_list, function(x)rownames(order_atac_mtx)[x] )
genes_clusters <- lapply(row_order_list, function(x) rownames(rna_mtx[fdf$gene, ])[x]) 

merge_info <- Map(list, peaks_clusters, genes_clusters)
fname_merge_info <- glue("plots{plots_opt}/merge_info_{TimePoint}_{number.clusters}.Rds")
saveRDS(merge_info, file=fname_merge_info)



fname_atac <- glue("plots{plots_opt}/final_mtx_atac_{TimePoint}_{number.clusters}.Rds")
fname_rna <- glue("plots{plots_opt}/final_mtx_rna_{rna_name}_{number.clusters}.Rds")
saveRDS(order_atac_mtx[all_peaks, col_order_], file=fname_atac)
saveRDS(rna_mtx[fdf$gene[row_order_], col_order_i], file=fname_rna)

