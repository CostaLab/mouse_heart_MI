mute <- suppressPackageStartupMessages
mute(library(Matrix))
mute(library(Seurat))
mute(library(optparse))
mute(library(stringr))
mute(library(dplyr))
mute(library(RANN))
mute(library(data.table))
mute(library(foreach))
mute(library(glue))
mute(library(doParallel))
registerDoParallel(cores=1)

suppressPackageStartupMessages(library(future))
plan("multiprocess", workers = 100)
options(future.globals.maxSize = 100000 * 1024^2)

AllOptions <- function(){
    parser <- OptionParser()
    parser <- add_option(parser, c("-t", "--TimePoint"), type="character", default="Healthy",
                    help="sample name for ATAC integration [default %default]",
                    metavar="character")
    return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)
TimePoint  = pa$TimePoint
dir.create(glue("save/{TimePoint}"))


ref_atac_file <- file.path("save", TimePoint, "atac_with_promoter.Rds")
ref_atac <- readRDS(file=ref_atac_file)


atac_seurat <- readRDS(file=glue("save/{TimePoint}/seurat_scale_atac_mtx.Rds"))
atac_seurat@meta.data <- ref_atac@meta.data[colnames(atac_seurat), ]

atac_seurat$celltype <- as.character(atac_seurat$celltype)
Idents(atac_seurat) <- "celltype"

de.df <- FindAllMarkers(atac_seurat)
cluster.de <- split(de.df, de.df$cluster)
saveRDS(cluster.de, file=glue("save/{TimePoint}/ScaleMtx_Peaks_differential.Rds"))

