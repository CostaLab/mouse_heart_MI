###{r}
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
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100000 * 1024^2)


mute(library(Rcpp))



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


rna_name <- file.path("save", TimePoint,"aggred_a_rna.Rds")
rna_mtx <- readRDS(file=rna_name)


message("rna ", date())
rna_seurat <- CreateSeuratObject(counts = rna_mtx)
rna_seurat <- NormalizeData(rna_seurat)
rna_seurat <- FindVariableFeatures(rna_seurat)
rna_seurat <- ScaleData(rna_seurat, features=rownames(rna_seurat))
scaled_rna_mtx <- as.matrix(rna_seurat@assays$RNA@scale.data)


message("Scaling atac ", date())
atac_name <- file.path("save", TimePoint,"aggred_a_atac.Rds")
atac_mtx <- readRDS(file=atac_name)
atac_mtx <- atac_mtx[, colnames(rna_mtx)]




atac_seurat <- CreateSeuratObject(counts = atac_mtx)
atac_seurat <- NormalizeData(atac_seurat)
atac_seurat <- FindVariableFeatures(atac_seurat)
atac_seurat <- ScaleData(atac_seurat, features=rownames(atac_seurat))
scaled_atac_mtx <- as.matrix(atac_seurat@assays$RNA@scale.data)


saveRDS(rna_seurat, file=glue("save/{TimePoint}/seurat_scale_rna_mtx.Rds"))
saveRDS(atac_seurat, file=glue("save/{TimePoint}/seurat_scale_atac_mtx.Rds"))
