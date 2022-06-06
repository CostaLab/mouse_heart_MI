library(Seurat)
library(optparse)
library(glue)
library(Signac)


integrated_name <- function(nms, time.ident="day3"){
    dict <-  c("sham_1"  = 1,
               "sham_2"  = 6,
               "day3_1"  = 2,
               "day3_2"  = 3,
               "day10_1" = 4,
               "day10_2" = 5)

    sub1 = dict[glue("{maps[TimePoint]}_{1}")]
    sub2 = dict[glue("{maps[TimePoint]}_{2}")]
    nms <- gsub("-2", glue("-{sub2}"), nms)
    nms <- gsub("-1", glue("-{sub1}"), nms)

    return(nms)
}




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

maps <- c(

    "Healthy" = "sham",
    "Day3" = "day3",
    "Day10" = "day10"
)


dir_maps <- c(

    "Healthy" = "sham",
    "Day3" = "day3_myo",
    "Day10" = "day10"
)




promo_file = glue("../../Seurat3/PromoterAccessability_timePoint/{TimePoint}/PromoterAccessability.Rds")
fmtx = glue("../../2Round_Peakcalling/data/{maps[TimePoint]}/SharedPeaks_matrix.Rds")
meta_file  <- "save/meta.csv"


meta <- read.csv(meta_file, row.names=1)

mtx <- readRDS(file=fmtx)
nnms <- integrated_name(colnames(mtx), maps[TimePoint])
colnames(mtx) <- nnms


valid_columns <- intersect(rownames(meta),  colnames(mtx) )
atac <- CreateSeuratObject(counts=mtx[, valid_columns],
                           assay = "peaks", 
                           meta.data = meta[valid_columns, ])



atac$condition = atac$time.ident
atac$condition[atac$condition == "Sham"]  = "Healthy"




##UMAP from integration timepoint
f_ref_atac <- glue("../../Seurat3/IntegrationWithscOpen_timePoint/{dir_maps[TimePoint]}/umap_integrated_annotated.Rds")
ref_atac <- readRDS(file=f_ref_atac)




pca_redu = ref_atac@reductions$pca_integrated@cell.embeddings
nnms <- integrated_name(rownames(pca_redu), maps[TimePoint])
rownames(pca_redu) <- nnms
atac[["pca_integrated"]] <- CreateDimReducObject(embeddings = pca_redu[valid_columns, ],
                         key = "PC_", assay = DefaultAssay(atac))


umap_redu <- ref_atac@reductions$umap_pca_integrated@cell.embeddings
nnms <- integrated_name(rownames(umap_redu), maps[TimePoint])
rownames(umap_redu) <- nnms
atac[["umap_pca_integrated"]] <- CreateDimReducObject(embeddings = umap_redu[valid_columns, ],
                         key = "UMAP_", assay = DefaultAssay(atac))




atac <- NormalizeData(atac)
atac <- FindVariableFeatures(atac)
atac <- ScaleData(atac)
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac, reduction.name = "Peak_LSI")
atac <- RunUMAP(atac, reduction = "Peak_LSI", dims = 2:30, reduction.name="Peak_UMAP")


promoter <- readRDS(promo_file)
nnms <- integrated_name(colnames(promoter), maps[TimePoint])
colnames(promoter) <- nnms

promoter <- promoter[, colnames(atac)]


atac[["promoter"]] <- CreateAssayObject(counts = promoter)
DefaultAssay(atac) <- "promoter"
atac <- NormalizeData(atac)
atac <- FindVariableFeatures(atac)
atac <- ScaleData(atac)
atac <- RunPCA(atac, npcs = 30, verbose = T, reduction.name="Promoter_PCA")
atac <- RunUMAP(atac, reduction = "Promoter_PCA", dims = 1:30, reduction.name="Promoter_UMAP")

saveRDS(atac, glue("save/{TimePoint}/atac_with_promoter.Rds"))




