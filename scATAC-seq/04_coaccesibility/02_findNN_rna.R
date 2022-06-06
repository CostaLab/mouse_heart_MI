mute <- suppressPackageStartupMessages
mute(library(Matrix))
mute(library(Seurat))
mute(library(optparse))
mute(library(RANN))
mute(library(glue))
mute(library(data.table))
mute(library(doParallel))
registerDoParallel(cores=10)

FindNN <- function(
  object,
  cells1 = NULL,
  cells2 = NULL,
  grouping.var = NULL,
  dims = 1:10,
  reduction = "cca.l2",
  nn.dims = dims,
  nn.reduction = reduction,
  k = 300, 
  eps = 0, 
  integration.name = 'integrated',
  verbose = TRUE 
) {
  if (xor(x = is.null(x = cells1), y = is.null(x = cells2))) {
    stop("cells1 and cells2 must both be specified")
  }
  if (!is.null(x = cells1) && !is.null(x = cells2) && !is.null(x = grouping.var)) {
    stop("Specify EITHER grouping.var or cells1/2.")
  }
  if (is.null(x = cells1) && is.null(x = cells2) && is.null(x = grouping.var)) {
    stop("Please set either cells1/2 or grouping.var")
  }
  if (!is.null(x = grouping.var)) {
    if (nrow(x = unique(x = object[[grouping.var]])) != 2) { 
      stop("Number of groups in grouping.var not equal to 2.")
    }    
    groups <- names(x = sort(x = table(... = object[[grouping.var]]), decreasing = TRUE))
    cells1 <- colnames(x = object)[object[[grouping.var]] == groups[[1]]]
    cells2 <- colnames(x = object)[object[[grouping.var]] == groups[[2]]]
  }
  if (verbose) {
    message("Finding neighborhoods")
  }
  dim.data.self <- Embeddings(object = object[[nn.reduction]])[ ,nn.dims]
  dim.data.opposite <- Embeddings(object = object[[reduction]])[ ,dims]
  dims.cells1.self <- dim.data.self[cells1, ]
  dims.cells1.opposite <- dim.data.opposite[cells1, ]
  dims.cells2.self <- dim.data.self[cells2, ]
  dims.cells2.opposite <- dim.data.opposite[cells2, ]

  nnaa <- nn2( 
    data = dims.cells1.self,
    k = k + 1, 
    eps = eps
  )
  nnab <- nn2( 
    data = dims.cells2.opposite,
    query = dims.cells1.opposite,
    k = k, 
    eps = eps
  )
  nnbb <- nn2( 
    data = dims.cells2.self,
    k = k + 1, 
    eps = eps
  )
  nnba <- nn2( 
    data = dims.cells1.opposite,
    query = dims.cells2.opposite,
    k = k, 
    eps = eps
  )
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'neighbors',
    new.data = list('nnaa' = nnaa, 'nnab' = nnab, 'nnba' = nnba, 'nnbb' = nnbb, 'cells1' = cells1, 'cells2' = cells2)
  )
  return(object)
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


atac_file <- file.path("save", TimePoint,"atac_with_promoter.Rds")
atac <- readRDS(file=atac_file)
atac$condition <- atac$time.ident
atac$condition[atac$condition == "Sham"] = "Healthy"

Idents(atac) <- "condition"

a_atac <- subset(atac, ident = TimePoint)


DefaultAssay(a_atac) <- "peaks"



###TODO:
### 1. Need compute same celltype 
### 2. Maybe can be changed to umap_pac_integreted


## ------ neighbours without caring about celltypes
a <- colnames(a_atac)
b <- colnames(a_atac)
#rst <- FindNN(a_atac, a, b, reduction="umap_pca_integrated", k=60, dims=1:30)
rst <- FindNN(a_atac, a, b, reduction="umap_pca_integrated", k=60, dims=1:2)
nns <- rst@tools$integrated@neighbors


peak_nearest_cell_list = list()
for(i in seq_along(colnames(a_atac))){
  vec <- b[nns$nnbb$nn.idx[i, 1:51]]
  peak_nearest_cell_list[[i]] <- vec 
}


fname = file.path("save/", TimePoint ,"umap_a_atac_aggre_all.csv")
write.table(transpose(peak_nearest_cell_list), sep = "\t", 
            file = fname, row.names = T,col.names = F,
            quote=FALSE)

## ------ neighbours without caring about celltypes




### aggregate cells of same celltype
all_celltypes <- unique(a_atac$celltype)
all_list <- list()
for(ct in all_celltypes){
  
  a <- colnames(a_atac[, a_atac$celltype == ct])
  b <- colnames(a_atac[, a_atac$celltype == ct])
  rst <- FindNN(a_atac, a, b, reduction="umap_pca_integrated", k=60, dims=1:2)
  nns <- rst@tools$integrated@neighbors
  a_peak_nearest_list <- list()
  for(i in seq_along(a)){
    vec <- b[nns$nnbb$nn.idx[i, 1:51]]
    a_peak_nearest_list[[a[i]]] <- vec
  }
  all_list[[ct]] <- a_peak_nearest_list
}

combine_list <- do.call(c, all_list)

fname = file.path("save/", TimePoint,"a_rna_aggre_celltype.csv")
write.table(transpose(combine_list), sep = "\t", 
            file = fname, row.names = T,col.names = F,
            quote=FALSE)



mtx_a_rna <- GetAssayData(a_atac, slot="counts", assay="promoter")
all_mtx_list <- lapply(all_celltypes, function(ct) mtx_a_rna[, names(all_list[[ct]])])
names(all_mtx_list) <- all_celltypes


aggr_mtx_list <- foreach(i=1:length(all_celltypes)) %dopar%{
  ct <- all_celltypes[i]
  cell_num <- length(all_list[[ct]])
  mtx <- all_mtx_list[[ct]]
  lst_peak_nearest_cell  <- all_list[[ct]]
  message(date(), " ", ct)
  aggre_a_atac <- foreach(j =1:cell_num, .combine = cbind)%dopar%{
                if( j%%100 == 99){
                  message(date(), "  ", j)
                }
                one <- lst_peak_nearest_cell[[j]]
                cell_name <- one[1]
                rowSums(mtx_a_rna[, one[1:51]])
  }
  message(date(), " ok. ", ct)
  aggre_a_atac

}

names(aggr_mtx_list) <- all_celltypes
for(ct in all_celltypes){
   col_names <- colnames(all_mtx_list[[ct]])
   colnames(aggr_mtx_list[[ct]]) <- col_names
}
aggre_a_atac_mtx <- do.call(cbind, aggr_mtx_list)


print(paste("end ", date()))
saveRDS(aggre_a_atac_mtx, file=glue("save/{TimePoint}/aggred_a_rna.Rds"))

