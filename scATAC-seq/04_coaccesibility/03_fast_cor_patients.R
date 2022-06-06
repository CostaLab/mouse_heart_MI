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


Rcpp::sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX)-1 > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY)-1 > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }

    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);

    const int nx = X.ncol();
    const int ny = Y.ncol();

    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }

    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }

    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }

    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }

    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);    }

    return(cor);

  }'
)


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


saveRDS(scaled_rna_mtx, file=glue("save/{TimePoint}/scale_rna_mtx.Rds"))
saveRDS(scaled_atac_mtx, file=glue("save/{TimePoint}/scale_atac_mtx.Rds"))


#scaled_rna_mtx <- readRDS(file=glue("save/{TimePoint}/scale_rna_mtx.Rds"))
#scaled_atac_mtx <- readRDS(file=glue("save/{TimePoint}/scale_atac_mtx.Rds"))


new_ann <- data.frame(peak=character(), gene=character(), 
                      distance=numeric(), #peak_type=character(), 
                      corr=numeric(), stringsAsFactors = F)

ann <- read.csv(file=glue("save/{TimePoint}/peak2gene_putative-simple.tsv"), sep="\t", stringsAsFactors = F)
###

###{r}
#future_sapply()
#scaled_atac_mtx <- as.matrix(atac_mtx)
#scaled_rna_mtx <- as.matrix(rna_mtx)

#rna_genes <-rownames(scaled_rna_mtx)
#names(rna_genes) <- rna_genes


`%notin%` <- Negate(`%in%`) 

atac_rownames <- rownames(scaled_atac_mtx)
atac_rownames <- gsub("-", "_", atac_rownames)
rownames(scaled_atac_mtx) <- atac_rownames

ann <- ann %>% filter(gene %in% rownames(scaled_rna_mtx))
ann <- ann %>% filter(peak %in% rownames(scaled_atac_mtx))


peaks <- ann$peak
genes <- ann$gene

peak_idx <- match(peaks, rownames(scaled_atac_mtx))
gene_idx <- match(genes, rownames(scaled_rna_mtx))

message(date(), "start calculate correlation")
corrs <- rowCorCpp(peak_idx, gene_idx, scaled_atac_mtx, scaled_rna_mtx)
new_ann <- data.frame(peaks=ann$peak, gene=ann$gene, distance=ann$distance, corr=corrs, stringsAsFactors=F)


order_v <- factor(sapply(new_ann[["peaks"]], function(x) str_split(x, ":")[[1]][1] ),
                    levels= paste0("chr", c(1:19, "X", "Y")) )

new_ann <- new_ann[order(order_v) ,]


message( paste0(TimePoint,  "  finished calculate: " , dim(ann)[1], "  ", date()))
fname = glue("save/{TimePoint}/corr.tsv")
write.table(new_ann, file = fname, sep="\t",  quote=FALSE)
message( paste0(TimePoint,  "  saved: " , dim(ann)[1], "  ", date()))


###
