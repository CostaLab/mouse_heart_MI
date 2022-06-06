mute <- suppressPackageStartupMessages
mute(library(Matrix))
mute(library(Seurat))
mute(library(optparse))
mute(library(RANN))
mute(library(data.table))
mute(library(glue))
mute(library(foreach))
mute(library(doParallel))
registerDoParallel(cores=3)
set.seed(1)

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


atac_mtx_undf <- readRDS(file=glue("save/{TimePoint}/scale_atac_mtx.Rds"))
rna_mtx_undf <- readRDS(file=glue("save/{TimePoint}/scale_rna_mtx.Rds"))

atac_rownames <- rownames(atac_mtx_undf)
atac_rownames <- gsub("-", "_", atac_rownames)
rownames(atac_mtx_undf) <- atac_rownames



corr_df <- read.csv(glue("save/{TimePoint}/corr.tsv"), sep="\t", stringsAsFactors=F)
chromosome <- read.csv(glue("save/{TimePoint}/chromosome_range.txt"),
                            sep = " ", stringsAsFactors = F)
nr <- nrow(chromosome)
maxidx <- chromosome$end[nr]
whole_seq = seq(1, maxidx)


set.seed(1)
chr_all_corrs <- foreach(i = 1:nrow(chromosome)) %dopar% {
  a_chr <- chromosome[i, ]
  tseq <- setdiff(whole_seq, a_chr$start:a_chr$end)
  peaks_idx <- sample(tseq)[1:1000]
  peaks <- corr_df[peaks_idx, ]$peak
  
  genes <- strsplit(a_chr$genes, ",")[[1]]
  genes <- intersect(genes, rownames(rna_mtx_undf))
    

  peak_mtx_idx = match(peaks, rownames(atac_mtx_undf)) 
  gene_mtx_idx = match(genes, rownames(rna_mtx_undf))

  grid <- expand.grid(peak_mtx_idx, gene_mtx_idx)
  all_corrs <- rowCorCpp(as.integer(grid[,1]), as.integer(grid[,2]), atac_mtx_undf, rna_mtx_undf)

  all_corrs = as.vector(all_corrs)
  all_corrs = na.omit(all_corrs)
}


names(chr_all_corrs) <-  chromosome$chromosome

message(paste0("  finished calculate: ", date()))
fname = glue("save/{TimePoint}/null_hypothesis.Rds")
saveRDS(chr_all_corrs, file=fname)
message(paste0("  written: ", date()))
#write.table(new_ann, file = fname, sep=",",  quote=FALSE)
