library(ggplot2)
library(stringr)
library(magrittr)
library(WriteXLS)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(Seurat)
library(cluster)
library(clustree)
library(mclust)
library(cowplot)
library(gridExtra)
library(ggrastr)
library(viridis)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Mmusculus.UCSC.mm10)
#library(EnsDb.Hsapiens.v86)
library(data.table)
library(patchwork)
library(foreach)
library(doParallel)
library(optparse)
library(Matrix)



AllOptions <- function(){
    parser <- OptionParser()
    parser <- add_option(parser, c("-s", "--SAMPLE"), type="character", default="day3",
                    help="sample name for ATAC integration [default %default]",
                    metavar="character")
    return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)
sample <- pa$SAMPLE


dir.out <- paste0("../data/", sample)
if(!dir.exists(dir.out)){
    dir.create(dir.out)
}


## processing ATAC-seq data

fragment.path <- file.path(dir.out,"Fragments/fragments.tsv")
cell_df <- read.csv(file.path(dir.out, "singlecell.csv"), row.names=1)
cells <- rownames(cell_df)


## peak calling for each predicted label

fragment_dir <- paste0(dir.out, "/Fragments")
peaks_dir <- paste0(dir.out, "/Peaks")

if(!dir.exists(fragment_dir)){
	dir.create(fragment_dir)
}

if(!dir.exists(peaks_dir)){
    dir.create(peaks_dir)
}


fragment_file_filtered <- paste0(fragment_dir, "/", sample, ".tsv")
#FilterFragments(fragment.path = fragment.path,
#                    cells = cells,
#                    output.path = fragment_file_filtered)
    
command <- sprintf("macs2 callpeak -g mm --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all", 
	    sample, 
	    paste0(fragment_file_filtered, ".bgz"), 
	    peaks_dir)
    message("Running Macs2...")
	message(command)
	#system(command, intern = TRUE)



sessionInfo()

