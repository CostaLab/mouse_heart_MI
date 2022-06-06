library(tidyverse)
library(Seurat)
library(glue)
library(Rsamtools)
library(data.table)

celltype_name = "macro"

output_dir  = file.path("../save", celltype_name )
dir.create(output_dir)


celltype_group = list(
                    "fib" =  c("Fibroblasts", "Myofibroblasts", "Pericytes/vSMC"),
                    "macro" = c("Macrophages", "Anti-inflammatory macrophages", "inflammatory macrophages")
) 

###########################################################

compress_n_idx <- function(fname){
  outf <- bgzip(file = fname)
  message("indexing", date())
  index.file <- indexTabix(file = paste0(outf),
                   format = 'bed',
                   zeroBased = TRUE)

}


ref_atac_file <- "../../../Seurat3/IntegrationWithscOpen/obj.integrated_annotated_day3_myo.Rds" 
ref_atac <- readRDS(file=ref_atac_file)
meta <- ref_atac@meta.data

barcodes <- rownames(meta[meta$celltype %in% celltype_group[[celltype_name]], ])

fragment.path = "../../../1Round_Peakcalling/data/fragments_ordered.tsv"
tib <- read_tsv(fragment.path, col_names=F)


sub_tib <- tib[tib$X4 %in% barcodes,  ]
write_tsv(sub_tib, file=file.path(output_dir, "fragments.tsv"), col_names=FALSE)
compress_n_idx(file.path(output_dir, "fragments.tsv"))



