library(ggplot2)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)

metadata <- read.csv(file = "../../1Round_Peakcalling/data/singlecell.csv",
						 header = TRUE, row.names = 1)
fragment.path <- "../../1Round_Peakcalling/data/fragments_ordered.tsv.bgz"

# extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

# create a gene by cell matrix
gene.activities <- FeatureMatrix(fragments = fragment.path,
                                 features = genebodyandpromoter.coords,
                                 cells = rownames(metadata))

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
gene.activities <- gene.activities[rownames(gene.activities) != "", ]

saveRDS(gene.activities, file = "GeneActivity.Rds")
