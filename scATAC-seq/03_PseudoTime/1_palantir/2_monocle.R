suppressPackageStartup(library(monocle))
suppressPackageStartup(library(reshape2))
suppressPackageStartup(library(Seurat))
suppressPackageStartup(library(cicero))
suppressPackageStartup(library(glue))


orig_atac <- readRDS(file="../../Seurat3/IntegrationWithscOpen/obj.integrated_annotated_day3_myo.Rds")
#pd <- new("AnnotatedDataFrame", data = meta)
#cds <- newCellDataSet(as.matrix(mtx), phenoData = pd)

topcelltypes <- c("Fibroblasts", "Macrophages")
topcelltypes <- c("Macrophages")


for(topcelltype in topcelltypes){

  mtx <- readRDS(file=glue("../save/{topcelltype}_matrix.Rds"))
  meta <- orig_atac@meta.data[colnames(mtx), ]

  
  df <- melt(as.matrix(mtx))
  names(df) <- c("Peak", "Cell", "Count")
  

  input_cds <- make_atac_cds(df, binarize = TRUE)
  
  
  pData(input_cds) <- cbind(pData(input_cds), meta[row.names(pData(input_cds)),"time.ident"])
  pData(input_cds)$cell <- NULL
  

  agg_cds <- aggregate_nearby_peaks(input_cds, distance = 10000)

  agg_cds <- detectGenes(agg_cds)

  agg_cds <- estimateSizeFactors(agg_cds)

  agg_cds <- estimateDispersions(agg_cds)

  
  agg_cds@phenoData@data$timepoint <- agg_cds@phenoData@data$`meta[row.names(pData(input_cds)), "time.ident"]` 
  
  

  
  diff_timepoint <- differentialGeneTest(agg_cds,
                        fullModelFormulaStr="~timepoint + num_genes_expressed")
  

  # We chose a very high q-value cutoff because there are
  # so few sites in the sample dataset, in general a q-value
  # cutoff in the range of 0.01 to 0.1 would be appropriate
  ordering_sites <- row.names(subset(diff_timepoint, qval < .5))
  length(ordering_sites)
  

  agg_cds <- setOrderingFilter(agg_cds, ordering_sites)

  agg_cds <- reduceDimension(agg_cds, max_components = 2,
            residualModelFormulaStr="~num_genes_expressed",
            reduction_method = 'DDRTree')

  agg_cds <- orderCells(agg_cds)

  saveRDS(agg_cds, file=glue("{topcelltype}_agg_cds.Rds"))
}

