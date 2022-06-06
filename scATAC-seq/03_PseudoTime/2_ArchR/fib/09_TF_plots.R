library(ggplot2)
library(stringr)
library(magrittr)
library(WriteXLS)
library(tidyr)
library(dplyr)
library(plotly)
library(cluster)
library(cowplot)
library(gridExtra)
library(viridis)
library(GenomicRanges)
library(GenomeInfoDb)
library(data.table)
library(ArchR)
library(glue)

## set parameters
set.seed(42)
addArchRThreads(threads =100)
addArchRGenome("mm10")


cols.celltypes <- c(
    'Fibroblasts' = '#f032e6',
    'Myofibroblasts' = '#911eb4',
    'Pericytes/vSMC' = '#e6beff')



Add_A_Matrix <- function(
    project = proj,
    mtx = mtx,
    vname = "z",
    name = "AMatrix"
){
  ArrowFiles <- getArrowFiles(proj)
  nms <- names(ArrowFiles)
  for(nm in nms){
    a_cellnm <- colnames(mtx)[grepl(nm, colnames(mtx))]
    message(nm, " ", length(a_cellnm))
    a_mtx <- mtx[, a_cellnm]
    ArrowFile <- ArrowFiles[nm]
    featureDF <-data.frame(seqnames = vname, 
                           idx = seq_len(nrow(a_mtx)), 
                           name = rownames(a_mtx),
                           stringsAsFactors = FALSE)
    
    #Initialize
   o <- ArchR:::.initializeMat(
    ArrowFile = ArrowFile,
    Group = name,
    Class = "double",
    Units =  vname,
    cellNames = a_cellnm,
    params = vname,
    featureDF = featureDF,
    force = TRUE)
    
   o <- ArchR:::.addMatToArrow(
        mat = as(a_mtx,"dgCMatrix"), 
        ArrowFile = ArrowFile, 
        Group = paste0(name, "/", vname), 
        binarize = FALSE,
        addColSums = TRUE,
        addRowSums = TRUE,
        addRowVarsLog2 = TRUE)   
    o <- h5write(obj = "Finished", file = ArrowFile, name = paste0(name,"/Info/Completed"))
  }  
  return(proj)
}


proj <- loadArchRProject(path = "./Fib", showLogo = FALSE)
proj <- addImputeWeights(proj)

### Nfe2l2 == Nfr2
pdf("../figure/fib/Nfe2l2_zscore_diffusion_map.pdf")
p <- plotEmbedding(ArchRProj = proj,
                       colorBy = "MotifZMatrix",
                       name = "Nfe2l2_101",
                       continuousSet = "horizonExtra",
                       embedding = "palantir",
                       size = 1,
                       quantCut = c(0.01, 0.95))
print(p)
dev.off()
