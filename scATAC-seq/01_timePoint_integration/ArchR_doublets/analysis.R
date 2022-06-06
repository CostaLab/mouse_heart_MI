library(ggplot2)
library(stringr)
library(magrittr)
library(WriteXLS)
library(tidyr)
library(dplyr)
library(ArchR)


## set working directory
setwd("data")


## Creating Arrow Files

addArchRThreads(threads = 12)
addArchRGenome("mm10")

inputFiles <- c(paste0(dirIn, "/DS1/outs/fragments.tsv.gz"),
                paste0(dirIn, "/DS3/outs/fragments.tsv.gz"),
                paste0(dirIn, "/DS4/outs/fragments.tsv.gz"),
                paste0(dirIn, "/DS5/outs/fragments.tsv.gz"),
                paste0(dirIn, "/DS6/outs/fragments.tsv.gz"),
                paste0(dirIn, "/DS7/outs/fragments.tsv.gz"))

sampleNames <- c("sham_female", "day3_female", "day3_male",
                 "day10_female", "day10_male", "sham_male")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 5, 
  filterFrags = 3000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


## Inferring Doublets

ArrowFiles <- c(
	"day10_female.arrow",
	"day10_male.arrow",
	"day3_female.arrow",
	"day3_male.arrow",
	"sham_female.arrow",
	"sham_male.arrow")


doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)



## Creating an ArchRProject

# With our Arrow files in hand, we are now ready to create an ArchRProject. An ArchRProject is associated with a set of Arrow files and is the backbone of nearly all ArchR analyses.
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchROutput",
  showLogo = FALSE,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

# Now we can filter putative doublets based on the previously determined doublet scores using the filterDoublets() function. This doesn’t physically remove data from the Arrow files but rather tells the ArchRProject to ignore these cells for downstream analysis.


cellColData <- proj@cellColData
saveRDS(cellColData, "cellColData_before_rm_doublets.Rds")

proj <- filterDoublets(ArchRProj = proj)

cellColData <- proj@cellColData
saveRDS(cellColData, "cellColData_after_rm_doublets.Rds")


## Dimensionality Reduction and Clustering

# ArchR implements an iterative LSI dimensionality reduction via the addIterativeLSI() function.
proj <- addIterativeLSI(ArchRProj = proj, 
                        useMatrix = "TileMatrix", 
                        name = "IterativeLSI",
                        varFeatures = 200000)

proj <- addClusters(input = proj, 
                    reducedDims = "IterativeLSI",
                    resolution = 0.6)

table(proj$Clusters)

## Visualizing in a 2D UMAP Embedding

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", 
                    name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", 
                    name = "Clusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")


sessionInfo()

