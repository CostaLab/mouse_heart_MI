#```{r}
mute <- suppressPackageStartupMessages
mute(library(cicero))
mute(library(Seurat))
mute(library(data.table))
mute(library(Matrix))
mute(library(GenomicRanges))
mute(library(magrittr))
mute(library(SummarizedExperiment))
mute(library(optparse))
mute(library(Signac))
mute(library(yaml))
mute(library(Rcpp))
mute(library(stringr))
mute(library(WriteXLS))
mute(library(glue))



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






#gtf_file = "/mnt/data/data_files/Reference/refdata-cellranger-atac-mm10-1.2.0/genes/genes.gtf"
gtf_file = "../../../ref_cellranger/mm10-1.2/genes.gtf"

atac_file <- glue("save/{TimePoint}/atac_with_promoter.Rds")
atac <- readRDS(file=atac_file)


set.seed(1)

getGeneGTF <- function(file){
  #Import
  message("Reading in GTF...")
  importGTF <- rtracklayer::import(file)
  #Exon Info
  message("Computing Effective Exon Lengths...")
  exonGTF <- importGTF[importGTF$type=="exon",]
  exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
  exonReduced <- unlist(exonList, use.names=TRUE)
  mcols(exonReduced)$gene_id <- names(exonReduced)
  mcols(exonReduced)$widths <- width(exonReduced)
  exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
  exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>% 
    unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
  #Gene Info
  message("Constructing gene GTF...")
  geneGTF1 <- importGTF[importGTF$type=="gene",]
  geneGTF2 <- GRanges(
      #seqnames=#paste0("chr",seqnames(geneGTF1)),
      seqnames=seqnames(geneGTF1),
      ranges=ranges(geneGTF1),
      strand=strand(geneGTF1),
      gene_name=geneGTF1$gene_name,
      gene_id=geneGTF1$gene_id
    #) %>% keepFilteredChromosomes %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
    ) %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
  mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
  return(geneGTF2)
}


fixA <- "center"
fixB <- "start"
associationWindow <- 2 * 250*10^3 + 1 #+-250 Kb
corCutOff <- 0.35 #Pearson Correlation Cutoff
fdrCutOff <- 0.1 #FDR Cutoff
distCutOff <- 2500 #Min Dist to TSS


#load(file="save/atac_ZKY095.Rdata")
#load(file="save/atac_a_ZKY095.Rdata")


gtf <- getGeneGTF(gtf_file)

DefaultAssay(atac) <- "promoter"
counts <- GetAssayData(atac, slot="counts", assay="promoter")
inter_gene_name <- intersect(rownames(atac), gtf@elementMetadata$gene_name)
counts <- counts[inter_gene_name, ] 

seRNA <- SummarizedExperiment(
  assays = SimpleList(counts= counts)
)

#rm(atac)
#gc()


gtfMatch <- gtf[na.omit(match(rownames(seRNA), gtf$gene_name))]
names(gtfMatch) <- rownames(seRNA)
rowRanges(seRNA) <- gtfMatch



DefaultAssay(atac) <- "peaks"

atac_mtx <- GetAssayData(atac, slot = "counts", assay = "peaks")
#atac_mtx <- readRDS(file="../../2Round_Peakcalling/data/all_SharedPeaks_matrix.Rds")
#atac_mtx <- atac_mtx[, colnames(atac)]
#atac[["raw_peaks"]] <- CreateAssayObject(counts=atac_mtx)
#saveRDS(atac, "save/atac_with_promoter.Rds")


peaks <- rownames(atac_mtx)
peak.ranges <- StringToGRanges(regions =peaks , sep = c("-", "-"))
seATAC <- SummarizedExperiment(
                        assays <- SimpleList(atac_mtx),
                        colData = atac@meta.data,
                        rowRanges = peak.ranges)


#rm(atac_mtx)
#gc()
#assay(seATAC) <- log2(edgeR::cpm(assay(seATAC))/100+1)
#assay(seRNA) <- log2(edgeR::cpm(assay(seRNA))/100+1)


seRNAWindow <- resize(rowRanges(seRNA), width = 1, fix = fixB) %>%
  {suppressWarnings(resize(., width = associationWindow, fix = "center"))} %>% trim(.)

#Keep only seATAC within association window
seATAC <- seATAC[unique(queryHits(findOverlaps(resize(seATAC,1,fixA), seRNAWindow, ignore.strand = TRUE)))]

#Getting distances
message("Getting Distances...")
o <- findOverlaps(seRNAWindow, resize(rowRanges(seATAC),1,fixA), ignore.strand = TRUE)


#Get Distance from Fixed point A B correct for minus stranded
mcols(o)$distance <- start(resize(rowRanges(seATAC),1,fixA))[subjectHits(o)] - start(resize(rowRanges(seRNA),1,fixB))[queryHits(o)]
mcols(o)$distance[which(as.character(strand(rowRanges(seRNA)))[queryHits(o)]=="-")] <- -1*mcols(o)$distance[which(as.character(strand(rowRanges(seRNA)))[queryHits(o)]=="-")]

#Add other info
o <- DataFrame(o)
colnames(o) <- c("B","A","distance")
o <- o[,c("A","B","distance")]

#```




#```{r}
#Get GTF
#gtf <- getGeneGTF(gtf_file)
tssRNA <- resize(gtf, 1, "start")
strand(tssRNA) <- "*" 
peakLinks <- rowRanges(seATAC)[o[,1]]
geneLinks <- rowRanges(seRNA) %>% resize(1, "start") %>% {.[o[,2]]}

mcolsLinks <- data.frame(geneLinks)[,c("seqnames","start","strand","gene_name","gene_id","exonLength")]
colnames(mcolsLinks) <- c("gene_chr","gene_start","gene_strand","gene_name","gene_id","exonLength")
mcolsLinks <- cbind(mcolsLinks, data.frame(o))
mcolsLinks$nearestGene <- tssRNA$gene_name[subjectHits(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))]
mcolsLinks$nearestGeneDistance <- mcols(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))$distance
mcols(peakLinks) <- mcolsLinks
peakLinks$peakName <- paste(seqnames(peakLinks), start(peakLinks), end(peakLinks), sep = "_")

rm(seRNA, seATAC)
gc()


#peakName geneName peak_chr peak_start peak_end gene_chr gene_start gene_strand gene_name gene_id exonLength distance
out_list <- peakLinks@elementMetadata@listData
out_list$geneName <- out_list$gene_name
splited_list <- str_split(out_list$peakName, "_")
out_list$peak_chr <- sapply(splited_list, function(x) x[1]) 
out_list$peak_start <- sapply(splited_list, function(x) x[2])
out_list$peak_end <- sapply(splited_list, function(x) x[3])
#out_list$peakName <- sub("_", ":",  out_list$peakName)
#out_list$peakName <- sub("_", "-",  out_list$peakName)


out_list <- out_list[c("peakName", "geneName", "peak_chr",
                          "peak_start", "peak_end", "gene_chr",
                          "gene_start", "gene_strand", "gene_name",
                          "gene_id", "exonLength", "distance")]

#peak    gene    distance    peak_type
tsv_list <- out_list[c("peakName", "geneName", "distance")]
tsv_df <- as.data.frame(tsv_list)
tsv_df <- rename(tsv_df, c("peakName" = "peak",
                           "geneName" = "gene"))
write.table(tsv_df, file=glue("save/{TimePoint}/peak2gene_putative-simple.tsv"), sep="\t", quote=F, row.names = F)

write.table(as.data.frame(out_list), file=glue("save/{TimePoint}/peak2gene_putative.tsv"), sep="\t", quote=F, row.names = F)



#WriteXLS(as.data.frame(out_list),
#         file.path("data", "peak2gene_putative.xlsx"),
#         SheetNames = "peak2gene_putative")



#```
