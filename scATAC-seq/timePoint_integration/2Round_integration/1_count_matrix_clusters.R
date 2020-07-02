mute <- suppressPackageStartupMessages
mute(library(ggplot2))
mute(library(stringr))
mute(library(magrittr))
mute(library(WriteXLS))
mute(library(tidyr))
mute(library(dplyr))
mute(library(plotly))
mute(library(Signac))
mute(library(Seurat))
mute(library(cluster))
mute(library(clustree))
mute(library(mclust))
mute(library(cowplot))
mute(library(gridExtra))
mute(library(ggrastr))
mute(library(viridis))
mute(library(GenomicRanges))
mute(library(GenomeInfoDb))
mute(library(BSgenome.Mmusculus.UCSC.mm10))
mute(library(EnsDb.Hsapiens.v86))
mute(library(data.table))
mute(library(patchwork))
mute(library(foreach))
mute(library(doParallel))
mute(library(Matrix))
mute(library(optparse))

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

atac_dir <- sprintf("../../Seurat3/1Round_IntegrationWithscOpen_timePoint/%s/", sample)
heart.atac <- readRDS(file.path(atac_dir, "obj.umap_integrated.Rds"))

file_path = sprintf("../../1Round_Peakcalling/data/%s/Fragments", sample)
fragment.path = file.path(file_path, "fragments.tsv")


## peak calling for clusters 

Idents(heart.atac) <- "seurat_clusters"
#heart.atac <- RenameIdents(heart.atac, cell.names)
heart.atac$seurat_clusters <- Idents(heart.atac)

print(unique(heart.atac$seurat_clusters))

fragment_dir <- paste0(dir.out, "/Fragments")
peaks_dir <- paste0(dir.out, "/Peaks")
if(!dir.exists(fragment_dir)){
    dir.create(fragment_dir)
}

if(!dir.exists(peaks_dir)){
    dir.create(peaks_dir)
}

for(a_cluster in unique(heart.atac$seurat_clusters)) {
    fragment_file_filtered <- paste0(fragment_dir, "/", a_cluster, ".tsv")
    if(file.exists(fragment_file_filtered)){
        file.remove(fragment_file_filtered)
    }
    
    cells <- colnames(subset(heart.atac, subset = seurat_clusters == a_cluster))
    FilterFragments(fragment.path = fragment.path,
                    cells = cells,
                    output.path = fragment_file_filtered)
    
    command <- sprintf("macs2 callpeak -g mm --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all", 
	    a_cluster, 
	    paste0(fragment_file_filtered, ".bgz"), 
	    peaks_dir)
    message("Running Macs2...")
	message(command)
	system(command, intern = TRUE)
}



## define helper function

readSummits <- function(file){
    df <- data.frame(readr::read_tsv(file,
                                     col_names = c("chr","start","end","name","score")))
    df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
    return(makeGRangesFromDataFrame(df = df,
                                    keep.extra.columns = TRUE,
                                    starts.in.df.are.0based = TRUE))
}

clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE){
      gr <- sort(sortSeqlevels(gr))
      r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
      o <- findOverlaps(gr,r)
      mcols(gr)$cluster <- subjectHits(o)
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      gr <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(gr))
      mcols(gr)$cluster <- NULL
      return(gr)
}

nonOverlappingGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
    stopifnot(by %in% colnames(mcols(gr)))
    i <-  0
    gr_converge <- gr
    while(length(gr_converge) > 0){
      if(verbose){
        message(".", appendLF = FALSE)
      }
      i <-  i + 1
      gr_selected <- clusterGRanges(gr = gr_converge, filter = TRUE, by = by, decreasing = decreasing)
      gr_converge <- subsetByOverlaps(gr_converge ,gr_selected, invert=TRUE) #blacklist selected gr
      if(i == 1){ #if i=1 then set gr_all to clustered
        gr_all <- gr_selected
      }else{
        gr_all <- c(gr_all, gr_selected)
      } 
    }
    if(verbose){
      message("\nSelected ", length(gr_all), " from ", length(gr))
    }
    gr_all <- sort(sortSeqlevels(gr_all))
    return(gr_all)

}



## Make Non-Overlapping Peak Set

blacklist_file = "../../Blacklists/mm10-blacklist.v2.bed.gz"

blacklist <- rtracklayer::import.bed(blacklist_file)

chromSizes <- GRanges(names(seqlengths(BSgenome.Mmusculus.UCSC.mm10)), 
                      IRanges(1, seqlengths(BSgenome.Mmusculus.UCSC.mm10)))
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes,
                                                    pruning.mode = "coarse")
      

peaks_files <- list.files(peaks_dir, 
                          pattern = "\\_summits.bed", 
                          full.names = TRUE)

gr_list <- GenomicRangesList(lapply(peaks_files, function(x){
        extended_summits <- readSummits(x) %>%
            resize(., width = 2 * 250 + 1, fix = "center") %>%
            subsetByOverlaps(.,chromSizes,type="within") %>%
              subsetByOverlaps(.,blacklist,invert=TRUE) %>%
              nonOverlappingGRanges(., by="score", decreasing=TRUE)
        extended_summits <- extended_summits[order(extended_summits$score, 
                                                   decreasing=TRUE)]
        extended_summits <- head(extended_summits, 100000)
        mcols(extended_summits)$scoreQuantile <-trunc(rank(mcols(extended_summits)$score)) / length(mcols(extended_summits)$score)
        extended_summits
        }))


unionPeaks <- nonOverlappingGRanges(unlist(gr_list), 
                                    by = "scoreQuantile", 
                                    decreasing = TRUE)
unionPeaks <- sort(sortSeqlevels(unionPeaks))

unionPeaks <- unionPeaks[seqnames(unionPeaks) %in% paste0("chr",c(1:19,"X"))]
unionPeaks <- keepSeqlevels(unionPeaks, paste0("chr",c(1:19,"X")))

df <- data.frame(seqnames=seqnames(unionPeaks),
                 starts=start(unionPeaks)-1,
                 ends=end(unionPeaks),
                 name=paste0("peaks_", 1:length(unionPeaks)),
                 score=score(unionPeaks))

write.table(df, file = paste0(dir.out, "/unionPeaks.bed"), 
            sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)



## Session information

sessionInfo()

