library(Gviz)
library(glue)
library(cicero)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(genomation)

TimePoint <- "Day3"
TimePoint <- "Healthy"

# load in your data using rtracklayer
#gene_model <- readGFF("../../../Reference/refdata-cellranger-atac-GRCh38-1.1.0/genes/genes.gtf")
gene_model <- readGFF("../../../../ref_cellranger/mm10-1.2/genes.gtf")
gene_model$chromosome <- gene_model$seqid
gene_model$gene <- gene_model$gene_id
gene_model$transcript <- gene_model$transcript_id
gene_model$symbol <- gene_model$gene_name
gene_model$exon <- gene_model$exon_id
gene_model$width <- gene_model$end - gene_model$start + 1
gene_model$feature <- gene_model$transcript_type
gene_model <- subset(gene_model, !is.na(transcript) & !is.na(exon))

df_links <- read.csv(glue("../save/{TimePoint}/detail_corr_pval_final_0.6_0.05_with_annotation.tsv"),
                        header = TRUE, sep = "\t")
df_links$gene <- paste(df_links$gene_chr, 
                        df_links$gene_start, 
                        df_links$gene_start + 1,
                        sep = "_")

df_links$pval_adjust_log10 <- -log10(df_links$padj)

#df_links_sub <- subset(df_links, select = c(peakName, gene, pval_adjust_log10, 
#                                          gene_name, Day10_clusters)) ## Day3
#
#df_links_sub <- subset(df_links, select = c(peakName, gene, pval_adjust_log10, 
#                                          gene_name, Day3_clusters)) ## Day3

df_links_sub <- subset(df_links, select = c(peakName, gene, pval_adjust_log10, 
                                          gene_name, Healthy_clusters)) ## Day3



colnames(df_links_sub) <- c("Peak1", "Peak2", "coaccess", "gene", "cluster")


df_links_fib <- subset(df_links_sub, cluster == "Fibroblasts")


#a <- glue("../../../BigWig/{TimePoint}_Cardiomyocytes.bw")
data_track_cm <- DataTrack(range = as.character(glue("../../../BigWig/{TimePoint}_Cardiomyocytes.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Cardiomyocytes", 
                            col = "#911eb4",
                            ylim = c(0, 4))

data_track_peri <- DataTrack(range = as.character(glue("../../../BigWig/{TimePoint}_Pericytes.vSMC.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Pericytes.vSMC", 
                            col = "#568198",
                            ylim = c(0, 4))

data_track_endo <- DataTrack(range = as.character(glue("../../../BigWig/{TimePoint}_Endothelial cells.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Endothelial", col = "#800000",
                            ylim = c(0, 4))

data_track_mac <- DataTrack(range = as.character(glue("../../../BigWig/{TimePoint}_Macrophages.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Macrophages", col = "#e6194B",
                            ylim = c(0, 4))

data_track_fib <- DataTrack(range = as.character(glue("../../../BigWig/{TimePoint}_Fibroblasts.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Fibroblasts", col = "#f58231",
                            ylim = c(0, 4))

#data_track_3 <- DataTrack(range = "../../ChIPSeq/../BigWig/SRR6426182.bw", 
#                            genome = "mm10", type = "h", 
#                            name = "CM3", col = "#377eb8",
#                          ylim = c(0, 1.5))

df_tf <- read.table(glue("../../../MotifMatch/{TimePoint}_Fibroblasts_mpbs_RUNX2.bed"))
colnames(df_tf) <- c("chromosome", "start", "end", "name", "score", "strand")
df_tf$symbol <- "RUNX2"

tf_track <- GeneRegionTrack(range = df_tf,
                            genome = "mm10",
                            name = "RUNX2",
                            shape = "box",
                            col = "blue",
                            fill = "blue",
                            fontsize = 12,
                            fontcolor = "black")

for(up_gene in c("Col1a1")){
    message(sprintf("plotting links for gene %s", up_gene))
    df_links_sub_fib <- subset(df_links_fib, gene == up_gene)
    #df_links_sub_fib$Peak1 <- str_replace_all(df_links_sub_fib$Peak1, ":", "-")
    df <- as.data.frame(str_split_fixed(df_links_sub_fib$Peak1, "_", 3))
    df$V1 <- as.character(df$V1)
    df$V2 <- as.numeric(as.character(df$V2))
    df$V3 <- as.numeric(as.character(df$V3))
    chr <- unique(df$V1)

    viewpoint <- unique(df_links_sub_fib$Peak2)
    view_start <- stringr::str_split_fixed(viewpoint, "_", 3)[, 2]
    view_end <- stringr::str_split_fixed(viewpoint, "_", 3)[, 3]
    
    minbp <- as.numeric(view_start) - 10000
    maxbp <- as.numeric(view_end) + 10000
    
    gr_peaks <- GRanges(seqnames = df$V1, IRanges(start = df$V2, end = df$V3))
    
    coaccess_cutoff = -log10(0.05)
    comparison_coaccess_cutoff = -log10(0.05)    
    ymax <- max(df_links_sub_fib$coaccess) + 1
    
    link_tracks <- plot_connections(connection_df = df_links_fib, 
                                    gene_model = gene_model,
                                    gene_model_shape = "smallArrow",
                                    collapseTranscripts = "longest",
                     chr = chr, 
                     minbp = minbp, 
                     maxbp = maxbp,
                     alpha_by_coaccess = FALSE,
                     viewpoint = viewpoint,
                     connection_ymax = ymax,
                     coaccess_cutoff = coaccess_cutoff,
                     connection_width = 1,
                     connection_color = "#f58231",
                     include_axis_track = FALSE,
                     return_as_list = TRUE)

    
#    link_track_fib <- link_tracks@trackList[[1]]
#    peak_track_fib <- link_tracks@trackList[[2]]
#    gene_model_track <- link_tracks@trackList[[3]]

    link_track_fib <- link_tracks[[1]]
    peak_track_fib <- link_tracks[[2]]
    gene_model_track <- link_tracks[[3]]
    gene_axis_track <- GenomeAxisTrack(fontsize = 4)


    
    link_track_fib@dp@pars$fontsize <- 12
    data_track_cm@dp@pars$fontsize <- 12
    data_track_peri@dp@pars$fontsize <- 12
    data_track_endo@dp@pars$fontsize <- 12
    data_track_mac@dp@pars$fontsize <- 12
    data_track_fib@dp@pars$fontsize <- 12
    
    
    
    trackList <- list(link_track_fib,
                      peak_track_fib,
                      data_track_fib,
                      data_track_cm,
                      data_track_peri,
                      data_track_endo,
                      data_track_mac,
                      gene_model_track, 
                      tf_track,
                      gene_axis_track)
    

    gr <- GRanges(seqnames = chr, 
                  ranges = IRanges(start = minbp, end = maxbp))
    gr_tf <- GRanges(seqnames = df_tf$chromosome,
                     ranges = IRanges(start = df_tf$start,
                                      end = df_tf$end))
    gr_overlap <- findOverlaps(query = gr_tf, 
                               subject = gr, 
                               type = "within")
    gr_tf_sub <- gr_tf[gr_overlap@from, ]
    
    gr_overlap <- findOverlaps(query = gr_tf_sub, 
                               subject = gr_peaks, 
                               type = "any",
                               maxgap = 1000)
    
    gr_tf_sub <- gr_tf_sub[gr_overlap@from, ]
    
    gr_promoter <- GRanges(seqnames = chr,
                           ranges = IRanges(start = as.numeric(view_start),
                                            end = as.numeric(view_end)))
    if(isEmpty(gr_tf_sub)){
        gr_view <- gr_promoter
    } else{
        gr_tf_sub <- resize(gr_tf_sub, width = 3000, fix = "center")
        gr_view <- c(gr_tf_sub, gr_promoter)
    }
    
    trackList <- HighlightTrack(trackList = trackList, 
                                range = gr_view,
                                chromosome = chr,
                                col = "#F0544F", 
                                fill = "#EFD8D7", 
                                inBackground = FALSE, 
                                alpha = 0.3)
    
    pdf(glue("{TimePoint}_Fib_{up_gene}_V2.pdf"), height = 8, width = 8)
    plotTracks(trackList, title.width = 0.5, showTitle = TRUE, 
               from = minbp, to = maxbp, chromosome = chr, 
               sizes = c(0.4, 0.1, 0.2, 0.2, 0.2, 0.2,0.2,0.2,0.2, 0.2),#0.1, 0.2), 
               transcriptAnnotation = "symbol", background.title = "white", 
               col.border.title = "transparent", lwd.border.title = "transparent", 
               col.axis = "black", fontcolor.legend = "black",
               innerMargin = 3)
    dev.off()
}
