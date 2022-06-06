library(Gviz)
library(glue)
library(cicero)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(genomation)


celltype <- "Macrophages"
tf <- "RUNX1" 
#celltype <- "Cardiomyocytes"
#celltype <- "Pericytes/vSMC"




gene_to_factor = list("Col1a1"= "RUNX2",
                    "Pdgfra" = "CTCF",
                    "Il1b" = c("IRF1", "IRF4", "IRF5", "IRF9", "Mafb", "MAF")) #CTCF CEBPG
#celltype <- "Fibroblasts"
#tf <- "RUNX2" 
#tf <- "CEBPG"



condition_colors <- c("Sham"="#66c2a5", 
                      "Healthy"="#66c2a5", 
                      "Day3"="#fc8d62",
                      "Day10"="#800000")


# load in your data using rtracklayer
#gene_model <- readGFF("../../../Reference/refdata-cellranger-atac-GRCh38-1.1.0/genes/genes.gtf")
gene_model <- readGFF("/data/sz753404/project/ref_cellranger/mm10-1.2/genes.gtf")
gene_model$chromosome <- gene_model$seqid
gene_model$gene <- gene_model$gene_id
gene_model$transcript <- gene_model$transcript_id
gene_model$symbol <- gene_model$gene_name
gene_model$exon <- gene_model$exon_id
gene_model$width <- gene_model$end - gene_model$start + 1
gene_model$feature <- gene_model$transcript_type
gene_model <- subset(gene_model, !is.na(transcript) & !is.na(exon))

TimePoints <- c("Healthy", "Day3", "Day10")

dict <- list(
    "Fibroblasts" = c("Col1a1", "Pdgfra", "Rgs5","Runx1","Bach2", "Smad3", "Fosl1", "Nfkb1"), 
    "Pericytes/vSMC" = c("Col1a1", "Runx1","Bach2", "Smad3", "Fosl1", "Nfkb1"), 
    "Macrophages" = c("Il1b"), 
    "Cardiomyocytes" = c("Mef2c") 
)


list_macro <- c(
    "Macrophages", 
    "inflammatory macrophages", 
    "Anti-inflammatory macrophages"
)



up_genes <- dict[[celltype]] 


df_links_list <- list()
for(TimePoint in TimePoints){
    df_links <- read.csv(glue("../save/{TimePoint}/detail_corr_pval_final_0.6_0.05_with_annotation.tsv"),
                            header = TRUE, sep = "\t")
    df_links$gene <- paste(df_links$gene_chr, 
                            df_links$gene_start, 
                            df_links$gene_start + 1,
                            sep = "_")
    df_links$pval_adjust_log10 <- -log10(df_links$padj)
    
    df_links_sub <- df_links[ c("peakName", "gene", "pval_adjust_log10", 
                                "gene_name", glue("{TimePoint}_clusters"))]

    colnames(df_links_sub) <- c("Peak1", "Peak2", "coaccess", "gene", "cluster")
    df_links_list[[TimePoint]] <- df_links_sub
}


df_links_Healthy <- subset(df_links_list[["Healthy"]], cluster == celltype)
df_links_Day3 <- subset(df_links_list[["Day3"]], cluster == celltype)
df_links_Day10 <- subset(df_links_list[["Day10"]], cluster == celltype)

if(celltype == "Macrophages"){
    df_links_Healthy <- subset(df_links_list[["Healthy"]], cluster %in% list_macro)
    df_links_Day3 <- subset(df_links_list[["Day3"]], cluster %in% list_macro)
    df_links_Day10 <- subset(df_links_list[["Day10"]], cluster %in% list_macro)

}
df_links_day_all_ <- rbind( df_links_Healthy, df_links_Day3, df_links_Day10)


scelltype <- str_replace_all(celltype, "/", ".")
data_track_Healthy <- DataTrack(range = as.character(glue("../../../BigWig/Healthy_{scelltype}.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Healthy", col = condition_colors['Healthy'],
                            ylim = c(0, 4))

data_track_Day3 <- DataTrack(range = as.character(glue("../../../BigWig/Day3_{scelltype}.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Day3", col = condition_colors['Day3'],
                            ylim = c(0, 4))

data_track_Day10 <- DataTrack(range = as.character(glue("../../../BigWig/Day10_{scelltype}.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Day10", col = condition_colors['Day10'],
                            ylim = c(0, 4))


message("loading motifs healthy...", date() )
df1__ <- read.table(glue("../../../MotifMatch/Healthy_{scelltype}_mpbs.bed"))
message("loading motifs day3...", date() )
df2__ <- read.table(glue("../../../MotifMatch/Day3_{scelltype}_mpbs.bed"))
message("loading motifs day10...", date() )
df3__ <- read.table(glue("../../../MotifMatch/Day10_{scelltype}_mpbs.bed"))


for(up_gene in up_genes){
    tfs <- gene_to_factor[[up_gene]]

    for(tf in tfs){

           if(is.na(tf) && celltype == "Fibroblasts"){
               tf <- "RUNX2"
           }
           
           if(is.na(tf) && celltype == "Macrophages"){
               tf <- "RUNX1"
           }


           df1 <- df1__[which(grepl(glue("\\.{tf}$"), df1__$V4)), ]
           df2 <- df2__[which(grepl(glue("\\.{tf}$"), df2__$V4)), ]
           df3 <- df3__[which(grepl(glue("\\.{tf}$"), df3__$V4)), ]
 
           df_tf <- rbind(df1, df2, df3)
           df_tf <- unique(df_tf)
           colnames(df_tf) <- c("chromosome", "start", "end", "name", "score", "strand")
           df_tf$symbol <- glue("{tf}")
           
           write.table(df_tf, file = glue("{tf}.bed"), col.names = FALSE, row.names = FALSE,
                       sep = "\t", quote = FALSE)
           
           tf_track_Day3 <- GeneRegionTrack(range = df_tf,
                                       genome = "mm10",
                                       name = glue("{tf}"),
                                       shape = "box",
                                       col = "blue",
                                       fill = "blue",
                                       fontsize = 12,
                                       fontcolor = "black")
           
           # If I want specific color scheme, I can make a column of color names
           df_links_day_all_$conn_color <- "#800000"
           #df_links_Healthy$conn_color[df_links_Healthy$cluster == "Cardiomyocytes 2"] <- "#9A6324"
           
       
       
           message(sprintf("plotting links for gene %s", up_gene))
           df_links_day_all <- subset(df_links_day_all_, gene == up_gene)
           if(nrow(df_links_day_all) == 0){
              message(glue("no link: {up_gene}"))
              next 
           }
           df <- as.data.frame(str_split_fixed(df_links_day_all$Peak1, "_", 3))
           df$V1 <- as.character(df$V1)
           df$V2 <- as.numeric(as.character(df$V2))
           df$V3 <- as.numeric(as.character(df$V3))
           chr <- unique(df$V1)
           minbp <- min(df$V2) - 10000
           maxbp <- max(df$V3) + 10000
           viewpoint <- unique(df_links_day_all$Peak2)
           
           gr_peaks <- GRanges(seqnames = df$V1, IRanges(start = df$V2, end = df$V3))
           
           
           coaccess_cutoff = -log10(0.05)
           comparison_coaccess_cutoff = -log10(0.05)    
           ymax <- max(df_links_day_all$coaccess) + 1
           
           link_tracks <- plot_connections(connection_df = df_links_day_all, 
                                           gene_model = gene_model,
                                           gene_model_shape = "smallArrow",
                                           collapseTranscripts = "longest",
                            chr = chr, 
                            minbp = minbp, 
                            maxbp = maxbp,
                            comparison_track = df_links_Day10,
                            alpha_by_coaccess = FALSE,
                            # viewpoint = viewpoint,
                            connection_ymax = ymax,
                            comparison_ymax = ymax,
                            coaccess_cutoff = coaccess_cutoff,
                            comparison_coaccess_cutoff = comparison_coaccess_cutoff,
                            connection_width = 1,
                            comparison_connection_width = 1,
                            connection_color = "conn_color",
                            comparison_connection_color = "#9A6324",
                            include_axis_track = FALSE,
                            return_as_list = TRUE)
       
           
           gene_model_track <- link_tracks[[5]]
       
       
           df_links_Day3$conn_color <- condition_colors["Day3"]
           link_tracks <- plot_connections(connection_df = df_links_Day3, 
                                           gene_model = gene_model,
                                           gene_model_shape = "smallArrow",
                                           collapseTranscripts = "longest",
                            chr = chr, 
                            minbp = minbp, 
                            maxbp = maxbp,
                            comparison_track = df_links_Day10,
                            alpha_by_coaccess = FALSE,
                            # viewpoint = viewpoint,
                            connection_ymax = ymax,
                            comparison_ymax = ymax,
                            coaccess_cutoff = coaccess_cutoff,
                            comparison_coaccess_cutoff = comparison_coaccess_cutoff,
                            connection_width = 1,
                            comparison_connection_width = 1,
                            connection_color = "conn_color",
                            comparison_connection_color = condition_colors['Day3'],
                            include_axis_track = FALSE,
                            return_as_list = TRUE)
           link_track_Day3 <- link_tracks[[1]]
           peak_track_Day3 <- link_tracks[[2]]
       
       
           df_links_Day10$conn_color <- condition_colors["Day10"]
           link_tracks <- plot_connections(connection_df = df_links_Day10, 
                                           gene_model = gene_model,
                                           gene_model_shape = "smallArrow",
                                           collapseTranscripts = "longest",
                            chr = chr, 
                            minbp = minbp, 
                            maxbp = maxbp,
                            comparison_track = df_links_Day10,
                            alpha_by_coaccess = FALSE,
                            # viewpoint = viewpoint,
                            connection_ymax = ymax,
                            comparison_ymax = ymax,
                            coaccess_cutoff = coaccess_cutoff,
                            comparison_coaccess_cutoff = comparison_coaccess_cutoff,
                            connection_width = 1,
                            comparison_connection_width = 1,
                            connection_color = "conn_color",
                            comparison_connection_color = condition_colors['Day10'],
                            include_axis_track = FALSE,
                            return_as_list = TRUE)
           link_track_Day10 <- link_tracks[[1]]
           peak_track_Day10 <- link_tracks[[2]]
       
       
       
       
           gene_axis_track <- GenomeAxisTrack(fontsize = 4)
           
           link_track_Day3@dp@pars$fontsize <- 12
           link_track_Day10@dp@pars$fontsize <- 12
           data_track_Day3@dp@pars$fontsize <- 12
           data_track_Day10@dp@pars$fontsize <- 12
           
           trackList <- list(link_track_Day3,
                             peak_track_Day3,
                             link_track_Day10,
                             peak_track_Day10,
                             data_track_Healthy,
                             data_track_Day3,
                             data_track_Day10,
                             gene_model_track, 
                             tf_track_Day3,
                             gene_axis_track)
           
           view_start <- stringr::str_split_fixed(viewpoint, "_", 3)[, 2]
           view_end <- stringr::str_split_fixed(viewpoint, "_", 3)[, 3]
        
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
           
           gr_tf_sub <- resize(gr_tf_sub, width = 3000, fix = "center")
           gr_promoter <- GRanges(seqnames = chr,
                            ranges = IRanges(start = as.numeric(view_start),
                                             end = as.numeric(view_end)))
           
           gr_view <- c(gr_tf_sub, gr_promoter)
           
           trackList <- HighlightTrack(trackList = trackList, 
                                       range = gr_view,
                                       chromosome = chr,
                                       col = "#F0544F", 
                                       fill = "#EFD8D7", 
                                       inBackground = FALSE, 
                                       alpha = 0.3)
           
           pdf(glue("Macrophages_tfs/condition_{celltype}_{up_gene}_{tf}.pdf"), height = 8, width = 8)
           plotTracks(trackList, title.width = 0.5, showTitle = TRUE, 
                      from = minbp, to = maxbp, chromosome = chr, 
                      sizes = c(0.4, 0.1, 0.4, 0.1, 0.3, 0.3, 0.3, 0.3, 0.2, 0.05), 
                      transcriptAnnotation = "symbol", background.title = "white", 
                      col.border.title = "transparent", lwd.border.title = "transparent", 
                      col.axis = "black", fontcolor.legend = "black",
                      innerMargin = 3)
           dev.off()
    }   
}
