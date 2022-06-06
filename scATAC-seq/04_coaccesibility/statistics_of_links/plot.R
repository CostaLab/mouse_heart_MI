library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(ggrepel)
library(tidyr)
library(readxl)

#df_deg <- read.table("../../scRNA_filtered/DGEA/CK160.CM1.vs.CM2.dge.txt", header = TRUE)

rna_dir <- "../../../Mouse_scATAC_scRNAPipeLine/"
fname = "charts/Sham.vs.Day3.de_celltype.xlsx"
fname = "save/de_cond_cluster.Rds"
#workbook <- loadWorkbook("mmc5.xlsx")

#workbook <- loadWorkbook("mmc5.xlsx")

#x <- excel_sheets(file.path(rna_dir, fname))
#lst <- list()
#for(nm in x){ 
#  lst[[nm]] <- read_excel(file.path(rna_dir,  fname), sheet=nm)
#}
#


vsRest_lists <- readRDS(file=file.path(rna_dir, "save", "all_de_list.Rds"))
cond_lists <- readRDS(file = file.path(rna_dir, "save", "de_cond.Rds")) 
cond_cluster_lists <- readRDS(file = file.path(rna_dir, "save", "de_cond_cluster.Rds")) 

df_deg <- cond_cluster_lists[["Day3"]][['Fibroblasts.vs.Pericytes/vSMC']]


df_links <- read.csv("../save/Day3/detail_corr_pval_final_0.6_0.05_with_annotation.tsv",
                     header = TRUE, sep = "\t") %>%
    subset(Day3_clusters %in% c("Fibroblasts", "Pericytes/vSMC")) %>%
    group_by(geneName, Day3_clusters) %>% summarise(num_links = n()) %>%
    spread(Day3_clusters, num_links)

colnames(df_links) <- c("gene", "num_links_fib", "num_links_pericytes")
df_links[is.na(df_links)] <- 0

df_links$de_gene <- ifelse(df_links$gene %in% df_deg$gene, "YES", "NO")
df_links <- merge.data.frame(df_links, df_deg, by = "gene", all = TRUE)
df_links$CellType <- ifelse(df_links$avg_logFC > 0, "fib", "pericytes")

df_links <- df_links[!is.na(df_links$num_links_fib) & 
                         !is.na(df_links$num_links_pericytes), ]

df_links$diff_num_links <- df_links$num_links_fib - df_links$num_links_pericytes 


df_text <- subset(df_links, abs(diff_num_links) > 10 & de_gene == "YES")

p <- ggplot(data = df_links, aes(x = reorder(gene, diff_num_links), y = diff_num_links, group=1)) +
    geom_point(size = 1) +
    geom_text_repel(data = df_text, aes(x = reorder(gene, diff_num_links), y = diff_num_links, 
                                        label = gene, color = CellType)) +
    xlab("Ranked gene list") + 
    ylab("Difference of number of correlated peaks") +
    scale_color_manual(values = c("fib" = "#800000", 
                                  "pericytes" = "#9A6324")) +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank()) +
    guides(
        fill = guide_legend(
            title = "Legend Title",
            override.aes = aes(label = "")
        )
    )

pdf("gene_ranked.pdf", height = 8, width = 8)
print(p)
dev.off()

write.table(df_links, file = "gene_ranked.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)





p <- ggplot(data = df_links, aes(x = reorder(gene, diff_num_links), y = diff_num_links, group=1)) +
    geom_point(size = 1) +
    #geom_text_repel(data = df_text, aes(x = reorder(gene, diff_num_links), y = diff_num_links, 
    #                                    label = gene, color = CellType)) +
    xlab("Ranked gene list") + 
    ylab("Difference of number of correlated peaks") +
    scale_color_manual(values = c("fib" = "#800000", 
                                  "pericytes" = "#9A6324")) +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank()) +
    guides(
        fill = guide_legend(
            title = "Legend Title",
            override.aes = aes(label = "")
        )
    )

pdf("gene_ranked_noname.pdf", height = 8, width = 8)
print(p)
dev.off()

#write.table(df_links, file = "gene_ranked.txt",
#            sep = "\t", quote = FALSE, row.names = FALSE)
