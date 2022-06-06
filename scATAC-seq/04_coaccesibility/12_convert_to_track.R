library(stringr)
library(glue)
library(data.table)

TimePoint_list <- c("Healthy", "Day3", "Day10")

for(TimePoint in TimePoint_list){
    filename <- glue("save/{TimePoint}/detail_corr_pval_final_0.6_0.05_with_annotation.tsv")
    df <- read.csv(file = filename, header = TRUE, sep = "\t")
    df$peak_center <- (df$peak_start + df$peak_end) / 2
    df_links <- df[, c(3, 18, 7, 17, 13)]
    colnames(df_links) <- c("peak_chr", "peak_center", "gene_start", "celltype", "corr")
    df_links$celltype <- str_replace_all(df_links$celltype, " ", "_")
    header <- sprintf("track name=\"%s\" description=\"%s\" graphType=junctions\n", TimePoint, TimePoint)
    cat(header, file = glue("Links_{TimePoint}.bed"))
    fwrite(x = df_links, file = glue("tracks/Links_{TimePoint}.bed"), sep = "\t", 
           row.names = FALSE, col.names = FALSE, append = TRUE)
}
