library(dplyr)
library(glue)



plot_opt = "_0.6_0.05"



TimePoint = "Healthy"
f <- glue("save/{TimePoint}/detail_corr_pval_final{plot_opt}_with_annotation.tsv")
healthy_df <- read.csv(f, sep="\t")
healthy_links_df <- healthy_df[, c("peakName", "geneName", "Healthy_clusters")] 
healthy_links = sapply(1:nrow(healthy_links_df), function(x) paste(as.character(healthy_links_df[x, ]), collapse = ",")) 

TimePoint = "Day3"
f <- glue("save/{TimePoint}/detail_corr_pval_final{plot_opt}_with_annotation.tsv")
day3_df <- read.csv(f, sep="\t")
day3_links_df <- day3_df[, c("peakName", "geneName", "Day3_clusters")] 
day3_links = sapply(1:nrow(day3_links_df), function(x) paste(as.character(day3_links_df[x, ]), collapse = ",")) 



TimePoint = "Day10"
f <- glue("save/{TimePoint}/detail_corr_pval_final{plot_opt}_with_annotation.tsv")
day10_df <- read.csv(f, sep="\t")
day10_links_df <- day10_df[, c("peakName", "geneName", "Day10_clusters")] 
day10_links = sapply(1:nrow(day10_links_df), function(x) paste(as.character(day10_links_df[x, ]), collapse = ",")) 

Healthy_only <- healthy_links %>% setdiff( day3_links) %>% setdiff(day10_links) %>% list
Day3_only <- day3_links  %>% setdiff(healthy_links) %>% setdiff(day10_links) %>% list
Day10_only <- day10_links  %>% setdiff(healthy_links) %>% setdiff(day3_links) %>% list 

shared_links <- Reduce(intersect, list(healthy_links, day3_links, day10_links)) %>% list




