library(glue)
library(dplyr)
library(tidyverse)
library(WriteXLS)



TimePoint="Day3"
plots_opt = "_0.6_0.05"

TimePoints <- c("Healthy", "Day3", "Day10")




dict <- c("Myofibroblasts" = "Fibroblasts",
         "Fibroblasts" = "Fibroblasts",
         "Pericytes.vSMC"= "Pericytes.vSMC",
         "Anti-inflammatory macrophages" = "Macrophages",
         "inflammatory macrophages" = "Macrophages",
         "Endothelial cells" = "Endothelial cells",
         "Cardiomyocytes" = "Cardiomyocytes")





list_df <- list()

for( TimePoint in TimePoints){

    f = glue("save/{TimePoint}/detail_corr_pval_final{plots_opt}_with_annotation.tsv")
    
    df <- read.csv(f, sep="\t")
    
    df.splited <- split(df, df[, glue('{TimePoint}_clusters')])
    
    
    
    
    df.splited$Fibroblasts %>% group_by(geneName) %>% summarise(counts= length(geneName)) 
    
    
    df_statistic <- lapply(df.splited, function(df) df %>% group_by(geneName) %>% summarise(counts= length(geneName)))
    
    df_statistic <- lapply(df_statistic, function(x) x[order(x$counts, decreasing=T), ])
    
    names(df_statistic) <- stringr::str_replace_all(names(df_statistic), "/", ".")
    list_df[[TimePoint]] <- df_statistic
    
    WriteXLS(
               df_statistic,
               glue("save/{TimePoint}/links_statistics_{TimePoint}{plots_opt}.xlsx"),
               SheetNames = names(df_statistic))
    
}

list_ct <- list()
for(TimePoint in TimePoints){
    for(ct in names(dict)){
        if(ct %in% names(list_df[[TimePoint]])) {
            if(is.null(list_ct[[glue("{dict[ct]}")]])){
                list_ct[[glue("{dict[ct]}")]] <- list()
             }
             list_ct[[glue("{dict[ct]}")]][[glue("{TimePoint}_{ct}")]] <- list_df[[TimePoint]][[ct]]
        }     

    } 
}



full_join_by_gene_name = function(x, y) full_join(x, y, by="geneName")


xlsx_df_list <- list()

for(nm in names(list_ct)){
    xlsx_df_list[[nm]] <- Reduce(full_join_by_gene_name, list_ct[[nm]]) 
    colnames(xlsx_df_list[[nm]]) <- c("geneName", names(list_ct[[nm]]))
}


WriteXLS(xlsx_df_list, 
         glue("save/links_statistics{plots_opt}.xlsx"),
         SheetNames = names(xlsx_df_list))

filtered_xlsx_df_list <- lapply(xlsx_df_list, function(x) x[complete.cases(x), ])
WriteXLS(filtered_xlsx_df_list, 
         glue("save/links_statistics_completecases_{plots_opt}.xlsx"),
         SheetNames = names(filtered_xlsx_df_list))




Reduced_all <- Reduce(full_join_by_gene_name, xlsx_df_list) 
colnames(Reduced_all) <- c("geneName", sapply(names(xlsx_df_list), function(x) names(xlsx_df_list[[x]])[-1]))


WriteXLS(Reduced_all, 
         glue("save/links_statistics_all_{plots_opt}.xlsx"),
         SheetNames = "ALL_LINKS")





