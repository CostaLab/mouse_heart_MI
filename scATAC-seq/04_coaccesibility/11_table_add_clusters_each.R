suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(glue))
source("config_heatmap_orders.R")

parser <- OptionParser()
parser <- add_option(parser, c("-p", "--plots"), type="character", default="",
                    help="plots for different pval_option [default %default]",
                    metavar="character")

pa = parse_args(parser)

plots_opt <-  ifelse(pa$plots == "",  "",  paste0("_", pa$plots))
print(plots_opt)
dir.create(sprintf("plots%s", plots_opt))

condition_vec <- c("Healthy", "Day3", "Day10")  


pb = txtProgressBar(min = 0, max = length(condition_vec), initial = 0) 
for(a_condition in condition_vec){
	a_name = glue("save/{a_condition}/detail_corr_pval_final{plots_opt}.tsv")
	print(a_name)
	df <- read.csv(file=a_name , sep="\t")
    df[, paste0(a_condition, "_clusters")] <- NA

	celltype.list <- get_cluster_anno(plots_opt, a_condition)
	number.clusters <- length(celltype.list)

    fname <- glue("plots{plots_opt}/merge_info_{a_condition}_{number.clusters}.Rds")
    c_list <- readRDS(fname)
    setTxtProgressBar(pb, which(condition_vec == a_condition) )
    message(a_condition, "\t",date())
    for(nm in names(c_list)){
        a_list <- c_list[[nm]]
        nums <- length(a_list[[1]])
        message(nm, date())
        for(i in 1:nums){
           peak <- a_list[[1]][i]
           gene <- a_list[[2]][i]
		   celltype = celltype.list[as.character(nm)]
           df[df$peakName==peak & df$geneName == gene,][, paste0(a_condition, "_clusters")] <-celltype
        }
    }
    fname <- glue("save/{a_condition}/detail_corr_pval_final{plots_opt}_with_annotation.tsv")
    write.table(df, file = fname, sep="\t", row.names=F,quote=FALSE)
}



