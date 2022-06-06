library(dplyr)
library(optparse)
library(glue)
library(foreach)
library(doParallel)
registerDoParallel(cores=3)



#2*pnorm(-abs(((o$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))

#$FDR <- p.adjust(o$Pval, method = "fdr")

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




message("processing ", TimePoint, date())
df <- read.csv(file=glue("save/{TimePoint}/corr.tsv"), sep="\t", stringsAsFactors = F)
fn = glue("save/{TimePoint}/null_hypothesis.Rds")
null_dist_list <- readRDS(file=fn)
stt_vars <- lapply(null_dist_list, function(lst) return(c("mean"=mean(lst), "sd"=sd(lst))))
df$pval <- sapply(1:nrow(df), function(x){
          chromosome <- strsplit(df[x,]$peak, "_")[[1]][1]
          mean_null_corr <- stt_vars[[chromosome]]["mean"]
          sd_null_corr <- stt_vars[[chromosome]]["sd"]
          2*pnorm(-abs(((df[x,]$corr - mean_null_corr) / sd_null_corr)))
  })  

df$padj <- p.adjust(df$pval, method="BH")

write.table(df, file=glue("save/{TimePoint}/corr_pval.tsv"), sep = "\t", quote = F, row.names = F)

message("finished", TimePoint, date())

