mute <- suppressPackageStartupMessages
mute(library(dplyr))
mute(library(glue))
#mute(library(VennDiagram))
mute(library(tidyverse))
mute(library(optparse))
#futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")



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


filter_pval <- 0.05

df <- read.csv(file=glue("save/{TimePoint}/corr_pval.tsv"), sep = "\t")

fdf <- df %>% filter(corr>0)    
fdf <- fdf %>% filter(abs(distance) > 2000)
fdf$Sig <- ifelse(fdf$padj<0.05, "TRUE", "FALSE")


for(x in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)){

    fdf <- fdf %>% filter(corr > x)
    write.table(fdf, file=glue("save/{TimePoint}/corr_pval_final_{x}.tsv"),
                        sep = "\t", quote = F, row.names = F)
    
    write.table(fdf %>% filter(Sig == TRUE), file=glue("save/{TimePoint}/corr_pval_final_{x}_0.05.tsv"),
                        sep = "\t", quote = F, row.names = F)


}



#write.table(fdf, file=glue("save/{TimePoint}/corr_pval_final.tsv"),

#                    sep = "\t", quote = F, row.names = F)

fdf <- fdf %>% filter(fdf$padj<filter_pval)
write.table(fdf, file=glue("save/{TimePoint}/corr_pval_final_{filter_pval}.tsv"),
                    sep = "\t", quote = F, row.names = F)


#write.table(fdf, file=glue("save/{TimePoint}/corr_pval_final.tsv"),
#                    sep = "\t", quote = F, row.names = F)

fdf <- fdf %>% filter(fdf$padj<filter_pval)
write.table(fdf, file=glue("save/{TimePoint}/corr_pval_final_{filter_pval}.tsv"),
                    sep = "\t", quote = F, row.names = F)


