library(stringr)
library(glue)

file="/data/scATA/Mouse_MI_ATAC/Myh11/ArchR/data/cellColData_after_rm_doublets.Rds"

ccd <- readRDS(file=file)
rnms <- ccd@rownames

sham_myh11_name <- rnms[grepl("sham_myh11", ccd@rownames)]
day3_myh11_name <- rnms[grepl("day3_myh11", ccd@rownames)]
day10_myh11_name <- rnms[grepl("day10_myh11", ccd@rownames)]


sham_myh11 <- sapply(str_split(sham_myh11_name, "#"), function(x) return(x[2]))
day3_myh11 <- sapply(str_split(day3_myh11_name, "#"), function(x) return(x[2]))
day10_myh11 <- sapply(str_split(day10_myh11_name, "#"), function(x) return(x[2]))


#inputFiles <- c(
#        "../../data/DS8/outs/fragments.tsv.gz",
#        "../../data/DS9/outs/fragments.tsv.gz",
#        "../../data/DS10/outs/fragments.tsv.gz"
#
#)
#
#sampleNames <- c("day3_myh11", "day10_myh11", "sham_myh11")




dir_dict <- c(
	"sham" = "DS10",
	"day3" = "DS8",
	"day10" = "DS9")

names_dict <- list(
	"sham"   =sham_myh11,
	"day3"   =day3_myh11,
	"day10"  =day10_myh11)


for(condition in c("sham", "day3", "day10")){
    ds_dir <-glue("/data/scATA/Mouse_MI_ATAC/Myh11/data/{dir_dict[condition]}/outs")
	file <- file.path(ds_dir, "singlecell.csv")
	df <- read.csv(file, row.names=1)
	df <- df[names_dict[[condition]], ]
	df$is__cell_barcode = 1
	dir.create(glue("../data/{condition}"))
	fout_name <- file.path("../data", condition, glue("{condition}_myh11_singlecell.csv"))
	write.csv(df, file=fout_name, quote=F)
}

