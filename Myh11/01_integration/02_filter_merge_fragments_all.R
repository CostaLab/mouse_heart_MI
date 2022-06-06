library(data.table)
library(doParallel)
library(foreach)
library(future.apply)

registerDoParallel(6)


dict <- list(
	"sham" = c("DS10"),
	"day3" = c("DS8"),
	"day10" = c("DS9")
)

bc_tail_dict <- list(
    "DS8" = "-7",
    "DS9" = "-8",
    "DS10" = "-9")

time_dict <- list(
    "DS8" = "day3",
    "DS9" = "day10",
    "DS10" = "sham")

sex_dict <- list(
    "DS8" = "myh11",
    "DS9" = "myh11",
    "DS10" = "myh11")


data_vec <- names(bc_tail_dict)


timepoint <- c("sham", "day3", "day10")

cmtx_dir <- "/data/scATA/Mouse_MI_ATAC/Myh11/data/"

cell_dir <- "../data/"


fragments_file_vec <- file.path(cmtx_dir, data_vec, "/outs/fragments.tsv.gz")
fragments_dt_list <- future_lapply(fragments_file_vec,
            function(x){
                fread(x, col.names = c('chr', 'start', 'end', 'cell', 'count'))})


cell_file_vec <- sapply(1:3,
            function(x) file.path(cell_dir, time_dict[[x]],
                    sprintf("%s_%s_singlecell.csv", time_dict[[x]], sex_dict[[x]]))  )

cell_df_list <- future_lapply(cell_file_vec,
            function(x)  read.csv(file=x, row.names=1) )

cell_df_list <- lapply(1:3,
            function(x) {
               df <- cell_df_list[[x]]
               rownames(df) <- gsub("-1",bc_tail_dict[[x]], rownames(df))
               return(df)
            })


filtered_fragments_dt_list <- lapply(1:3,
            function(x){
                dt <- fragments_dt_list[[x]]
                df <- cell_df_list[[x]]
                dt$cell <- gsub("-1",bc_tail_dict[[x]], dt$cell)
                dt <- dt[cell %in% rownames(df) ]
            })



dt <- do.call(rbind, filtered_fragments_dt_list)
fwrite(dt, file="../data/fragments.tsv", sep="\t", col.names=F)


library(Rsamtools)


outf <- bgzip(file = "../data/fragments.tsv")
index.file <- indexTabix(file = paste0(outf),
                 format = 'bed',
                 zeroBased = TRUE)




df <- do.call(rbind, cell_df_list)
write.csv(df, file="../data/singlecell.csv", quote=F)
