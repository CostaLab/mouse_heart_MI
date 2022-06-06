library(data.table)
library(doParallel)
library(foreach)

registerDoParallel(3)


dict <- list(
	"sham" = c("DS10"),
	"day3" = c("DS8" ),
	"day10" = c("DS9")
)

timepoint <- c("sham", "day3", "day10")

cmtx_dir <- "/data/scATA/Mouse_MI_ATAC/Myh11/data/"
cell_dir <- "../data/"


foreach (tm = timepoint) %dopar% {
	message("processing", tm, date())
	out_dir <- sprintf("../data/%s/Fragments", tm)
	if (!dir.exists(out_dir)){
	    dir.create(out_dir, recursive=T)
	}
	samples <- dict[[tm]]
	myh11_dir <- file.path(cmtx_dir, samples[1], "/outs/fragments.tsv.gz")

	# fragments
	myh11_dt <- fread(myh11_dir, col.names = c('chr', 'start', 'end', 'cell', 'count'))

	# filtered cells
	myh11_df <- read.csv(file=file.path(cell_dir, tm, sprintf("%s_myh11_singlecell.csv", tm)),
							row.names=1)
	myh11_df$sex="myh11"

	# filter fragments
	myh11_dt <- myh11_dt[cell %in% rownames(myh11_df)]

	#myh11_dt  <- myh11_dt[chr != "chrY"]

	fwrite(myh11_dt, file=file.path(out_dir, "fragments.tsv"), sep="\t", col.names=F)

	# genereate single cell file
	write.csv(myh11_df, file=file.path(cell_dir, tm, sprintf("singlecell.csv", tm)), quote=F)
}


