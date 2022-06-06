library(data.table)
library(doParallel)
library(foreach)

registerDoParallel(3)


dict <- list(
	"sham" = c("DS1", "DS7"),
	"day3" = c("DS3", "DS4"),
	"day10" = c("DS5", "DS6")
)

timepoint <- c("sham", "day3", "day10")

cmtx_dir <- "/hpcwork/izkf/projects/SingleCellOpenChromatin/exp/Mingbo/Mouse_scATAC/CountMatrix"
cell_dir <- "../data/"


foreach (tm = timepoint) %dopar% {
	message("processing", tm, date())
	out_dir <- sprintf("../data/%s/Fragments", tm)
	samples <- dict[[tm]]
	female_dir <- file.path(cmtx_dir, samples[1], "/outs/fragments.tsv.gz")
	male_dir <- file.path(cmtx_dir, samples[2], "/outs/fragments.tsv.gz")

	# fragments
	female_dt <- fread(female_dir, col.names = c('chr', 'start', 'end', 'cell', 'count'))
	male_dt <- fread(male_dir, col.names = c('chr', 'start', 'end', 'cell', 'count'))
	# filtered cells
	female_df <- read.csv(file=file.path(cell_dir, tm, sprintf("%s_female_singlecell.csv", tm)),
							row.names=1)
	female_df$sex="female"
	male_df <- read.csv(file=file.path(cell_dir, tm, sprintf("%s_male_singlecell.csv", tm)),
							row.names=1)

	male_df$sex="male"
	# change male barcodes
	male_dt$cell <-  gsub("-1", "-2", male_dt$cell)
	rownames(male_df) <- gsub("-1", "-2", rownames(male_df))

	# filter fragments
	female_dt <- female_dt[cell %in% rownames(female_df)] 
	male_dt <- male_dt[cell %in% rownames(male_df)]
	
	#female_dt  <- female_dt[chr != "chrY"]	
	#male_dt <- male_dt[chr != "chrY"]


	# bind female & male		
	dt <- rbind(female_dt, male_dt)
	fwrite(dt, file=file.path(out_dir, "fragments.tsv"), sep="\t", col.names=F)


	# genereate single cell file
	df <- rbind(female_df, male_df)
	write.csv(df, file=file.path(cell_dir, tm, sprintf("singlecell.csv", tm)), quote=F)
}


