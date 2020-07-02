library(stringr)

file="/work/sz753404/Mouse_scATAC/ArchR/data/cellColData_after_rm_doublets.Rds"

ccd <- readRDS(file=file)
rnms <- ccd@rownames

sham_female_name <- rnms[grepl("sham_female", ccd@rownames)]
sham_male_name <- rnms[grepl("sham_male", ccd@rownames)]

day3_female_name <- rnms[grepl("day3_female", ccd@rownames)]
day3_male_name <- rnms[grepl("day3_male", ccd@rownames)]


day10_female_name <- rnms[grepl("day10_female", ccd@rownames)]
day10_male_name <- rnms[grepl("day10_male", ccd@rownames)]


sham_female <- sapply(str_split(sham_female_name, "#"), function(x) return(x[2]))
sham_male <- sapply(str_split(sham_male_name, "#"), function(x) return(x[2]))

day3_female <- sapply(str_split(day3_female_name, "#"), function(x) return(x[2]))
day3_male <- sapply(str_split(day3_male_name, "#"), function(x) return(x[2]))

day10_female <- sapply(str_split(day10_female_name, "#"), function(x) return(x[2]))
day10_male <- sapply(str_split(day10_male_name, "#"), function(x) return(x[2]))



dir_dict <- c(
	"sham_female" = "DS1",
	"day3_female" = "DS3",
	"day3_male" = "DS4",
	"day10_female" = "DS5",
	"day10_male" = "DS6",
	"sham_male" = "DS7")

names_dict <- list(
	"sham_female"   =sham_female,
	"day3_female"   =day3_female,
	"day3_male"     =day3_male,
	"day10_female"  =day10_female,
	"day10_male"    =day10_male,
	"sham_male"     =sham_male)


ds_dir <- "/hpcwork/izkf/projects/SingleCellOpenChromatin/exp/Mingbo/Mouse_scATAC/CountMatrix/"
for(condition in c("sham", "day3", "day10")){
	gender <- "male"	
	key <- paste0(condition, "_male")
	file <- file.path(ds_dir, dir_dict[key], "outs/singlecell.csv")
	df <- read.csv(file, row.names=1)
	df <- df[names_dict[[key]], ]
	df$is__cell_barcode = 1
	fout_name <- file.path("../data", condition, paste0(key, "_singlecell.csv"))
	write.csv(df, file=fout_name, quote=F)


	gender <- "female"	
	key <- paste0(condition, "_female")
	file <- file.path(ds_dir, dir_dict[key], "outs/singlecell.csv")
	df <- read.csv(file, row.names=1)
	df <- df[names_dict[[key]], ]
	df$is__cell_barcode = 1
	fout_name <- file.path("../data", condition, paste0(key, "_singlecell.csv"))
	write.csv(df, file=fout_name, quote=F)
	
}

