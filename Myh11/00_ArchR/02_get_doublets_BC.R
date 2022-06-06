library(stringr)

before <- readRDS(file="../data/cellColData_before_rm_doublets.Rds")
after <- readRDS(file="../data/cellColData_after_rm_doublets.Rds")

rnm_before <- before@rownames
rnm_after <- after@rownames

doublets <- setdiff(rnm_before, rnm_after)
splits <- str_split(doublets, "#")

doublets_list = list()

for (x in 1:length(splits)){
    element <- splits[[x]]
    key <- element[1]
    val <- element[2]
    doublets_list[[key]] <- c(doublets_list[[key]], val)
}

saveRDS(doublets_list, file="../data/doublets_list.Rds")

