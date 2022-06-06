mute <- suppressPackageStartupMessages
mute(library(Seurat))
mute(library(stringr))
mute(library(ArchR))


#obj <- readRDS("../../Seurat3_myh11/IntegrationWithscOpen/save_res/annotation_1.0.Rds")
obj <- readRDS("../../Seurat3_myh11/IntegrationWithscOpen_final_rm_res0.1_456/save_res/trans_add_tdtomato_annotation_res1.0_remove_res0.1_456.Rds")

Idents(obj) <- "annotation"
fib <- subset(obj, idents=c("Fibroblasts"))
myo <- subset(obj, idents=c("Myofibroblasts"))
a_sub <- subset(obj, idents=c("Fibroblasts", "Myofibroblasts"))



df <- a_sub@meta.data
df$Barcode <- rownames(df)

df$Barcode <-stringr::str_replace(df$Barcode, "7", "1")
df$Barcode <-stringr::str_replace(df$Barcode, "8", "1")
df$Barcode <-stringr::str_replace(df$Barcode, "9", "1")

df$Sample <- df$time.ident
df$celltype <- df$annotation


df <- df[, c("Sample", "Barcode", "celltype")]

df_Sham <- subset(df, Sample == "Sham")
df_Day3 <- subset(df, Sample == "Day3")
df_Day10 <- subset(df, Sample == "Day10")

df_Sham$Sample <- NULL
df_Day3$Sample <- NULL
df_Day10$Sample <- NULL


write.table(df_Sham, file = "Sham.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(df_Day3, file = "Day3.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(df_Day10, file = "Day10.txt", row.names = FALSE, sep = "\t", quote = FALSE)



