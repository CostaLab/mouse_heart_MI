suppressPackageStartupMessages(library(cicero))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(glue))


topcelltype <- "Macrophages"
agg_cds <- readRDS(file=glue("{topcelltype}_agg_cds.Rds"))

orig_atac <- readRDS(file="../../Seurat3/IntegrationWithscOpen/obj.integrated_annotated_day3_myo.Rds")
meta <- orig_atac@meta.data 
meta <- meta[colnames(agg_cds), ] 

agg_cds@phenoData@data$gender <- meta$gender.ident
agg_cds@phenoData@data$celltype <- meta$celltype

tbl <- table(agg_cds@phenoData@data$celltype, agg_cds@phenoData@data$State)
pct <- round(sweep(100*tbl, 2, colSums(tbl), `/`))  
print(pct) ### State 6 account for the most 


pdf(glue("{topcelltype}_pseudotime.pdf"))
plot_cell_trajectory(agg_cds, color_by = "timepoint", cell_size=0.3)
plot_cell_trajectory(agg_cds, color_by = "gender", cell_size=0.3)
plot_cell_trajectory(agg_cds, color_by = "celltype", cell_size=0.3)
plot_cell_trajectory(agg_cds, color_by = "State", cell_size=0.3)
agg_cds <- orderCells(agg_cds, root_state = 6)
plot_cell_trajectory(agg_cds, color_by = "Pseudotime", cell_size=0.3)

dev.off()



topcelltype <- "Fibroblasts"
agg_cds <- readRDS(file=glue("{topcelltype}_agg_cds.Rds"))
meta <- orig_atac@meta.data 
meta <- meta[colnames(agg_cds), ] 

agg_cds@phenoData@data$gender <- meta$gender.ident
agg_cds@phenoData@data$celltype <- meta$celltype

tbl <- table(agg_cds@phenoData@data$celltype, agg_cds@phenoData@data$State)
pct <- round(sweep(100*tbl, 2, colSums(tbl), `/`))  
print(pct) ### State 2 account for the most 


pdf(glue("{topcelltype}_pseudotime.pdf"))
plot_cell_trajectory(agg_cds, color_by = "timepoint", cell_size=0.3)
plot_cell_trajectory(agg_cds, color_by = "gender", cell_size=0.3)
plot_cell_trajectory(agg_cds, color_by = "celltype", cell_size=0.3)
plot_cell_trajectory(agg_cds, color_by = "State", cell_size=0.3)
agg_cds <- orderCells(agg_cds, root_state = 2)
plot_cell_trajectory(agg_cds, color_by = "Pseudotime", cell_size=0.3)

dev.off()


