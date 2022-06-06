library(ggplot2)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(optparse)
library(glue)
colors <- c(
    Female = "#9A6324", 
    Male =  "#800000"
)
celltypes <- c("Cardiomyocytes",
                 "Endothelial cells",
                 "Fibroblasts",
                 "Macrophages",
                 "Pericytes.vSMC")


for(celltype in celltypes){
    tbl_txt <- glue("../Diff/{celltype}/{celltype}_statistics.txt")
    df <- read.csv(tbl_txt, sep="\t", header = TRUE)
    df.plot <- subset(df, Num > 1000 & P_values < 0.05)
    df.plot <- subset(df, Num > 1000 & P_values < 0.05)
    df.plot$TF_Activity = -1 * df.plot$TF_Activity
    df.plot$Specific <- ifelse(df.plot$TF_Activity>0, "Female", "Male")
    df.plot$Specific <- factor(df.plot$Specific, levels=c("Female", "Male"))

    p <- ggplot(df.plot, aes(x = reorder(Motif, TF_Activity), y = TF_Activity, fill = Specific)) +
        geom_bar(stat = "identity")+
        scale_fill_manual(values = colors) +
        coord_flip() +
        ggtitle(ifelse(celltype=="Pericytes.vSMC", "Pericytes/vSMC", celltype)) + 
        xlab("Motifs") +
        ylab("TF activity change") +
        theme_cowplot() +
        theme(legend.title = element_blank())

    pdf(file = glue("{celltype}_diff_TF_Female.vs.Male.pdf"), width = 6, height = 6)
    print(p)
    dev.off()
    df.Female<- subset(df.plot, Specific == "Female")
    df.Male<- subset(df.plot, Specific == "Male")

    write.table(df.Female$Motif, file = glue("{celltype}_Female.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(df.Male$Motif, file = glue("{celltype}_Male.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

}

