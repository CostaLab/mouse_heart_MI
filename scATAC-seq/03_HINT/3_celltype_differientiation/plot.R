library(ggplot2)
library(cowplot)
library(gridExtra)
library(reshape2)
library(RColorBrewer)


colors = c('Cardiomyocytes' = '#800000',
           'Fibroblasts' = '#911eb4',
           'Myofibroblasts' = '#f032e6',
           'Endothelial.cells' = '#000075',
           'Macrophages' = '#e6194B',
           'Anti.inflammatory.macrophages' = '#568198',
           'inflammatory.macrophages' = '#469990',
           'Pericytes' = '#f58231',
           'Adipocytes' = '#000000',
           'Neuronal cells' = '#42d4f4',
           'Erythrocytes' = '#999999',
           'Proliferating cells' = '#999999',
           'Damaged endothelial cells' = '#999999',
           'vSMC' = '#aaffc3',
           'Pericytes/vSMC' = '#aaffc3')


tm="Day3"



smtm <- function (vector, n=5, r=5){
    k=length(vector)
    z=((vector-mean(vector))/var(vector))
    vector[rank(z)>(k-r)]=NA
    res=c()
    for(i in 1:length(vector)){
        a=max(1,i-round(n/2))
        b=min(length(vector),ceiling(i+n/2))
        res=c(res,
              mean(vector[a:b],na.rm=TRUE))
    }
    return(res)
}


tf_0.2 <- read.csv(sprintf("../heatmap/%s_TF_filter_0.2.csv", tm), header = TRUE, row.names=1)

tf_list <- rownames(tf_0.2)

for (tf in tf_list) {
    filename <- paste0("../Diff/Lineplots/", tf, ".txt")
    df <- read.csv(filename, header = TRUE, sep="\t")
    df$Pro_cells <- NULL
    df$Pos <- seq(from = -99, to = 100)
    for (i in 1:(length(df)-1)){
        df[,i]=smtm(df[,i])
    } 
    df <- melt(df, id.var = c("Pos"))
    
    df$variable <- stringr::str_replace_all(df$variable, c("Pericytes.vSMC" = "Pericytes/vSMC"))
    
    p <- ggplot(data = df, aes(x = Pos, y = value, color = variable)) +
        geom_line(size=1) + xlab("") + ylab("") + 
        scale_color_manual(values = colors) +
        ggtitle(tf) +
        theme_cowplot() + 
        theme(legend.title = element_blank())
    
    pdf(paste0(tf, ".pdf"), width = 8, height = 4)    
    print(p)
    dev.off()
}
