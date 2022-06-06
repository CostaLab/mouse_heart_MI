library(ggplot2)
library(glue)
library(cowplot)
library(gridExtra)
library(reshape2)
library(RColorBrewer)


colors = c('Cardiomyocytes' = '#800000',
           'Fibroblasts' = '#911eb4',
           'Endothelial cells' = '#000075',
           'Macrophages' = '#e6194B',
           'Pericytes.vSMC' = '#aaffc3')
celltypes=names(colors)


condition_colors <- c(
    Healthy  = "#66c2a5",
    Day3  = "#fc8d62",
#    Day10  = "#fc8d62"
    Day10 = "#8da0cb"
)




#smtm <- function (vector, n=5, r=5){
#    k=length(vector)
#    z=((vector-mean(vector))/var(vector))
#    vector[rank(z)>(k-r)]=NA
#    res=c()
#    for(i in 1:length(vector)){
#        a=max(1,i-round(n/2))
#        b=min(length(vector),ceiling(i+n/2))
#        res=c(res,
#              mean(vector[a:b],na.rm=TRUE))
#    }
#    return(res)
#}


for(celltype in celltypes){
    message("celltype:", celltype)
    dir.create(celltype)

    tf <- read.csv(sprintf("../heatmap/%s_TF_filter_0.02.csv", celltype), header = TRUE, row.names=1)
    #tf_pval <- read.csv(sprintf("../heatmap/pval_%s_TF_filter.csv", celltype), header = TRUE, row.names=1)
    #same_name <- intersect(rownames(tf_0.02), rownames(tf_pval))
    #if (length(same_name) > 0){
    #    idx <- which(rownames(tf_pval) %in% same_name)
    #    tf_pval <- tf_pval[-idx, ]
    #}
    #tf <- rbind(tf_0.02, tf_pval)

    tf_list <- rownames(tf)
    tf_list <- unique(tf_list)

    for (tf in tf_list) {
        filename <- glue("../Diff/{celltype}/Lineplots/{tf}.txt")
        df <- read.csv(filename, header = TRUE, sep="\t")
        df <- df[, c("Healthy", "Day3", "Day10")]
        df$Pro_cells <- NULL
        df$Pos <- seq(from = -99, to = 100)
        #for (i in 1:(length(df)-1)){
        #    df[,i]=smtm(df[,i])
        #}
        df <- melt(df, id.var = c("Pos"))

        df$variable <- stringr::str_replace_all(df$variable, c("Pericytes.vSMC" = "Pericytes/vSMC"))

        p <- ggplot(data = df, aes(x = Pos, y = value, color = variable)) +
            geom_line(size=1) + xlab("") + ylab("") +
            scale_color_manual(values = condition_colors) +
            ggtitle(tf) +
            theme_cowplot() +
            theme(legend.title = element_blank())

        pdf(file.path(celltype, paste0(tf, ".pdf")), width = 8, height = 4)
        print(p)
        dev.off()
}





}

