library(glue)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(reshape2)
library(RColorBrewer)


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
  message("celltype: ", celltype)

  fname <- glue("../barplot/{celltype}_Female.txt")

  if (file.size(fname) > 0){
    x <- read.csv(fname, header = F, stringsAsFactors = F)  
  }

  fname <- glue("../barplot/{celltype}_Male.txt") 
  if (file.size(fname) > 0){
    y <- read.csv(fname, header = F, stringsAsFactors = F)  
  }

  TFs <- c(x$V1, y$V1)

  TFs <- stringr::str_replace_all(TFs, c("\\("="_", "\\)"=""))


  for (tf in TFs) {
      filename <- glue("../Diff/{celltype}/Lineplots/{tf}.txt")
      df <- read.csv(filename, header = TRUE, sep="\t")
      df$Pro_cells <- NULL
      df$Pos <- seq(from = -99, to = 100)
      df <- melt(df, id.var = c("Pos"))
      
      df$variable <- stringr::str_replace_all(df$variable, c("Pericytes.vSMC" = "Pericytes/vSMC"))
      
      p <- ggplot(data = df, aes(x = Pos, y = value, color = variable)) +
          geom_line(size=1) + xlab("") + ylab("") + 
          scale_color_manual(values = colors) +
          ggtitle(tf) +
          theme_cowplot() + 
          theme(legend.title = element_blank())
     
      dir.create(celltype)
      pdf(file.path(celltype, paste0(tf, ".pdf")), width = 8, height = 4)    
      print(p)
      dev.off()
  }
}

