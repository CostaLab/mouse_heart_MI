library(ArchR)


 proj <- loadArchRProject(path = "./Fib", showLogo = FALSE)
 ## write table for HINT

 df <- as.data.frame(proj@cellColData)
 rownames(df) <- sapply(stringr::str_split(rownames(df), "#"), function(x)x[2])
 write.table(df, "meta.tsv",sep="\t", quote = F)

