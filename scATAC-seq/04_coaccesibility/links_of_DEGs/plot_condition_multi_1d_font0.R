library(Gviz)
library(glue)
library(rlist)
library(cicero)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(genomation)
library(optparse)


Which_Day = "Day3"


#celltype <- "Cardiomyocytes"
#celltype <- "Macrophages"
#tf <- "RUNX1" 

#celltype <- "Cardiomyocytes"

#celltype <- "Pericytes/vSMC"

#celltype <- "Fibroblasts"

##!!!ATTENTION
##!!!
## Smad3 if use grep, you can get MA0795.1.SMAD3	Smad2::Smad3, MA0795.1.SMAD3, SMAD2::SMAD3::SMAD4, etc.....
## FOSL1 ----> FOSL1::JUN, FOSL1::JUND ...... 
##!!!
##!!!

#Fibroblasts, collagen1 <- BACH1, BACH2, Smad3, Fosl1 and NFKB1
#"IRF9", "MAF", "MAFB", "CEBPB", "CEPBD"


gene_to_group_factor <- list("S100a8"= c("IRF9", "MAF", "Mafb", "CEBPB", "CEBPD"),
                             "S100a9" = c("IRF9", "MAF", "Mafb", "CEBPB", "CEBPD")) 




AllOptions <- function(){
    parser <- OptionParser()
    parser <- add_option(parser, c("-c", "--celltype"), type="character", default="Macrophages",
                    help="sample name for ATAC integration [default %default]",
                    metavar="character")
    parser <- add_option(parser, c("-g", "--group"), type="character", default="TRUE",
                    help="if use gene_to_group_factor [default %default]",
                    metavar="character")


    return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)
celltype  = pa$celltype
isgroup <- pa$group




condition_colors <- c("Sham"="#66c2a5", 
                      "Healthy"="#66c2a5", 
                      "Day3"="#fc8d62",
                      "Day10"="#800000")

### celltype related tfs
celltype_tfs <- c(

    "Cardiomyocytes" = c(),
    "Macrophages" = c(),
    "Pericytes/vSMC" = c(),
    "Fibroblasts" = c("BACH1, BACH2, Smad3, Fosl1 and NFKB1" )
)



# load in your data using rtracklayer
#gene_model <- readGFF("../../../Reference/refdata-cellranger-atac-GRCh38-1.1.0/genes/genes.gtf")
gene_model <- readGFF("../../../ref_cellranger/mm10-1.2/genes.gtf")
gene_model$chromosome <- gene_model$seqid
gene_model$gene <- gene_model$gene_id
gene_model$transcript <- gene_model$transcript_id
gene_model$symbol <- gene_model$gene_name
gene_model$exon <- gene_model$exon_id
gene_model$width <- gene_model$end - gene_model$start + 1
gene_model$feature <- gene_model$transcript_type
gene_model <- subset(gene_model, !is.na(transcript) & !is.na(exon))

TimePoints <- c("Healthy", "Day3", "Day10")

dict <- list(
    "Fibroblasts" = c("Col1a2"), 
    #"Fibroblasts" = c("Col1a1","Col1a2",  "Pdgfra", "Rgs5","Runx1","Bach2", "Smad3", "Fosl1", "Nfkb1"), 
    #  "Fibroblasts"  = c("Dcn", "Col5a1", "Medag", "Pdgfra", "Col6a1", "Col15a1", "Islr", "Loxl1", "Oaf", "Col5a3", "Acvrl1", "Col6a3", "Ptgir", "Sox9", "Gm17268", "Mrgprf", "S1pr3", "Col1a1", "Col6a5", "Eln", "Errfi1", "Gng8", "Dpep1", "Egfr", "Slc1a5", "Smoc2", "Adamts14", "Axl", "Ccl2", "Ddr2", "Has2", "Lum", "Adamts12", "Col6a2", "Fads2", "Fbln2", "Gfpt2", "Gli2", "Optc", "Spon2", "C1qtnf7", "Cygb", "Msc", "Pde9a", "Prpsap1", "Serping1", "Tnni1", "Tsku", "Cilp", "Dennd2a", "Enpp2", "Gas1", "Hmcn2", "Kcnt2", "Kirrel", "Ngf", "Pdgfrb", "Rtn4rl2", "Tcf21", "Antxr1", "Clec3b", "Col1a2", "Epha1", "Gm10974", "Itgb8", "Mgp", "Mmp2", "Osr1", "Rem1", "Rhoj", "Cdrt4", "Col3a1", "Fmo3", "Plekhh2", "Ror2", "Slc43a3", "Tex26", "Tril", "Ube2l6", "Ccdc170", "Glis2", "Il16", "Mrc2", "Olfml3", "Ppt2", "Sertad4", "Tmem240", "Vcan", "Adamts2", "Adamts5", "Akap12", "Col8a1", "Dclk1", "Dpt", "Gdf6", "Heg1", "Igfbp7", "Mfap4", "Omd", "Phlda3", "Postn", "Ppp1r14a", "Ramp2", "Rcn3", "Rgs17", "Syt17", "Vgll3", "Aldh1a1", "Aldh1a2", "Apcdd1", "Col16a1", "Dact2", "Enpp3", "Fap", "Fibin", "Gpm6b", "Matn2", "Mmp14", "Sod3", "Ssc5d", "Tek", "Tmem45a", "Abi3bp", "Adamts1", "Anpep", "Aspn", "Bicc1", "Cd34", "Clip4", "Cxcl1", "Dio3", "Emilin1", "Entpd2", "Gpc6", "Hrct1", "Il1r1", "Layn", "Lingo1", "Lmo1", "Negr1", "Nfatc4", "Pagr1a", "Pcdh19", "Pcdh9", "Ptgis", "Ptk7", "Serpina3n", "Thy1", "Tshr", "Twist2"),

    "Pericytes/vSMC" = c("Rgs5", "Rgs4", "Gprc5c", "Trpc3", "Pdgfrb", "Agap2", "Gucy1a3", "Pde8b", "Dtx3", "Gja4", "Gucy1b3", "Notch3", "Cox4i2", "P2ry14", "Sdr42e1", "Arhgef25", "Adra2b", "Fxyd2", "Rem1", "Higd1b", "S1pr3", "Tpm2", "Arhgef17", "Col25a1", "Dlc1", "Egflam", "Errfi1", "Tbx2", "Mrvi1", "Nr2f2", "Oaf", "Sspo", "Adcy6", "Apobec4", "Colec11", "Dgkb", "Dock6", "Ebf1", "Esam", "Foxs1", "Gpr87", "Iqsec1", "Kcna5", "Mark1", "Myh11", "Ngf", "Nostrin", "Nr2f1", "Pabpn1l", "Pou2f3", "Prickle2", "Rarres2", "Zfp467"), 


    "Macrophages" = c("S100a8", "S100a9"), 

    "Cardiomyocytes" = c("Ttn", "Fhl2", "Ryr2", "Myh6", "Myh7", "Gm17428", "Rbm20", "Actn2", "Cmtm5", "Ldb3", "Myl2", "Myocd", "Obscn", "Speg", "Prdm16", "Slc25a34", "Myl3", "Tnnt2", "Cacna1c", "Gm9955", "Des", "Fhod3", "Adprhl1", "Pln", "Scn5a", "Bmpr1a", "Fabp3", "Kcnk3", "Tuba8", "Hspb7", "Atp2a2", "Gm10570", "Kcnj11", "Mybpc3", "Nppb", "Rbm24", "Abcc8", "Klhl38", "Lmod2", "Plin5", "Myoz2", "Plin4", "Trim63", "Kcnh2", "Lrrc10", "Myh14", "Mylk3", "Ppara", "Actc1", "Ckm", "Gzmm", "Sema4b", "Tpm1", "Csrp3", "Irx4", "Jph2", "Nkx2-5", "Pip5k1b", "Pnmt", "Smyd2", "Als2cl", "Kcnc3", "Kcnd3", "Otof", "Ppargc1a", "Ccdc141", "Kcnip2", "Ntsr2", "Rapsn", "Rbm38", "Rcor2", "Tcap", "Adrb1", "Chrna2", "Corin", "Diras1", "Kcnj12", "Rilp", "Sh3rf2", "Xirp1", "Bik", "Car14", "Fgf1", "Gsg1l", "Lrrc52", "Myh7b", "Myom1", "Nmrk2", "Ppargc1b", "Sgcg", "Slc25a4", "Tnnc1", "Trpm4", "Ube2ql1", "Dsp", "Hrc", "Lingo3", "Lrrc2", "Orm2", "Orm3", "Slc26a6", "Slc38a3", "Tango2", "Tmem120b", "Tmem89", "Ankrd1", "Arhgef7", "Bcat2", "Cnksr1", "Dmpk", "Doc2g", "Fsd2", "Gpt", "Homer2", "Klhl31", "Klhl40", "Myom2", "Nhlrc4", "Rcan2", "Art5", "Arvcf", "Carns1", "Cdh2", "Cdk5rap1", "Cox6a2", "Eef1a2", "Fitm1", "Fitm2", "Hes6", "Inha", "Lpl", "Lrtm1", "Miip", "Mknk2", "Palld", "Pard6b", "Pigq", "Prox1", "Rpa1", "Rrad", "Rtn4r", "Tnni3", "Wnk2", "Dnaja2", "Dok7", "Dysf", "Gipc2", "Gipr", "Hhatl", "Hspb2", "Klhl30", "Lsp1", "Myo18b", "Popdc2", "Ppfia4", "Scx", "Tnnt1", "Trim72", "Asb18", "Bves", "Cryab", "Ctnnal1", "Frem2", "Ggt1", "Hcn4", "Ldhd", "Limch1", "Nrap", "Pgam2", "Pkdrej", "Prob1", "Rbbp8nl", "Samd4", "Slc27a1", "Slc44a2", "Acacb", "Alpk2", "Art1", "Cited4", "Coro6", "Cpt1b", "Esrrb", "Fam174b", "Foxo6", "Ip6k3", "Ky", "Macrod1", "Mypn", "Neo1", "Nsmce4a", "Pcdh7", "Pkp2", "Sh3bgr", "Slc25a13", "Syde2", "Asb10", "Btbd16", "C2cd4b", "Crip2", "Crybb1", "Dsg2", "Hand2", "Hspb8", "Igfals", "Il22ra1", "Irx3", "Mlxipl", "Nfib", "Obsl1", "Rhbdl1", "Rhbdl2", "Rims4", "Slc22a23", "Smim7", "Tkt", "Trabd", "Ankrd23", "Apobec2", "Asb2", "Baiap2l1", "Casq2", "Cdh16", "Cmya5", "Cryba4", "Cxadr", "Dhrs7c", "Ephx2", "Fbxo40", "Fign", "Flnc", "Guk1", "Hfe2", "Kcnd2", "Lmod3", "Mybphl", "Nrtn", "Nsg2", "Nxnl1", "Pdgfa", "Pgm5", "Pik3ip1", "Ppip5k2", "Prlh", "Rnf207", "Rpl3l", "Sbk3", "Sgca", "Slc1a7", "Slc38a7", "Sox6", "Synpo2l", "Tacc2", "Tbx5", "Tcea3", "Tfcp2l1", "Thrb", "Tmc7", "Tnfrsf18", "Tnni3k", "Trim54", "Ucn", "Wfdc1", "Abra", "Agt", "Alpk3", "Atp1b1", "Chrm2", "Ckmt2", "Ddo", "Fam189a2", "Fndc5", "Gpa33", "Itgb3bp", "Itgb6", "Kcnn2", "Lama5", "Ldhb", "Lrrc3b", "Lypd2", "Mlf1", "Mrps24", "Murc", "Nr0b2", "Oxct1", "Pank1", "Parm1", "Polq", "Proser3", "Sec14l5", "Slc8a2", "Smco1", "Smim6", "Sorbs1", "Sorcs2", "Sv2a", "Tmem108", "Wipf3", "Adamtsl5", "Bmp7", "Ces1d", "Crhr2", "Dtna", "Elp5", "Fam160a1", "Fam84a", "Gata6", "Hspb3", "Irx1", "Klhl34", "Lcn5", "Lrrc24", "Lsmem1", "Med12l", "Mettl24", "Mreg", "Mrpl39", "Nkain1", "Nt5c1a", "Paqr9", "Pdk4", "Perm1", "Ppp1r3a", "Prune2", "Scn4b", "Slc2a4", "Sptb", "Srpk3", "Stc2", "Tmem182", "Trdn", "Trim55", "Ucp3", "Wnk4", "Yod1", "Acsm5", "Actr3b", "Asb12", "Cav3", "Ccdc162", "Cdh13", "Clu", "Crnkl1", "Dsc2", "Erbb4", "Fam19a5", "Fbxl22", "Fbxo16", "Fgf13", "Filip1", "Filip1l", "Fn3k", "Fxyd1", "Fxyd7", "Gcom1", "Gja3", "Gja5", "Gmpr", "Gpd1", "Gpr75", "Hopx", "Hrasls", "Itga7", "Itgad", "Kcnj3", "Kcnn1", "Lrrc14b", "Lrrc39", "Mical3", "Myadml2", "Myo5c", "Myom3", "Myrip", "Myzap", "Npm2", "Opn5", "Pcp4l1", "Pde4d", "Perp", "Pla2g2c", "Pla2g4e", "Ppp2r4", "Prkcg", "Proser2", "Rarb", "Rgs6", "Rnf39", "RP24-559D1.1", "Shisa6", "Slc25a47", "Slco5a1", "Smpx", "Snrpn", "St5", "Stox2", "Stxbp5l", "Tfap2a", "Tmem171", "Tmem179", "Upb1", "Vldlr", "Wdr6", "Wif1", "Xirp2", "Yipf7", "Zbtb32", "Zfp445", "Zfp791")#,

#    "Endothelial cells" = c('Kdr", "Egfl7", "Flt1", "Podxl", "Prkch", "Ptprb", "Sh3tc2", "Aplnr", "Robo4", "Sox18", "Cav1", "Fabp4", "Gimap6", "Rapgef3", "Gimap5", "Pcdh12", "C1qtnf9", "Gimap4", "Ly6c1", "Sox17", "Dll4", "Eng", "Sema6a", "Esam", "Hspa12b", "Tcf15", "Aqp1", "Fgd5", "Timp4", "Gpihbp1", "Meox2", "Rassf9", "Shank3", "Sox7", "Tie1", "Wscd1", "Btnl9", "Myct1", "Nova2", "Rph3al", "Sdpr", "Tmem40", "Vsig2", "Cd300lg", "Cdh5", "Hyal2", "Syne1", "Ang2", "Cldn5", "Gimap8", "Ifitm5", "Slc43a3", "Sparcl1", "Tmem211", "Bcl6b", "Ccm2l", "Eml1", "Hrct1", "Pdgfb", "Rasip1", "Rtp3", "Sh3tc1", "She", "Slc16a13", "Tmem204", "Arhgef3", "Babam1", "Emcn", "Lmo7", "Lrrc32", "Nos3", "Nxpe2", "Pxdc1", "Rab3a", "Rasgrp3", "Slfn5", "Spry4", "Ssu2", "Tmem255b", "Vwa1", "Arhgef15", "Ctla2a", "Ctla2b", "Cttnbp2nl", "Dgke", "Ebf1", "Evl", "Fam124b", "Fam171a1", "Glb1l3", "Gm11037", "Grap", "Pcdh1", "Rtn4rl2", "Slc28a2", "Slc7a7", "Slc9a3r2", "Tead4", "Tox2", "Tshz2", "Amotl1", "Cldn15", "Dennd5b", "Ecscr", "Extl3", "Fchsd2", "Gng11", "Gpr4", "Grrp1", "Heg1", "Hint3", "Il20rb", "Med10", "Meox1", "Mgll", "Mmp11", "Nhsl1", "Nrarp", "Olfr282", "Prex2", "Psd3", "Rangap1", "Rhebl1", "Smco4", "Snn", "Spcs2", "Stk3", "Thsd1", "Tinagl1", "Tjp1", "Usp30", "Xiap", "Zdhhc1", "Zfp458")
)


list_macro <- c(
    "Macrophages", 
    "inflammatory macrophages", 
    "Anti-inflammatory macrophages"
)


get_gr_tf_sub <- function(minbp, maxbp, df_tf){
    gr <- GRanges(seqnames = chr, 
                  ranges = IRanges(start = minbp, end = maxbp))
    gr_tf <- GRanges(seqnames = df_tf$chromosome,
                     ranges = IRanges(start = df_tf$start,
                                      end = df_tf$end))
    gr_overlap <- findOverlaps(query = gr_tf, 
                              subject = gr, 
                              type = "within")
    gr_tf_sub <- gr_tf[gr_overlap@from, ]

    gr_overlap <- findOverlaps(query = gr_tf_sub, 
                               subject = gr_peaks, 
                               type = "any",
                               maxgap = 1000)

    gr_tf_sub <- gr_tf_sub[gr_overlap@from, ]

    gr_tf_sub <- resize(gr_tf_sub, width = 3000, fix = "center")

    gr_tf_sub
}

get_a_tf_df <- function(df1__, df2__, df3__, a_tf){
      df1 <- df1__[which(grepl(glue("\\.{a_tf}$"), df1__$V4)), ]
      df2 <- df2__[which(grepl(glue("\\.{a_tf}$"), df2__$V4)), ]
      df3 <- df3__[which(grepl(glue("\\.{a_tf}$"), df3__$V4)), ]

      df_tf <- rbind(df1, df2, df3)
      df_tf <- unique(df_tf)
      colnames(df_tf) <- c("chromosome", "start", "end", "name", "score", "strand")
      df_tf$symbol <- glue("{a_tf}")
      
      write.table(df_tf, file = glue("{a_tf}.bed"), col.names = FALSE, row.names = FALSE,
                  sep = "\t", quote = FALSE)

    df_tf
}

wrap_colors <- function(gr_tf_sub_list){
    colors <- c("#EFD8D7", "purple", "cyan", "brown", "skyblue")
    i = 1 
    color_list <- list()
    for(L in sapply(gr_tf_sub_list, length)){
         color_list[[i]] <- rep(colors[i], L)
         i = i + 1 
    }
    color_list[[i]] <- "red"

    return(unlist(color_list))
}


up_genes <- dict[[celltype]] 


df_links_list <- list()
for(TimePoint in TimePoints){
    df_links <- read.csv(glue("../save/{TimePoint}/detail_corr_pval_final_0.6_0.05_with_annotation.tsv"),
                            header = TRUE, sep = "\t")
    df_links$gene <- paste(df_links$gene_chr, 
                            df_links$gene_start, 
                            df_links$gene_start + 1,
                            sep = "_")
    df_links$pval_adjust_log10 <- -log10(df_links$padj)
    
    df_links_sub <- df_links[ c("peakName", "gene", "pval_adjust_log10", 
                                "gene_name", glue("{TimePoint}_clusters"))]

    colnames(df_links_sub) <- c("Peak1", "Peak2", "coaccess", "gene", "cluster")
    df_links_list[[TimePoint]] <- df_links_sub
}


df_links_Healthy <- subset(df_links_list[["Healthy"]], cluster == celltype)
df_links_Day3 <- subset(df_links_list[["Day3"]], cluster == celltype)
df_links_Day10 <- subset(df_links_list[["Day10"]], cluster == celltype)

if(celltype == "Macrophages"){
    df_links_Healthy <- subset(df_links_list[["Healthy"]], cluster %in% list_macro)
    df_links_Day3 <- subset(df_links_list[["Day3"]], cluster %in% list_macro)
    df_links_Day10 <- subset(df_links_list[["Day10"]], cluster %in% list_macro)

}
df_links_day_all_ <- rbind( df_links_Healthy, df_links_Day3, df_links_Day10)


scelltype <- str_replace_all(celltype, "/", ".")
data_track_Healthy <- DataTrack(range = as.character(glue("../../../BigWig/Healthy_{scelltype}.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Healthy", col = condition_colors['Healthy'],
                            ylim = c(0, 4))

data_track_Day3 <- DataTrack(range = as.character(glue("../../../BigWig/Day3_{scelltype}.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Day3", col = condition_colors['Day3'],
                            ylim = c(0, 4))

data_track_Day10 <- DataTrack(range = as.character(glue("../../../BigWig/Day10_{scelltype}.bw")), 
                            genome = "mm10", type = "h", 
                            name = "Day10", col = condition_colors['Day10'],
                            ylim = c(0, 4))


message("loading motifs healthy...", date() )
df1__ <- read.table(glue("../../../MotifMatch/Healthy_{scelltype}_mpbs.bed"))
message("loading motifs day3...", date() )
df2__ <- read.table(glue("../../../MotifMatch/Day3_{scelltype}_mpbs.bed"))
message("loading motifs day10...", date() )
df3__ <- read.table(glue("../../../MotifMatch/Day10_{scelltype}_mpbs.bed"))

for(up_gene in up_genes){
    tf <- gene_to_group_factor[[up_gene]]
    if(is.na(tf) && celltype == "Fibroblasts"){
        tf <- "RUNX2"
    }
    
    if(is.na(tf) && celltype == "Pericytes/vSMC"){
        tf <- "RUNX2"
    }
 
    if(is.na(tf) && celltype == "Macrophages"){
        tf <- "RUNX1"
    }
    if(is.na(tf) && celltype == "Cardiomyocytes"){
        tf <- "RUNX1"
    }

     if(is.na(tf) && celltype == "Endothelial cells"){
        tf <- "RUNX1"
    }
    tf_info_list <- list() 
    for(idx in 1:length(tf)){
        a_tf <- tf[idx] 
        colors <- c("blue", "purple", "cyan", "brown", "skyblue")         
        df_tf <- get_a_tf_df(df1__, df2__, df3__, a_tf)                  
        tf_track_1 <- GeneRegionTrack(range = df_tf,
                                     genome = "mm10",
                                     name = glue("{a_tf}"),
                                     shape = "box",
                                     col = colors[idx],
                                     fill = colors[idx],
                                     fontsize.group=0,
                                     fontcolor = "black")

         tf_info_list[[a_tf]] <- list(df_tf, tf_track_1)
    }
    
    # If I want specific color scheme, I can make a column of color names
    df_links_day_all_$conn_color <- "#800000"
    #df_links_Healthy$conn_color[df_links_Healthy$cluster == "Cardiomyocytes 2"] <- "#9A6324"
    


    message(sprintf("plotting links for gene %s", up_gene))
    df_links_day_all <- subset(df_links_day_all_, gene == up_gene)
    if(nrow(df_links_day_all) == 0){
       message(glue("no link: {up_gene}"))
       next 
    }
    df <- as.data.frame(str_split_fixed(df_links_day_all$Peak1, "_", 3))
    df$V1 <- as.character(df$V1)
    df$V2 <- as.numeric(as.character(df$V2))
    df$V3 <- as.numeric(as.character(df$V3))
    chr <- unique(df$V1)
    minbp <- min(df$V2) - 10000
    maxbp <- max(df$V3) + 10000
    viewpoint <- unique(df_links_day_all$Peak2)
    
    gr_peaks <- GRanges(seqnames = df$V1, IRanges(start = df$V2, end = df$V3))
    
    
    coaccess_cutoff = -log10(0.05)
    comparison_coaccess_cutoff = -log10(0.05)    
    ymax <- max(df_links_day_all$coaccess) + 1

    jumpout <- FALSE 
    link_tracks  <- NULL
    df_links_Which_day <- switch(Which_Day, "Healthy"=df_links_Healthy,
                                            "Day3" = df_links_Day3,
                                            "Day10" = df_links_Day10)

    tryCatch({

        link_tracks <- plot_connections(connection_df = df_links_day_all, 
                                    gene_model = gene_model,
                                    gene_model_shape = "smallArrow",
                                    collapseTranscripts = "longest",
                                    chr = chr, 
                                    minbp = minbp, 
                                    maxbp = maxbp,
                                    comparison_track = df_links_Which_day,
                                    alpha_by_coaccess = FALSE,
                                    # viewpoint = viewpoint,
                                    connection_ymax = ymax,
                                    comparison_ymax = ymax,
                                    coaccess_cutoff = coaccess_cutoff,
                                    comparison_coaccess_cutoff = comparison_coaccess_cutoff,
                                    connection_width = 1,
                                    comparison_connection_width = 1,
                                    connection_color = "conn_color",
                                    connection_color_legend = F,
                                    comparison_connection_color = "#9A6324",
                                    include_axis_track = FALSE,
                                    return_as_list = TRUE)

    }, error = function(cond) {
        message("Error:", "All ",celltype, "-", up_gene, " ", cond)
        jumpout <<- TRUE
    })   

    if(jumpout){
        next
    }
    
    gene_model_track <- link_tracks[[5]]

    jumpout <- FALSE
    df_links_Which_day$conn_color <- condition_colors[Which_Day]
    tryCatch({
        link_tracks <- plot_connections(connection_df = df_links_Which_day, 
                                    gene_model = gene_model,
                                    gene_model_shape = "smallArrow",
                                    collapseTranscripts = "longest",
                     chr = chr, 
                     minbp = minbp, 
                     maxbp = maxbp,
                     comparison_track = df_links_Which_day,
                     alpha_by_coaccess = FALSE,
                     # viewpoint = viewpoint,
                     connection_ymax = ymax,
                     comparison_ymax = ymax,
                     coaccess_cutoff = coaccess_cutoff,
                     comparison_coaccess_cutoff = comparison_coaccess_cutoff,
                     connection_width = 1,
                     comparison_connection_width = 1,
                     connection_color = "conn_color",
                     connection_color_legend = F,
                     comparison_connection_color = condition_colors[Which_Day],
                     include_axis_track = FALSE,
                     return_as_list = TRUE)
    }, error = function(cond) {
        message("Error:", "Which_Day ", celltype, "-", up_gene)
        jumpout <<- TRUE
    })   
    if(jumpout){
        next
    }
    link_track_Which_Day <- link_tracks[[1]]
    peak_track_Which_Day <- link_tracks[[2]]


    gene_axis_track <- GenomeAxisTrack(fontsize = 4)
    
    link_track_Which_Day@dp@pars$fontsize <- 12
    data_track_Healthy@dp@pars$fontsize <- 12
    data_track_Day3@dp@pars$fontsize <- 12
    data_track_Day10@dp@pars$fontsize <- 12
    
    trackList <- list(link_track_Which_Day,
                      peak_track_Which_Day,
                      data_track_Healthy,
                      data_track_Day3,
                      data_track_Day10,
                      gene_model_track, 
                      lapply(tf_info_list, function(x)x[[2]]),
                      gene_axis_track)
    

    trackList <- list.flatten(trackList)

    view_start <- stringr::str_split_fixed(viewpoint, "_", 3)[, 2]
    view_end <- stringr::str_split_fixed(viewpoint, "_", 3)[, 3]
    

    gr_tf_sub_list <- lapply(tf_info_list, function(x) get_gr_tf_sub(minbp, maxbp, x[[1]]) )   

    gr_promoter <- GRanges(seqnames = chr,
                     ranges = IRanges(start = as.numeric(view_start),
                                      end = as.numeric(view_end)))
    
    gr_view <- Reduce(c,  list.flatten(list(gr_tf_sub_list, gr_promoter)))

    
    trackList <- HighlightTrack(trackList = trackList, 
                                range = gr_view,
                                chromosome = chr,
                                col = wrap_colors(gr_tf_sub_list), 
                                fill = c(rep("#EFD8D7", length(wrap_colors(gr_tf_sub_list)) -1), "red"), 
                                inBackground = FALSE, 
                                alpha = 0.3)
    dir.create(glue("{scelltype}"), recursive=T) 
    tfs <- paste0(tf, collapse="_")
    pdf(glue("{scelltype}/condition_{scelltype}_{up_gene}_{tfs}_noTFnames.pdf"), height = 8, width = 8)
    plotTracks(trackList, title.width = 0.5, showTitle = TRUE, 
               from = minbp, to = maxbp, chromosome = chr, 
               sizes = c(0.4, 0.1, 0.3, 0.3, 0.3, 0.3, rep(0.2, length(tf_info_list)), 0.05), 
               transcriptAnnotation = "symbol", background.title = "white", 
               col.border.title = "transparent", lwd.border.title = "transparent", 
               col.axis = "black", fontcolor.legend = "black",
               innerMargin = 3)
    dev.off()

}
