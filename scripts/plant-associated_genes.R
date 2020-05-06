rm(list=ls())

#load config 
source("/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/scripts/_config.R")
source("/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/scripts/_KO_enrichment_function.R")
source("/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/scripts/_ko_plot_function.R")
##


####plot
optns <- theme(
  axis.text.x = element_text(margin = margin(c(20, 0, 0, 0)), size = 60, face = "bold", family = "Helvetica", angle = 0, hjust = .5, vjust = .5 ),
  axis.text.y = element_text(margin = margin(c(0, 20, 0, 0)) ,size = 60, face = "bold", family = "Helvetica"), 
  axis.title = element_blank(),
  axis.ticks.y = element_line(size = 3),
  axis.ticks.x = element_line(size = 3),
  axis.ticks.length = unit(.5, "cm"),
  axis.line  = element_line(size = 3),
  panel.background = element_blank(),
  legend.text = element_text(size = 40, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 40, face = "bold", family = "Helvetica"),
  legend.key = element_rect(fill = NA , size= .1)
) 

####

strain <- c("Leaf1","Leaf130", "Leaf155","Leaf177", "Leaf187", "Leaf176", "Leaf404", "Root935", "Soil763", "Pto",  "D36E")
phyla <- c("Actinobacteria1", "Pseudomonas","Alphaproteobacteria","Burkholderiales","Bacillales"
           ,"Bacteroidetes","Bacteroidetes","Bacteroidetes","Actinobacteria1" , "Pseudomonas", "Pseudomonas")


#combined file directory 
dir_data <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/rnaseq_summary_each_strain/"
out_dir <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/PA_DEG/expression_with_PA/"
out_dir2 <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/PA_DEG/box_plot/"
out_dir3 <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/PA_DEG/ko_enrich/"




for (i in c(1:length(strain))){

  data <- read.delim(paste(dir_data, "A346_rnaseq_summary_", strain[i], ".txt", sep = "" ), header = T, row.names = 1)
  
  #====PA gene data (Levy et al 2017)
  
  genome_dir <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A310_psl_commensal_rnaseq_analysis/Levy_etal_data/"
  
  pa_threshold = 1 #significant in at least one statistical test
  
  dat <- read.delim(paste(genome_dir, "ko_", phyla[i], ".txt", sep = "" ), header = T, row.names = 1)
  rownames(dat) <- gsub("KO:", "", rownames(dat))
  dat <- dat[, -c(1:7)]
  
  
  pa_test <- apply(dat, 1, function(x){sum(x=="Y")}) #number of statistical tests that were significant
  pa <- names(pa_test)[pa_test > pa_threshold] %>% as.character()
  
  data[, ncol(data)+1] <- "nonPA"
  b <- grep(paste(pa, collapse = "|"), data[, 6]) %>% na.omit() %>% as.numeric()
  data[b,ncol(data)] <- "PA"
  
  colnames(data)[ncol(data)] <- "Plant Associated"
  
  data_for_save <- data[, c(5,6,7,1,2,3,4,8)]
  
  write.table(data_for_save, file=paste(out_dir, "A346_gene_expression_withPA_", strain[i] , ".txt", sep = ""), row.names=T, col.names=NA,sep="\t", quote=F)
  
  #=============#=============#=============#=============#=============#=============
  ###link PA and gene regulation
  data2 <- data[, c(6, 4, 8)]
  data2 <- data2[!data2[, 1]=="", ]
  
  out_pa <- table(data2[data2[,3]=="PA" ,2])
  out_nonpa <- table(data2[data2[,3]=="nonPA" ,2])
  
  out_pa2 <- out_pa/sum(out_pa)
  out_nonpa2 <- out_nonpa/sum(out_nonpa)
  names(out_pa2) <- c("PA down", "PA nonDEG","PA up")
  names(out_nonpa2) <- c("nonPA down", "nonPA nonDEG","nonPA up")
  
  print(out_pa2)
  print(out_nonpa2)
  
  #####drawing a box plot######
  x <- data[, c(3, 4, 6, 8)]
  x <- x[!x[, 3]== "", ]
  x <- x[order(x[, 1]), ] %>% na.omit()
  x <- x[order(x[, 4]), ] %>% na.omit()
  
  colnames(x) <- c("FC", "DEG", "KO","PA")
  
  q <- ggplot(x, aes(x = PA, y = FC)) + geom_violin( size = 1.5, fill = NA)
  q <- q + scale_fill_grey() + theme_classic()
  q <- q + optns
  
  #setting DEG colors
  c3 <- gsub("nonDEG", "black", x[, 2])
  c3 <- gsub("down", "deepskyblue3", c3)
  c3 <- gsub("up", "deeppink", c3)%>% na.omit() %>% as.character()
  
  q <- q + geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', position = position_jitter(.1), 
                        binwidth = .1 , alpha = 1, dotsize = 1,  fill = c3, color = NA)+ 
    xlab("") + ylab("")+ theme(legend.position = "none") 
  
  ggsave(paste(out_dir2, "Expression_boxplot_PA_DEG_", strain[i],".png"),  height = 10, width = 10, q)
  
  
  #####KO enrichment analysis for genes that are PA/nonPA and up/down-regulated#####
  infile <- data2 %>% na.omit()
  infile[, 2] <- paste(infile[, 2], infile[, 3], sep = "_")
  pt <- unique(infile[, 2]) %>% as.character()
  
  out <- data.frame("KO" = NA, "pval" = NA, "FoldEnrichment" = NA, "genes" = NA, "size" = NA, "cluster" = NA)
  
  #ko enrichment summary table
  for (id in c(1:length(pt))){
    
    indata <- infile[infile[, 2]==pt[id] ,1]
    indata <- as.data.frame(indata)
    
    hypr <- ko.enrich.strain(indata, strain[i])
    if (nrow(hypr)>0){
      hypr[, 6] <- pt[id]
      colnames(hypr)[6] <- "cluster"
      out <- rbind(out, hypr)
    }
  }
  out <- out[-1, ]
  
  write.table(out, file=paste(out_dir3, "A346_PA_DEG_",  strain[i] ,"_KOenrich.txt", sep = ""), row.names=F,  sep="\t", quote=F)
  
  
  ######hypergeometric test#####
  p <- phyper(q = out_pa[3], m = out_pa[3]+out_nonpa[3], n = sum(out_pa)+sum(out_nonpa)-out_pa[ 3]-out_nonpa[3], k = sum(out_pa), lower.tail = F, log.p = F )
  pp <- phyper(q = out_nonpa[1], m = out_pa[1]+out_nonpa[1], n = sum(out_pa)+sum(out_nonpa)-out_pa[ 1]-out_nonpa[1], k = sum(out_nonpa), lower.tail = F, log.p = F )
  pp <- c(p, pp)
  pp <- p.adjust(pp, method = "fdr") %>% as.data.frame()
  print(pp)
}
