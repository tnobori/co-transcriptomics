#This script is for performing KO enrichment analysis on genes whose expression
#is significantly changed in planta compared with in vitro in various bacterial strains.
#The results are used for generating heatmaps in Fig. 2A and Fig. S3A.

#created by Tatsuya Nobori
#tnobori@salk.edu

rm(list=ls())

#load config 
source("/scripts/_config.R")
source("/scripts/_ko_enrichment_function.R")

#directories
dir_data <- "/data/normalized_bacterial_RNA-seq_data/"
out_dir <- "/output/dir/"

strain <- c("Leaf1","Leaf130", "Leaf155","Leaf177", "Leaf187", "Leaf176", "Leaf404", "Root935", "Soil763", "Pto", "D36E", "MM")

for (i in c(1:length(strain))){
  
  data <- read.delim(paste(dir_data, "A346_rnaseq_summary_", strain[i], ".txt", sep = "" ), header = T, row.names = 1)
  
  up <- data[data[, 4]=="up"  ,]
  down <-data[data[, 4]=="down" ,]
  
  indata_up <- na.omit(up[, 6])
  indata_up <- as.data.frame(indata_up)
  
  #hypergeometric test
  hypr <- ko.enrich.strain(indata_up, strain[i])
  hypr[, 6] <- "up"
  colnames(hypr)[6] <- "updown"
  
  indata_down <- na.omit(down[, 6])
  indata_down <- as.data.frame(indata_down)
  
  hypr2 <-ko.enrich.strain(indata_down, strain[i])
  hypr2[, 6] <- "down"
  colnames(hypr2)[6] <- "updown"
  
  out <- rbind(hypr, hypr2)
  
  write.table(out, file=paste(out_dir, "DEG_KO_enrichment_", strain[i],"_category3.txt", sep = ""), row.names=F,  sep="\t", quote=F)
  
}
