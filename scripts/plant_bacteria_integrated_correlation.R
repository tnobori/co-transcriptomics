#This script is to generate the correlation heatmap in Fig. 6 of Nobori et al., 2020, bioRxiv.

#created by Tatsuya Nobori
#tnobori@salk.edu

rm(list=ls())

#load config 
source("/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/scripts/_config.R")

#=========#=========#=========#=========#=========#=========#=========#=========
#=========#=========#=========#=========#=========#=========#=========#=========
#strains to use
strains <- c("Leaf404", "Leaf130", "Leaf155", "Leaf177", "Root935" ,"Leaf176", "Leaf187","Leaf1","Soil763")
#=========#=========#=========#=========#=========#=========#=========#=========
#=========#=========#=========#=========#=========#=========#=========#=========

dir_data <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A310_psl_commensal_rnaseq_analysis/A310_OG_analysis_3/A310_OG_combined/"
out_dir <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A347_commensal_PLANT_rnaseq/out/"

# bacterial OG expression file
bac <- read.delim(paste(dir_data, "A310_tmm_OG-KO_FC_strain_combined.txt", sep = "" ), header = T, row.names = 1)
bac <- bac[, -c(10:13)]
bac2 <- na.omit(bac)

#plant RNA-seq data; genes that are DEG at least in one of the strains
path_plant <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A347_commensal_PLANT_rnaseq/out/A347_fitted_mean_ave1_DEG_atleast_one_strain.txt"
plant <- read.delim(path_plant, header = T, row.names = 1)


a <- grep("mock", colnames(plant)) #mock
plant2 <- sweep(plant, 2, plant[, a]) #fold changes 
plant2 <- plant2[, -a]
plant2 <- plant2[, colnames(bac)]

#============#============#============#============#============
#function to calculate correlation without one strain
make_corr <- function(strain){
  #selected strains
  plant3 <- plant2[, strains[-grep(strain, strains)]]
  bac3 <- bac2[, strains[-grep(strain, strains)]]
  
  out1 <- matrix(data=NA, nrow = nrow(bac3), ncol = nrow(plant3),
                 dimnames = list(rownames(bac3), rownames(plant3)))
  
  for (j in c(1:nrow(plant3))){
    x <- apply(bac3, 1, function(x){cor(x,as.numeric(plant3[j, ]) )})
    out1[, j] <- x
  }
  
  out1  
}

#calculating correlation 
message("...")
out1 <- make_corr(strains[1])
message("...")
out2 <- make_corr(strains[2])
message("...")
out3 <- make_corr(strains[3])
message("...")
out4 <- make_corr(strains[4])
message("...")
out5 <- make_corr(strains[5])
message("...")
out6 <- make_corr(strains[6])
message("...")
out7 <- make_corr(strains[7])
message("...")
out8 <- make_corr(strains[8])
message("...")
out9 <- make_corr(strains[9])

##+======##+======##+======##+======
#using all strains
plant3 <- plant2
bac3 <- bac2

out <- matrix(data=NA, nrow = nrow(bac3), ncol = nrow(plant3),
              dimnames = list(rownames(bac3), rownames(plant3)))


for (j in c(1:nrow(plant3))){
  x <- apply(bac3, 1, function(x){cor(x,as.numeric(plant3[j, ]) )})
  out[, j] <- x
}

##+======##+======##+======##+======
#for each pair of plant-bacterial genes, take the smallest correlation coefficients
out_integrated <- out
out_integrated[] <- 0

for (i in c(1:nrow(out))){
  for (j in c(1:ncol(out))){
    if (out[i,j]>0 & out1[i,j]>0 & out2[i,j]>0 & out3[i,j]>0 & out4[i,j]>0 & out5[i,j]>0 & out6[i,j]>0 & out7[i,j]>0 & out8[i,j]>0 & out9[i,j]>0){
      out_integrated[i, j] <- min(out[i,j], out1[i,j], out2[i,j], out3[i,j], out4[i,j], out5[i,j], out6[i,j], out7[i,j], out8[i,j], out9[i,j])
    } else if (out[i,j]<0 & out1[i,j]<0 & out2[i,j]<0 & out3[i,j]<0 & out4[i,j]<0 & out5[i,j]<0 & out6[i,j]<0 & out7[i,j]<0 & out8[i,j]<0 & out9[i,j]<0){
      out_integrated[i, j] <- max(out[i,j], out1[i,j], out2[i,j], out3[i,j], out4[i,j], out5[i,j], out6[i,j], out7[i,j], out8[i,j], out9[i,j])
    } else{
      out_integrated[i ,j] <- 0
    }
  }
}

write.table(out_integrated, file=paste(out_dir, "A310_plantDEG_bac_corr_bootstrap.txt", sep = ""), row.names=T, col.names=NA, sep="\t", quote=F)

