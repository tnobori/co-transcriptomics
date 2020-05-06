#ko.plot() is a function for generating boxplots showing expression fold changes of 
#commensal genes annotated with a given set of KOs.

#created by Tatsuya Nobori
#tnobori@salk.edu

#load config files
source("/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/scripts/_config.R")

#test data
target <- c("K10439", "K10440", "K10441")

#input is a vector of KOs
ko.plot <- function(target){
  target <- unique(target)
  
  dir_data <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/rnaseq_summary_each_strain/"
  dir_database <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A344_commensal_KO_enrichment_analysis/database/"
  
  strain2 <- c("Leaf1","Soil763","Leaf130", "Leaf155","Leaf177", "Leaf187", "Leaf176", "Leaf404", "Root935",   "D36E", "MM", "Pto"  )
  strain <- c("D36E", "Leaf1", "Leaf130", "Leaf155", "Leaf176", "Leaf177", "Leaf187", "Leaf404", "MM", "Pto", "Root935", "Soil763")

  ko_database <- read.delim(file = paste(dir_database, "class-pathway-ko-description.txt", sep = ""),  header = F)
  
  out_fc <- data.frame(matrix(vector(), 0, length(strain),
                       dimnames = list(c(), strain)),
                       stringsAsFactors = F)
  
  out_deg <- data.frame(matrix(vector(), 0, length(strain),
                              dimnames = list(c(), strain)),
                       stringsAsFactors = F)
  
  
  
  for (i in c(1:length(strain))){
    data1 <- read.delim(file = paste(dir_data, "A346_rnaseq_summary_", strain[i],".txt", sep = ""),  header = T)
    fc_target1 <- data1[grep(paste(target, collapse = "|"), data1[,7]), ]
    
    a <- fc_target1
    
    if (nrow(fc_target1)>0){
      # fc_target1[a[, 5]=="DEG" & a[, 4]>0 , 5] <- "up"
      # fc_target1[a[, 5]=="DEG" & a[, 4]<0 , 5] <- "down"
      fc_target1 <- fc_target1[order(fc_target1[,4]), ]
      fc <- fc_target1[,4]
      deg <- fc_target1[,5] %>% as.character()
      out_fc[c(1:length(fc)), i] <- fc
      out_deg[c(1:length(fc)), i] <- deg
      
    }
    
  }
  
  out_fc <- out_fc %>% as.matrix()
  out_deg <- out_deg %>% as.matrix()

  z_fc <- melt(out_fc)
  z_deg <- melt(out_deg)
  
  
  #====#====#====#====#====#====#====#====
  #plot
  #====#====#====#====#====#====#====#====
  
  #plotting options
  optns <- theme(
    axis.text.x = element_text(angle = 270,hjust = 0, vjust = 0.5,margin = margin(c(20, 0, 0, 0)), size = 50, face = "bold", family = "Helvetica"),
    axis.text.y = element_text(margin = margin(c(0, 20, 0, 0)) ,size = 50, face = "bold", family = "Helvetica"), 
    axis.title = element_text(size = 50, face = "bold", family = "Helvetica"),
    axis.title.x = element_text(margin = margin(c(0, 0, 0, 0))),
    axis.title.y = element_text(margin = margin(c(0, 0, 0, 0))),
    axis.ticks.y = element_line(size = 3),
    axis.ticks.x = element_line(size = 3),
    axis.ticks.length = unit(.5, "cm"),
    axis.line  = element_line(size = 3),
    panel.background = element_blank(),
    legend.text = element_text(size = 40, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
    legend.title = element_text(size = 40, face = "bold", family = "Helvetica"),
    plot.title = element_text(hjust = 0.5, size = 20),
    plot.margin = margin(20, 0, 0, 0)
  ) 
  
  
  #making strain-color list
  commensal <- read.delim(file = "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A338_plant_RNA-seq_ALL_combined/commensal_strains.txt",  header = F, row.names = 1)
  strains <- c("Leaf404", "Leaf130", "Leaf155", "Leaf177", "Root935" ,"Leaf176", "Leaf187","Leaf1","Soil763", "Pto", "D36E","MM") %>% as.data.frame()
  rownames(strains) <- strains[,1]
  strains[c(11:12), 1] <- "Pto"
  strains[,1] <- commensal[as.character(strains[,1]), 2]
  names(strains) <- "strains"
  list_strains <- list(strains = c(Actinobacteria ="firebrick3",
                                   Bacteroidetes = "blue",
                                   Firmicutes = "orange",
                                   Alphaproteobacteria = "darkolivegreen4",
                                   Betaproteobacteria = "darkolivegreen1",
                                   Gammaproteobacteria = "darkgreen" )) %>%data.frame()
  
  #for coloring boxes
  strain_exist <- z_fc[!is.na(z_fc[, 3]), 2] %>% as.character()
  strain_melt <-z_fc[, 2] %>% as.character()
  phyla_melt <- strains[strain_melt , 1]%>% as.factor()
  phyla_exist <- strains[strain_exist , 1]%>% as.factor()
  
  #etermine colors for phyla, in the order appears in the plot (alphabetic)
  strain_color2 <- list_strains[as.character(sort(unique(phyla_exist))), 1] %>% as.character()
  
  #for coloring DEGs
  deg_color <- z_deg[,3]
  deg_color <- gsub("nonDEG", "black", deg_color)
  deg_color <- gsub("down", "deepskyblue3", deg_color)
  deg_color <- gsub("up", "deeppink", deg_color)%>% na.omit() %>% as.character()
  
  #plot
  q <- ggplot(z_fc, aes(x = Var2, y = value,  color= phyla_melt)) + 
    geom_boxplot(outlier.color = NA, size = 1.5, fill = NA) + 
    scale_fill_grey() + 
    theme_classic() + 
    optns + 
    scale_color_manual(values = unique(strain_color2)) + 
    geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', position = position_jitter(.1), 
                 binwidth = .1 , alpha = 1, dotsize = 2,  fill = deg_color, color = NA) + 
    xlab("") + ylab("") + theme(legend.position = "none") + 
    scale_x_discrete(limits = strain2) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = ".")) #use the same n of decimal 
  
  q
  
}

#test plot
ko.plot(target)
