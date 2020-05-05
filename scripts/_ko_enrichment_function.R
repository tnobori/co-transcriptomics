

#loading the config file 
source("/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/scripts/_config.R")

####test data
test_ko_list <- read.delim("/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A344_commensal_KO_enrichment_analysis/test_ko_list/test_ko_list.txt",  header = T)

#========#========#========#========#========#========#========#========#========#========#========
#ko enrichment function
#input: list of KOs as a data frame & strain ID
#Enrichment analysis is performed for KO category 3
#========#========#========#========#========#========#========#========#========#========#========
ko.enrich.strain <- function(test_ko_list, strain){
  
  #loading a KO database for the target strain
  dir <- "/Users/tatsuyanobori/Desktop/MPIPZ/Projects/Bacterial_Transcriptome/A346_commensal_rnaseq_for_paper/ko_database/"
  ko_data <- read.delim(file = paste(dir, "KO_database_", strain,".txt", sep = ""),  header = T)
  
  f3 <- unique(ko_data[, 4]) %>% as.character() #all KO terms in the category 3
  f_source3 <- matrix(NA, nrow = length(f3), ncol = 1, dimnames = list(f3, "Ngenes"))
  
  #summarize the number of each KO term in the database
  for (i in c(1:length(f3))){
    count.f <- ko_data[ko_data[, 4]==f3[i], 6] %>% sum()
    f_source3[i, 1] <- count.f
  }
  
  f_source3 <- f_source3 %>% as.data.frame() %>% rownames_to_column()
  
  #summarize the number of each KO term in the input data 
  f_input <- matrix(NA, nrow = 1, ncol = ncol(ko_data), dimnames = list("", colnames(ko_data)))
  
  for (i in c(1:nrow(test_ko_list))){
    x <- ko_data[ko_data[, 1] == as.character(test_ko_list[i, ]), ]
    
    f_input <- rbind(f_input, x)
  }
  
  f_input <- f_input[-1, ]
  
  f_in_table <- table(f_input[, 4]) %>% as.data.frame()
  
  
  #hypergeometric test
  a <- f_source3
  b <- f_in_table
  
  out_hypr <- b 
  out_hypr[, 2] <- NA 
  
  for (i in (1:nrow(b))){
    q <- b[i, 2] #n of white balls drawn
    m <- a[b[i, 1], 2] # n of white balls in the box
    n <- sum(a[,2])-m #n of black balls in the box
    k <- nrow(f_input) #n of ball drawn
    
    pval <- phyper(q, m, n, k, log.p = FALSE, lower.tail = FALSE)
    fold_enrich <- (q/k)/(m/(m+n)) 
    
    out_hypr[i, 2] <- pval
    out_hypr[i, 3] <- fold_enrich
    
  }
  
  out_hypr[, 2] <- p.adjust(out_hypr[, 2], method = "fdr")
  
  #select significantly enriched KO terms
  out_significant <- out_hypr[out_hypr[, 2]<0.01 & out_hypr[, 3] > 2 ,]
  colnames(out_significant) <- c("KO", "pval", "FoldEnrichment")
  
  #summerize genes involved in significant KOs 
  for (i in c(1:nrow(out_significant))){
    z <- f_input[f_input[, 4] == out_significant[i, 1], c(1,5)]
    zz <- paste(z[, 1], z[,2], sep = "---")
    g4paste <- paste(zz, collapse = "|")
    out_significant[i, 4] <- g4paste
    out_significant[i, 5] <- nrow(z)
  }
  colnames(out_significant)[c(4, 5)] <- c("genes", "size")
  out_significant2 <- out_significant[out_significant[, 5] >4, ]
  out_significant2  
}




