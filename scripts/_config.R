library(tidyverse)
library(pheatmap)
require(RColorBrewer)
require(reshape2)
require(pheatmap)
require(scales)
require(viridis)
library(magrittr)
library(edgeR)
library(vegan)
library(labdsv)
library(splitstackshape)
library(rlist)
library(qvalue)
library(ggrepel)
#====#====#====#====#====#====#====#====
#plot
#====#====#====#====#====#====#====#====

optns <- theme(
  axis.text.x = element_text(angle = 270,hjust = 1, vjust = 0.5,margin = margin(c(20, 0, 0, 0)), size = 20, face = "bold", family = "Helvetica"),
  axis.text.y = element_text(margin = margin(c(0, 20, 0, 0)) ,size = 20, face = "bold", family = "Helvetica"), 
  axis.title = element_text(size = 20, face = "bold", family = "Helvetica"),
  axis.title.x = element_text(margin = margin(c(30, 0, 0, 0))),
  axis.title.y = element_text(margin = margin(c(0, 30, 0, 0))),
  axis.ticks.y = element_line(size = 3),
  axis.ticks.x = element_line(size = 3),
  axis.ticks.length = unit(.5, "cm"),
  axis.line  = element_line(size = 3),
  panel.background = element_blank(),
  legend.text = element_text(size = 20, face = "bold", family = "Helvetica", margin = margin(c(20, 0, 20, 0))),
  legend.title = element_text(size = 20, face = "bold", family = "Helvetica"),
  plot.title = element_text(hjust = 0.5, size = 20)
) 


#====#====#====#====#====#====#====#====
#clustering functions
#====#====#====#====#====#====#====#====


kplot <- function(input){
  wss <- function(k) {
    kmeans(input, k, nstart = 10 )$tot.withinss
  }
  k.values <- 1:30
  wss_values <- map_dbl(k.values, wss)
  plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE,
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
}

kclust <- function(input, n){
  km<- kmeans(input,n)
  m.kmeans<- cbind(input, km$cluster)
  colnames(m.kmeans)[ncol(m.kmeans)] <- "kmean"
  m.kmeans <- as.data.frame(m.kmeans)
  m.kmeans <- tibble::rownames_to_column(m.kmeans)
  o <-
    m.kmeans %>%
    arrange(kmean)
  
  data_km <- column_to_rownames(o)
  data_km
}

