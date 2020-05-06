#This script is for generating boxplots of commensal genes with various functional (KEGG) annotations
#shown in Fig. 2, 4, S5, S6, S9, S10, S11, and S12 in Nobori et al., 2020, bioRxiv

#created by Tatsuya Nobori
#tnobori@salk.edu

rm(list=ls())

#load config 
source("/scripts/_config.R")
source("/scripts/_ko_plot_function.R")

#
strain <- c("Leaf1","Leaf130", "Leaf155","Leaf177", "Leaf187", "Leaf176", "Leaf404", "Root935", "Soil763", "Pto", "D36E", "MM")

#combined file directory 
dir_database <- "/data/ko_database_for_each_strain/"
out_dir <- "/output/directory/"

#loading the KO database
ko_database <- read.delim(file = paste(dir_database, "_KO_master_database.txt", sep = ""),  header = F)

#search for KO pathway annotations
infile <- read.delim(file = paste(dir_database, "ko_pathway_annotation_list_for_plots.txt", sep = ""),  header = T)
ko2 <- infile[, 1] %>% unique()

for (i in c(1:length(ko2))) {
  l <- ko_database[grep(ko2[i], ko_database[, 4]), 1] %>% unique() %>% as.character()
  if (length(l)>0){
    q <- ko.plot(l)
  }
  ggsave(paste(out_dir, "Expression_boxplot_categ3_", ko2[i],".png"),  height = 10, width = 13, q)

}


#keyword search for KO BRITE annotations
target <- c("manganese",
            "mannose transport|mannose-specific",
            "fructose transp",
            "sucrose-specific|scrY",
            "arabinose transport",
            " galactonate ",
            " malate ",
            "malonate ",
            "xylose transport",
            "ubiquinol oxydase|ubiquinol-cytochrome",
            "flagellar",
            "chemotaxis",
            "accepting chemotaxis",
            "chemotaxis protein Mot",
            "urease",
            "DNA polymerase III",
            "DNA gyrase",
            "DNA-directed RNA polymerase",
            "glpK",
            "3-oxoacyl-\\[acyl-carrier protein] reductase",
            "polysaccharide biosynthesis",
            "catalase",
            "osmoprotectant",
            " glutamate synthase",
            "polysaccharide biosynthesis",
            "iron complex outermembrane",
            "zinc transport system",
            "putative tricarboxylic transport",
            "aquaporin",
            "pyoverdin transport",
            "stringent starvation protein", " cyanate transporter",
            "glycerol-3-phosphate transporter",
            "heme exporter protein",
            "cobalt-zinc-cadmium resistance",
            "Ca-activated chloride channel",
            "methionine transport system",
            "glycerol transport system",
            "nitrite transport system",
            "methylmethionine transporter",
            "aspartate transport system",
            "biopolymer transport ",
            "molybdate transport",
            "alpha-glucoside transport",
            "nickel transport",
            "multiple sugar transport",
            "hexuronate transporter",
            "enterobactin",
            "starch-binding outer membrane",
            "phosphate transport", 
            "homoserine lactone", 
            "preprotein translocase subunit",
            "PTS system",
            "iron transport|iron complex transport|iron\\(III) transport system",
            "lipopolysaccharide transport system",
            "starch\\-binding",
            "multidrug efflux system","type IV secretion", 
            "type VI secretion", 
            "flagellar", 
            "type III secretion",
            "preprotein translocase subunit", 
            "general secretion pathway protein",
            "large subunit ribosomal protein", 
            "alginate", 
            "ribose transport system", 
            "sulfonate transport system", 
            "urea transport system",
            "two-component system, OmpR family", 
            "two-component system, NarL family", 
            "two-component system, CitB family", 
            "two-component system, LytTR family", 
            "two-component system, NtrC family",
            "multidrug efflux system",
            "nuo",
            "ATPF")

for (i in c(1:length(target))){
  l <- ko_database[grep(target[i], ko_database[, 5]), 1] %>% unique() %>% as.character()
  
  q <- ko.plot(l)
  
  ggsave(paste(out_dir, "Expression_boxplot_", target[i],".png"),  height = 10, width = 13, q)
}


#for sulfur transport-related genes (Fig. 4C)
target <- c("K15553|K15554|K15555|K02047|K02048")
q <- ko.plot(target)
ggsave(paste(out_dir, "Expression_boxplot_sulfur_transport.png"),  height = 10, width = 13, q)




