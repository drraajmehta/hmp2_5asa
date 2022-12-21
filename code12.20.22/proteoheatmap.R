$grep 'Escherichia\|freundii\|amalonaticus\|_pneumoniae\|oxytoca\|morganii' /home/rsm34/5asabb3/archive/DNA_uniref/dna_merge_genefams.txt > protomgx.txt
$grep 'Escherichia\|freundii\|amalonaticus\|_pneumoniae\|oxytoca\|morganii' /home/rsm34/5asabb3/archive/RNA_uniref/rna_merge_genefams.txt > protomtx.txt



#### MGX in R ####

library(dplyr)
library(tidyr)

sixbugs <- read.table("/home/rsm34/5asabb3/reviews/protomgx.txt",header=F,sep="\t")
names(sixbugs )[names(sixbugs ) == 'V1'] <- 'uniref'
sixbugs$species <- gsub(".*\\|","",sixbugs$uniref)
sixbugs$uniref <- gsub("\\|.*","",sixbugs$uniref) 
 
#read in filtered uniref90 names (10e-8 etc filtration) 
filtuniref <- read.table("/home/rsm34/5asabb3/reviews/uniref90filtMGXnames.txt",header=T,sep="\t")
colnames(filtuniref) <- "uniref"


dnafiltbugs <- inner_join(sixbugs,filtuniref, by="uniref")

#### MTX in R ####

library(dplyr)
library(tidyr)

sixbugs <- read.table("/home/rsm34/5asabb3/reviews/protomtx.txt",header=F,sep="\t")
names(sixbugs )[names(sixbugs ) == 'V1'] <- 'uniref'
sixbugs$species <- gsub(".*\\|","",sixbugs$uniref)
sixbugs$uniref <- gsub("\\|.*","",sixbugs$uniref) 
 
#read in filtered uniref90 names (10e-8 etc filtration) 
filtuniref <- read.table("/home/rsm34/5asabb3/reviews/uniref90filtMTXnames.txt",header=T,sep="\t")
colnames(filtuniref) <- "uniref"


rnafiltbugs <- inner_join(sixbugs,filtuniref, by="uniref")

writeout <- rnafiltbugs %>% select(uniref) 
writeout$uniref <- gsub("UniRef90_","",writeout$uniref)


write.table(writeout, file="rnasixbugsnames.txt", sep="\t", quote=FALSE,row.names=FALSE,col.names=FALSE)

#there are no arylamine acetyltransferases in this selection at the mTX level, and there are none in the entire dataset for that matter.  


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ grep the species and make a proportion lollipop plot $$$$$$$$$$$$$$$$$$$$$$$#

$grep 'Escherichia_c\|freundii\|amalonaticus\|_pneumoniae\|oxytoca\|morganii\|bacterium_CAG_83\|sp_CAG_241\|bacterium_D16\|prausnitzii' /home/rsm34/5asabb3/archive/taxa/metaphlan_species_raw.txt > prototaxa.txt
$ grep 'R6CZ24\|R5CY66\|T5S060\|C7H1G6' /home/rsm34/5asabb3/archive/DNA_uniref/dna_merge_genefams.txt > bugsthatencodeuniref.txt

library(dplyr)
library(ggplot2) 
prototaxa <- read.table("/home/rsm34/5asabb3/reviews/prototaxa.txt",header=F,sep="\t")

prototaxa$prevalence <- (rowSums(prototaxa!=0)/ncol(prototaxa))*100

prototaxa$species <- gsub(".*s__","",prototaxa$V1) 

prototaxa <- prototaxa[5:10,]
 
protospec <-ggplot(prototaxa, aes(x=species, y=prevalence)) +
  geom_segment( aes(x=species, xend=species, y=0, yend=prevalence), color="dark green") +
  geom_point( color="#90EE90", size=4,alpha=0.8) +
  theme_light() +
  coord_flip() +
  theme_bw() 


#[ ] change name to Firm Cag 
#[ ] color by ours vs green
#[ ] sort by descending 

ggsave(filename='./protospec.png', plot=protospec , width = 6, height = 4, dpi = 600)


