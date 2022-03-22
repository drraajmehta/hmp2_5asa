library(dplyr)

#bring in species data
########################
#read in the SPECIES file after pre-processing script made it (1632 --> n=1464)
species_rotated <- read.table("/home/rsm34/5asabb3/archive/taxa/metaphlan_rotated.txt",header=F,sep="\t")
graphlanfile <- as.data.frame(t(species_rotated[1:2,]))
rownames(graphlanfile)=NULL
names(graphlanfile) <- c("taxa","no")
graphlanfile<- graphlanfile[-1,]

graphlanfile$taxa <- gsub(".*p__","",graphlanfile$taxa)
graphlanfile$taxa <- gsub("c__","",graphlanfile$taxa)
graphlanfile$taxa <- gsub("o__","",graphlanfile$taxa)
graphlanfile$taxa <- gsub("f__","",graphlanfile$taxa)
graphlanfile$taxa <- gsub("g__","",graphlanfile$taxa)
graphlanfile$taxa <- gsub("s__","",graphlanfile$taxa)
graphlanfile$taxa <- gsub("\\|",".",graphlanfile$taxa)

graphlanfile <- graphlanfile %>% select(1) %>% arrange(taxa)

write.table(graphlanfile, file="graphlanfile.txt", sep="\t", quote=FALSE,row.names=FALSE,col.names = F)


