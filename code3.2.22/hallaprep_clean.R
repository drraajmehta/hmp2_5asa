#### PURPOSE: GENERATE FILES FOR HALLA (RUN IN RC / FAS ) ###############

######## code #############

library(dplyr)
library(purrr)
library(tidyr)
library(magrittr)
library(sm)
library(ggplot2)
library(stringr) 

#prepare metabolite data
##########################
#read in the 683 selected metabolites with AUC >95, new users, subset to users only (done in ROCmetabs_clean) 
metab_select <- read.table("/home/rsm34/5asabb3/intermediate_files/metab_halla_input_pre.txt", header=T, sep="\t")  

#read in max of each metabolite group (gets of rid of co-correlated things) 

grouped <- read.table("/home/rsm34/5asabb3/intermediate_files/groupedmetabnames_halla.txt",header=T,sep="\t") %>% select(names)
groupednames <- rbind(c("SampleID"),grouped ) 


metab_select <- metab_select %>% select(any_of(groupednames$names))
 
dim(metab_select)
head(names(metab_select),20)
metab_halla_pre <- as.data.frame(t(metab_select))


#prepare species data
##########################
#read in the SPECIES file after pre-processing script made it (1632 --> n=1464)
species_rotated <- read.table("/home/rsm34/5asabb3/archive/taxa/metaphlan_rotated.txt",header=T,sep="\t")
species_rotated<-species_rotated[rowSums(species_rotated[,2:ncol(species_rotated)])!=0,,drop=FALSE] 	#there are a few rows that sum to zero (also a few that sum to <0.5
IDlist <- metab_select %>% select(SampleID)
species_subset <- inner_join(species_rotated,IDlist,by="SampleID")
dim(species_subset)
specnames <- names(species_subset)
specnameshort <- gsub(".*s__","",specnames)
names(species_subset) <- specnameshort


spec_ibd_nlme <- read.table("/home/rsm34/5asabb3/nlmeoutput/spec_ibd_nlme.txt",header=T,sep="\t") 
spec_ibd_nlme <- spec_ibd_nlme %>% select(!(feature))
spec_ibd_nlme2 <- spec_ibd_nlme %>% filter(FDR < 0.25) %>% rename(feature=names) %>% select(feature) 

spec_names <- rbind(c("SampleID"),spec_ibd_nlme2) 


species_subset2 <-  species_subset %>% select(any_of(spec_names$feature))

IDlist2 <- species_subset2 %>% select(SampleID) #used for later 
species_halla <- as.data.frame(t(species_subset2))
dim(species_halla)

#write out tables
##################

#species
#######
write.table(species_halla,"species_halla.txt",sep="\t",quote=FALSE,col.names=FALSE)

#metabs
#########
#write metab file that matches 
metab_select2 <- inner_join(IDlist2, metab_select,by="SampleID")
dim(metab_select2)
metab_halla <- as.data.frame(t(metab_select2))
write.table(metab_halla,"metab_halla.txt",sep="\t",quote=FALSE,col.names=F)

