############################################################### CODE ###################################################

library(dplyr)
library(purrr) 
#library(table1) 
library(tidyr)
library(magrittr)
#library(sm) 
library(ggplot2)
library(corrplot)
##### the purpose of this file is to determine add'l info about the new users for data in the text #####

#library(caret) 
library(vegan)
library(grid)
library(cowplot)
library(viridis) 
library(nlme)
library(stringr) 
library(pROC)
library(yardstick)
library(reshape2)  

##### STEP 1: Identify changes in metabolites before and after adminstration ####

#read in "light file"  
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt",header=T,sep="\t") 

#read in metadata file 
meta_data <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t")
meta_data$SampleID <- as.character(meta_data$SampleID)
cols <- c("oral5asa_lc","diagnosis","biologics","site_name","Antibiotics","bowelsurg","Participant.ID","steroids","pr5asa","bond5asa","sex","dysbiosis")
meta_data[cols] <- lapply(meta_data[cols], factor) #converts these to a factor variable)
meta_data <- meta_data %>% select(c("SampleID","diagnosis","Participant.ID","dysbiosis","site_name","visit_num","oral5asa_lc","sex","consent_age","Antibiotics","biologics","steroids")) %>% filter(diagnosis == "CD" | diagnosis == "UC")  

mesal_metab<- inner_join(metab_select,meta_data,by="SampleID") %>% drop_na(diagnosis)
   	#write.table(mesal_metab, "newusers_extra.txt", sep="\t", quote=FALSE,row.names=FALSE) #this file helps us find users who stopped 5-asa and who started biologics shown below(based on binary154)

########################### NEW BIOLOGICS #################################

#I then selected participants who were new users of biologics based on the file above
newbiologicids <-c("CSM67UB1",
"CSM79HIT",
"CSM79HKZ",
"CSM7KOKZ",
"ESM5MEDD",
"ESM718UH",
"HSM5MD73",
"HSM7CYZ7",
"HSM6XRQM",
"HSM7CYY3",
"HSM67VD6",
"HSM67VHF",
"HSM7J4LP",
"HSMA33NQ",
"MSM6J2OH",
"MSM6J2OP",
"MSM5LLF4",
"MSM6J2LL")

#subset
new_bio_iddata <- mesal_metab %>% filter(SampleID %in% newbiologicids)
new_bio_iddata%>% arrange(Participant.ID,binary154) %>% select(biologics,Participant.ID,visit_num,SampleID,X154.0502_3.83, binary154) #this shows the list n=9
   biologics Participant.ID visit_num SampleID X154.0502_3.83list n=9
1          0          C3016         8 CSM67UB1         268865
2          1          C3016        13 CSM79HIT          52320
3          0          C3028         6 CSM79HKZ          65062
4          1          C3028        15 CSM7KOKZ         200779
5          0          E5009         9 ESM5MEDD         122906
6          1          E5009        13 ESM718UH         211335
7          0          H4015         8 HSM5MD73         168868
8          1          H4015        19 HSM7CYZ7         275618
9          0          H4017         8 HSM6XRQM         127787
10         1          H4017        19 HSM7CYY3         252851
11         0          H4020         9 HSM67VD6         154486
12         1          H4020        13 HSM67VHF         160714
13         0          H4039         4 HSM7J4LP         152763
14         1          H4039        14 HSMA33NQ         367953
15         0          M2026        19 MSM6J2OH         189830
16         1          M2026        27 MSM6J2OP         319461
17         0          M2034         8 MSM5LLF4         155485
18         1          M2034        14 MSM6J2LL         211623

metab_rotated2 <- read.table("/home/rsm34/5asa/metab_data_nozero2.txt",header=T,sep="\t")

metabs_newids_bio <- inner_join(new_bio_iddata,metab_rotated2,by="SampleID") %>% arrange(Participant.ID)%>% arrange(biologics) %>% select(c("SampleID","nicotinuric.acid","nicotinate","Participant.ID","visit_num","biologics"))

wilcox.test(metabs_newids_bio$nicotinate ~ metabs_newids_bio$biologics, paired=TRUE) #p=0.91 (nicotinuric acid is undetectable in biologic users) 

library(ggpubr)
mesallabels <- c("pre", "post")

P2 <- ggpaired(metabs_newids_bio , x = "biologics", y = "nicotinate",
         color = "biologics", line.color = "gray", line.size = 0.4,
         palette = c("#3C5488FF", "#DC0000FF"),short.panel.labs = FALSE, xlab="Biologic use",ylab="rel.abund.",title = "Nicotinic Acid",point.size = 1.5,legend="none") +theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(labels=mesallabels)

ggsave(filename='/home/rsm34/5asabb3/reviews/nicotinbiologic.png', plot=P2, width = 3, height = 3, dpi = 600) 


########################### NEW STEROIDS #################################

#I then selected participants who were new users of biologics based on the file above
newsteroidids <-c(
"CSM79HKZ",
"CSM7KOKZ",
"ESM5MEDD",
"ESM718UH",
"HSM5MD4O",
"HSM6XRSN",
"HSM5MD6A",
"HSM6XRRV",
"HSM5MD73",
"HSM6XRS8",
"HSM6XRRD",
"HSM6XRV2",
"HSM7J4NE",
"HSMA33O1",
"HSM7J4K8",
"HSMA33OR",
"PSM6XBQY",
"PSM6XBSS",
"PSM6XBSK",
"PSM6XBUG",
"PSM7J17L",
"PSM7J158")

#subset
new_ster_iddata <- mesal_metab %>% filter(SampleID %in% newsteroidids)
new_ster_iddata<- new_ster_iddata%>% arrange(Participant.ID,binary154) %>% select(steroids,Participant.ID,visit_num,SampleID,X154.0502_3.83, binary154) %>% filter(binary154==0) %>% filter(!Participant.ID=="H4015") %>% filter(!Participant.ID=="P6010") #this shows the list n=3
  steroids Participant.ID visit_num SampleID X154.0502_3.83 binary154ipant.ID=="H4015") %>% filter(!Participan1        0          C3028         6 CSM79HKZ          65062         0
2        1          C3028        15 CSM7KOKZ         200779         0
3        0          E5009         9 ESM5MEDD         122906         0
4        1          E5009        13 ESM718UH         211335         0
5        0          H4014         8 HSM5MD6A         508022         0
6        1          H4014        13 HSM6XRRV         290011         0


metab_rotated2 <- read.table("/home/rsm34/5asa/metab_data_nozero2.txt",header=T,sep="\t")

metabs_newids_ster <- inner_join(new_ster_iddata,metab_rotated2,by="SampleID") %>% arrange(Participant.ID)%>% arrange(steroids) %>% select(c("SampleID","nicotinuric.acid","nicotinate","Participant.ID","visit_num","steroids","X2.aminoadipate"))

wilcox.test(metabs_newids_ster$nicotinate ~ metabs_newids_ster$steroids, paired=TRUE) #p=0.64 (nicotinuric acid is undetectable in biologic users) 
wilcox.test(metabs_newids_ster$X2.aminoadipate ~ metabs_newids_ster$steroids, paired=TRUE) #p=0.24 (nicotinuric acid is undetectable in biologic users) 

P3 <- ggpaired(metabs_newids_ster , x = "steroids", y = "nicotinate",
         color = "steroids", line.color = "gray", line.size = 0.4,
         palette = c("#3C5488FF", "#DC0000FF"),short.panel.labs = FALSE, xlab="Steroid use",ylab="rel.abund.",title = "Nicotinic Acid",point.size = 1.5,legend="none") +theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(labels=mesallabels)

ggsave(filename='/home/rsm34/5asabb3/reviews/nicotinsteroid.png', plot=P3, width = 3, height = 3, dpi = 600) 


