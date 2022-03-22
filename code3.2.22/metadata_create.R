############################################################### CODE ###################################################

library(dplyr)
library(purrr)
library(table1) 
library(tidyr)
library(magrittr)
library(sm)
library(caret) 
library(grid)

################################################################
#					#
#      PRE-PROCESSING: Generating the metadata file	#
#					#
################################################################


#read in file (this was lightly edited for ease of import into R) 
metadata <- read.csv("/home/rsm34/5asabb3/archive/metadata/hmp2_metadata_load.csv",stringsAsFactors=FALSE)
#make sure there is no difference between this and the other file 
#$comm -2 -3 <(sort /home/rsm34/5asabb3/archive/metadata/hmp2_metadata.csv) <(sort /home/rsm34/5asabb3/archive/metadata/hmp2_metadata_load.csv) > file3.txt

#select columns of interest
metadata2 <- metadata %>% select(1:10,"diagnosis","site_name","consent_age","Age.at.diagnosis",93,94,99,100,188,297,298,318,320,322,354,373,375,377,379,381,383,431,478, 334,337,342,348,350,352,358,361,366,369)

#divide up data into chunks
###########################
meta_WGS <- metadata2 %>% filter(data_type=="metagenomics")  ## this will be the master
meta_serology <- metadata2 %>% filter(data_type=="serology")  #has drug data
meta_hostgene <- metadata2 %>% filter(data_type=="host_genome") #has drug data
meta_metabolome <- metadata2 %>% filter(data_type=="metabolomics")
meta_metatranscript <- metadata2 %>% filter(data_type=="metatranscriptomics")  
meta_methylome <- metadata2 %>% filter(data_type=="methylome") #has drug data  
meta_drugdata <- metadata2 %>% filter(data_type=="methylome"|data_type=="host_genome"|data_type=="serology")

#bin 5asa and other drug use into right categories per the key below
####################################################
#bonded
#Dipentum (olsalazine) 
#Colozal (balasalizide)
#Sulfasalizine (Azulfidine) 

#PR
#Rowasa enemas (mesalamine enemas)
#Canasa suppositories (mesalamine suppositories)
 
#oral
#Asacol (mesalamine)
#Pentasa (mesalamine)
#Lialda (mesalamine)
#Apriso (mesalamine)

#steroids
#entecort (budesonide)
#solumedrol (medrol) 
#prednisone 

#biologics
#Remicade (infliximab)
#Humira (adalimumab)
#Cimzia (certilizumab)
#Tysabri (natalizumab) 

#immunomodulators
#Azathioprine (imuran)
#Methotrexate 
#Mercaptopurine (6MP) 

#current vs. taken since last visit 
###################################
#oral 5asa current or not (strict --> eg includes only "current")  
meta_drugdata <- meta_drugdata %>% 
	mutate(oral5asa_c = ifelse(Asacol..mesalamine. == "Current" | Pentasa..mesalamine. == "Current" | Lialda..mesalamine. == "Current" | Apriso..mesalamine.== "Current",1,
	                  ifelse(Asacol..mesalamine. == "Never taken" | Pentasa..mesalamine. == "Never taken" | Lialda..mesalamine. == "Never taken" | Apriso..mesalamine.== "Never taken"|
			Asacol..mesalamine. == "Taken since last visit" | Pentasa..mesalamine. == "Taken since last visit" | Lialda..mesalamine. == "Taken since last visit" | Apriso..mesalamine.== "Taken since last visit"|
			Asacol..mesalamine. == "Taken prior to baseline" | Pentasa..mesalamine. == "Taken prior to baseline" | Lialda..mesalamine. == "Taken prior to baseline" | Apriso..mesalamine.== "Taken prior to baseline",
			0,"")))

#oral 5asa current or not (loose --> eg includes "taken since last visit")
meta_drugdata <- meta_drugdata %>% 
	mutate(oral5asa_lc = ifelse(Asacol..mesalamine. == "Current" | Pentasa..mesalamine. == "Current" | Lialda..mesalamine. == "Current" | Apriso..mesalamine.== "Current"|
			Asacol..mesalamine. == "Taken since last visit" | Pentasa..mesalamine. == "Taken since last visit" | Lialda..mesalamine. == "Taken since last visit" | Apriso..mesalamine.== "Taken since last visit",1,
	                  ifelse(Asacol..mesalamine. == "Never taken" | Pentasa..mesalamine. == "Never taken" | Lialda..mesalamine. == "Never taken" | Apriso..mesalamine.== "Never taken"|
			Asacol..mesalamine. == "Taken prior to baseline" | Pentasa..mesalamine. == "Taken prior to baseline" | Lialda..mesalamine. == "Taken prior to baseline" | Apriso..mesalamine.== "Taken prior to baseline",
			0,"")))

#PR (loose) 
meta_drugdata <- meta_drugdata %>% 
	mutate(pr5asa = ifelse(Canasa.suppositories..mesalamine.suppositories. == "Current" | Rowasa.enemas..mesalamine.enemas. == "Current"|
		          Canasa.suppositories..mesalamine.suppositories. == "Taken since last visit" | Rowasa.enemas..mesalamine.enemas. == "Taken since last visit",1,
		   ifelse(Canasa.suppositories..mesalamine.suppositories. == "Never taken" | Rowasa.enemas..mesalamine.enemas. == "Never taken"|
		          Canasa.suppositories..mesalamine.suppositories. == "Taken prior to baseline" | Rowasa.enemas..mesalamine.enemas. == "Taken prior to baseline",
			0,"")))

#Bond5asa (loose) 
meta_drugdata <- meta_drugdata %>% 
	mutate(bond5asa  = ifelse(Dipentum..olsalazine. == "Current" | Colozal..balasalizide. == "Current" | Sulfasalizine..Azulfidine. == "Current" |
			Dipentum..olsalazine. == "Taken since last visit" | Colozal..balasalizide. == "Taken since last visit" | Sulfasalizine..Azulfidine. == "Taken since last visit" ,1,
	                  ifelse(Dipentum..olsalazine. == "Never taken" | Colozal..balasalizide. == "Never taken" | Sulfasalizine..Azulfidine. == "Never taken" |
			Dipentum..olsalazine. == "Taken prior to baseline" | Colozal..balasalizide. == "Taken prior to baseline" | Lialda..mesalamine. == "Taken prior to baseline",
			0,"")))

#Any 5asa use (loose)
meta_drugdata <- meta_drugdata %>% mutate(any5asa = ifelse(oral5asa_lc == 1 | pr5asa  == 1 | bond5asa  == 1,1,ifelse(oral5asa_lc == 0 | pr5asa  == 0 | bond5asa  == 0,0,"")))

#Other drugs 
##############
#biologics
meta_drugdata <- meta_drugdata %>% 
		mutate(biologics = ifelse(Remicade..Infliximab. == "Current" | Humira..Adalimumab. == "Current" | Cimzia..Certlizumab. == "Current" | Tysabri..Natalizumab.== "Current"|
				 Remicade..Infliximab. == "Taken since last visit" | Humira..Adalimumab. == "Taken since last visit" | Cimzia..Certlizumab. == "Taken since last visit" | Tysabri..Natalizumab.== "Taken since last visit",1,
			      ifelse(Remicade..Infliximab. == "Never taken" | Humira..Adalimumab. == "Never taken" | Cimzia..Certlizumab. == "Never taken" | Tysabri..Natalizumab.== "Never taken"|
				 Remicade..Infliximab. == "Taken prior to baseline" | Humira..Adalimumab. == "Taken prior to baseline" | Cimzia..Certlizumab. == "Taken prior to baseline" | Tysabri..Natalizumab.== "Taken prior to baseline",0,"")))

#steroids 
meta_drugdata <- meta_drugdata %>% 
		mutate(steroids = ifelse(Prednisone == "Current" | Entocort..Budesonide. == "Current" | Solumedrol..Medrol. == "Current"|
				Prednisone == "Taken since last visit" | Entocort..Budesonide. == "Taken since last visit" | Solumedrol..Medrol. == "Taken since last visit",1,
			     ifelse(Prednisone == "Never taken" | Entocort..Budesonide. == "Never taken" | Solumedrol..Medrol. == "Never taken"|
				Prednisone == "Taken prior to baseline" | Entocort..Budesonide. == "Taken prior to baseline" | Solumedrol..Medrol. == "Taken prior to baseline",0,"")))

#immunomodulators
meta_drugdata <- meta_drugdata %>% 
		mutate(immunomod = ifelse(Azathioprine..Imuran..Azasan. == "Current" | Methotrexate == "Current" | Mercaptopurine..Purinethol..6MP. == "Current" |
				Azathioprine..Imuran..Azasan. == "Taken since last visit" | Methotrexate == "Taken since last visit" | Mercaptopurine..Purinethol..6MP. == "Taken since last visit",1,
			      ifelse(Azathioprine..Imuran..Azasan. == "Never taken" | Methotrexate == "Never taken" | Mercaptopurine..Purinethol..6MP. == "Never taken"|
				Azathioprine..Imuran..Azasan. == "Taken prior to baseline" | Methotrexate == "Taken prior to baseline" | Mercaptopurine..Purinethol..6MP. == "Taken prior to baseline",0,"")))

#clean the 5asa dataset
########################
meta_drugdata$id <- paste(meta_drugdata$Participant.ID,meta_drugdata$visit_num,sep=".")
meta_drugdata2 <- meta_drugdata %>% distinct(id, .keep_all = TRUE) %>% arrange(Participant.ID, visit_num) %>% select(id,44:51,20,21,3,9) #gets rid of duplicate rows 

#now merge with all entries for metagenomics metadata (master)
##########################################################
meta_WGS$id <- paste(meta_WGS$Participant.ID,meta_WGS$visit_num,sep=".")
metaWGS2 <- meta_WGS[,c(2,3,9,11,12,13:18,32,33,44)]
meta_data_pre <- full_join(metaWGS2, meta_drugdata2,by="id") %>% 
		mutate(Participant.ID=coalesce(Participant.ID.x,Participant.ID.y), visit_num = coalesce(visit_num.x,visit_num.y)) %>%
			arrange(Participant.ID, visit_num) %>%
				select(-Participant.ID.x,-Participant.ID.y,-visit_num.x,-visit_num.y) %>%
				 	rename(age_dx=Age.at.diagnosis,hosp_2wks=X5..In.the.past.2.weeks..have.you.been.hospitalized.,bowelsurg=X6..Have.you.ever.had.bowel.surgery.,CRP=CRP..mg.L.,ESR=ESR..mm.hr.)

#carry forward within each participant so that non-measured data matches with stool data 
#this will be done first forward, then for those missing prior to baseline, will carry backwards
################################################################################
meta_data_pre <- meta_data_pre %>% mutate_at(13:20,as.numeric) #this gets rid of "" to NA

meta_data <- meta_data_pre %>% 
		group_by(Participant.ID) %>% 
			fill(c(diagnosis:bowelsurg,oral5asa_c:immunomod),.direction="down")%>%
				fill(c(diagnosis:bowelsurg,oral5asa_c:immunomod),.direction="up") 

#LOCF & NOCB, admittedly this is imperfect, but particpants were asked about meds at the start, middle, and end of the cohort, and most prescriptions for steroids (eg) are 6-12 weeks
#looked at the short med questionaire (any "immunosupp, eg steroids") but this seems unreliable because participants were on humira and entecort but said no immunosup 


#for non IBD patients with missing drug data, set to 0
meta_data[meta_data$diagnosis=="nonIBD",] %<>% mutate_at(vars(oral5asa_c:immunomod), ~replace(., is.na(.), 0))
#rename column 
meta_data <- meta_data %>% rename(SampleID = External.ID) 

cols <- c("oral5asa_lc","diagnosis","biologics","site_name","Antibiotics","bowelsurg",,"Participant.ID","steroids","pr5asa","bond5asa")
meta_data[cols] <- lapply(meta_data[cols], factor) #converts these to a factor variable)

meta_data$diagnosis <- relevel(meta_data$diagnosis, ref = "nonIBD") # Set nonIBD as reference


#this file is made by JLP for dysbiosis scoring (and also contains sex) 
dysbiosis_data <- read.table("/home/rsm34/5asa/hmp2_mgx_metadata.tsv", header=T,sep="\t")

dysbiosis_data <- dysbiosis_data[,c(1,4,5,9,11,19)]
dysbiosis_data <- dysbiosis_data %>% rename(SampleID = External_ID)
dysbiosis_data$SampleID <- as.character(dysbiosis_data$SampleID)

meta_data <- left_join(meta_data,dysbiosis_data,by="SampleID") 

meta_data$dysbiosis <- as.factor(meta_data$dysbiosis) 
meta_data$dysbiosis <- relevel(meta_data$dysbiosis, ref = "FALSE") # Set 'inactive' as reference


#print(meta_data,width=250)
write.table(meta_data, "meta_data.txt", sep="\t", quote=FALSE,row.names=FALSE)
#meta_data <- read.table("/home/rsm34/5asa/meta_data.txt",header=T,sep="\t")
