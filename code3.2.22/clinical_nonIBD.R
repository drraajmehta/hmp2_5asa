############################################################### CODE ###################################################

library(dplyr)
library(purrr)
library(table1) 
library(tidyr) 
library(nlme) 
library(ggplot2)
library(lme4)

#metadata prep 
##############

## step 1: find those who are "Ever users" of 5-ASA based on their metabolites
#read in metabolite prep file 
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt", header=T, sep="\t") %>% select(c("SampleID","binary154","X196.0609_2.81")) #saves time in the future. 
metab154 <- metab_select %>% filter(binary154 == 0) #178 samples with use 

#read in untouched metadata file & select columns of interest 
metadata <- read.csv("/home/rsm34/5asa/hmp2_metadata_load.csv",stringsAsFactors=FALSE)
metadata2 <- metadata %>% select(2,3,5,26:28,70:95,98,431,478,231,233,487) %>% rename(SampleID=External.ID) %>% filter(diagnosis=="nonIBD") 

userIDlist <- inner_join(metab154,metadata2,by="SampleID") %>% select(Participant.ID, binary154,SampleID,BMI,smoking.status) 
userIDlist2 <- userIDlist %>% distinct(Participant.ID) #yields 26 non-IBD people never using (not an enormous sample, does not count before mgx merge), not neccesarily with paired MGX data  

#merge with the metadata file to find ALL ever users
users_metadata <- inner_join(userIDlist2,metadata2,by="Participant.ID") #has all feature types here, therefore has some duplicate SampleIDs 
users_metadata2 <- users_metadata %>% arrange(SampleID) 

#now make sccai and hbi indices, ###******* might need to redefine these indices (<5, <3)  
users_metadata3 <- users_metadata2 %>%
     mutate(hbi_sccai = case_when(hbi < 4 ~ 0,
                                sccai <2 ~ 0,
                                 hbi > 3 ~ 1,
                          sccai > 1 ~ 1)) 

#now update smoking and BMI data, fill in missing with patient's data at timepoint one, otherwise set to median (43% missing smoking, 15% BMI) 
users_metadata3$smoking.status <- as.factor(users_metadata3$smoking.status)
users_metadata3 <- users_metadata3 %>% mutate_at(37:38,as.numeric) #this gets rid of "" to NA
users_metadata3$smoking.status[users_metadata3$smoking.status == 1] <- NA

users_metadata4 <- users_metadata3 %>% 
		group_by(Participant.ID) %>% 
			fill(c(BMI:smoking.status),.direction="down")%>%
				fill(c(BMI:smoking.status),.direction="up")
##numeric key 
#4=never smoked
#2=current smoker
#3=former smoker 

users_metadata4$smoking.status[users_metadata4$smoking.status == 2] <- "ever"
users_metadata4$smoking.status[users_metadata4$smoking.status == 3] <- "ever"
users_metadata4$smoking.status[users_metadata4$smoking.status == 4] <- "never"

users_metadata5 <- users_metadata4 %>% arrange(SampleID) %>% select(BMI,smoking.status,hbi_sccai,SampleID) %>%  distinct(SampleID,.keep_all=TRUE)

users_metadata5$smoking.status[is.na(users_metadata5$smoking.status)] <- "never" #43% missing smoking data
users_metadata5$BMI[is.na(users_metadata5$BMI)] <- median(users_metadata5$BMI,na.rm=TRUE) #15% missing BMI data

####now merge with my metadata file which has the other endpoints and metadata #### 
mymetadata <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t") 
mymeta2 <- inner_join(mymetadata, users_metadata5, by="SampleID") 
dim(mymeta2) #424, smaller because this is only among nonIBD patients 

#for purposes of pulling mgx files later
nonIBDmgxids <- mymeta2 %>% distinct(SampleID)
#write.table(nonIBDmgxids , file="nonIBDmgxids.txt", sep="\t", quote=FALSE,row.names=FALSE)


#some fine tuning

mymeta2$race[mymeta2$race == "American Indian or Alaska Native"] <- "nonwhite"
mymeta2$race[mymeta2$race == "Other"] <- "nonwhite"
mymeta2$race[mymeta2$race == "Black or African American"] <- "nonwhite"
mymeta2$race[is.na(mymeta2$race)] <- "nonwhite" #4% missing race data

mymeta2$sex[is.na(mymeta2$sex)] <- "Female" #4% missing sex data

cols <- c("smoking.status","diagnosis","bowelsurg","sex","race")
mymeta2[cols] <- lapply(mymeta2[cols], factor) #converts these to a factor variable)

mymeta2$smoking.status <- relevel(mymeta2$smoking.status, ref = "never") 
mymeta2$race <- relevel(mymeta2$race, ref = "White") 

mymeta2$bowelsurg[mymeta2$bowelsurg==""] <- NA

mymeta2<- mymeta2 %>%
     mutate(BMIcat = case_when(BMI < 25 ~ 0,
                                BMI >= 25 ~1) )

mymeta2$BMIcat <- as.factor(mymeta2$BMIcat)
mymeta2$Participant.ID.y <- as.factor(mymeta2$Participant.ID.y) 
 
################################################# now merge with MGX data ##########################################################

#extract DNA uniref90s (had to download an additional few unire90 files, see FAS)
######################################
### now extract  #####
#good to save workspace before leaving  
#done in UNIX -- Uniref DNA data in the shell b/c file is unwieldly
#$head -1 /home/rsm34/5asabb3/archive/DNA_uniref/DNAunirefproptable_clinc_nonIBD.txt > /home/rsm34/5asabb3/archive/DNA_uniref/uniref.header_clin_noIBD.txt 
#$grep -w -F -f /home/rsm34/5asabb3/intermediate_files/clinicauniref_match.txt /home/rsm34/5asabb3/archive/DNA_uniref/DNAunirefproptable_clinc_nonIBD.txt > clindnaselect_noIBD.txt 

#read in files and combine 
unirefselect <- read.table("/home/rsm34/5asabb3/intermediate_files/clindnaselect_noIBD.txt",header=F,sep="\t") 
unirefheader <- read.table("/home/rsm34/5asabb3/archive/DNA_uniref/uniref.header_clin_noIBD.txt",header=F,sep="\t") 
unirefselectcomb <- rbind(unirefheader,unirefselect)
colnames(unirefselectcomb) <- unirefselectcomb[1,]
dnauniselect <- unirefselectcomb[2:nrow(unirefselectcomb),] 

#transpose and sort again 
dim(dnauniselect) 
dnauniselect_rotated = setNames(data.frame(t(dnauniselect[,-1])), dnauniselect[,1])
dnauniselect_rotated<- cbind (SampleID=row.names(dnauniselect_rotated),dnauniselect_rotated)
row.names(dnauniselect_rotated) <- NULL
dim(dnauniselect_rotated)

dnauniselect_rotated_bn <- dnauniselect_rotated[,2:13]
dnauniselect_rotated_bn[dnauniselect_rotated_bn >0] <-1
dnauniselect_rotated_bn2 <- cbind(dnauniselect_rotated$SampleID,dnauniselect_rotated_bn) %>% rename(SampleID = "dnauniselect_rotated$SampleID")


mymeta3 <- inner_join(mymeta2, dnauniselect_rotated_bn2, by="SampleID") 
dim(mymeta3) 

mymeta3[,35:46] <- mymeta3 %>% select(35:46) %>% mutate_if(is.character,as.factor)


############################################################################################### testing gene DATA #######################################################################
#first create datasets 
df_steroids <- mymeta3 %>% drop_na(steroids,consent_age) 
df_steroids$steroids <- as.factor(df_steroids$steroids)

#### get proportions for users #### #n=384
prop.table(table(df_steroids$UniRef90_C7H1G6)) 
#0.3298429 0.6701571
prop.table(table(df_steroids$UniRef90_R6CZ24))
#0.7146597 0.2853403
prop.table(table(df_steroids$UniRef90_R5CY66))
#0.5628272 0.4371728
prop.table(table(df_steroids$UniRef90_T5S060))
#0.8769634 0.1230366
