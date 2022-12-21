################################################################
#					#
#	TABLE 1: 5asa users vs other (n=95)		#
#					#
################################################################

library(dplyr)
library(tidyr)
library(tableone) 
library(ggplot2)

#metadata prep  
##############

## step 1: find those who are "Ever users" of 5-ASA based on their metabolites
#read in metabolite prep file 
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt", header=T, sep="\t") %>% select(c("SampleID","binary154","X196.0609_2.81")) #saves time in the future. 
metab154 <- metab_select %>% filter(binary154 == 1) #178 samples with use, doesn't clarify IBD non IBD status (although should all be non IBD) 

#read in untouched metadata file & select columns of interest 
metadata <- read.csv("/home/rsm34/5asa/hmp2_metadata_load.csv",stringsAsFactors=FALSE)
metadata2 <- metadata %>% select(2,3,5,26:28,70:95,98,431,478,231,233,487) %>% rename(SampleID=External.ID) %>% filter(diagnosis=="CD"|diagnosis=="UC") 

##before we go further, let's get numbers... find how many IBD pts gave paired metagenomic and metabolomic samples and how many samples they gave
## this belongs in the text  

mgxnum <- metadata %>% filter(data_type=="metagenomics") %>% rename(SampleID=External.ID) %>% filter(diagnosis=="CD"|diagnosis=="UC")
mtxnum <- metadata %>% filter(data_type=="metatranscriptomics") %>% rename(SampleID=External.ID) %>% filter(diagnosis=="CD"|diagnosis=="UC")
mbxnum <- metadata %>% filter(data_type=="metabolomics") %>% rename(SampleID=External.ID) %>% filter(diagnosis=="CD"|diagnosis=="UC")
mbx_mgxnum <- inner_join(mbxnum,mgxnum,by="SampleID") #283
mbx_mgxnum %>% distinct(Participant.ID.x)#79
mtx_mbx_mgxnum <- inner_join(mtxnum,mbx_mgxnum,by="SampleID") #220

### okay now continue 

userIDlist <- inner_join(metab154,metadata2,by="SampleID") %>% select(Participant.ID, binary154,SampleID,BMI,smoking.status,data_type) 
userIDlist2 <- userIDlist %>% distinct(Participant.ID) #yields 46 IBD people ever users (but might have extra)

#THIS GENERATES NUMBERS OF USERS VS NONUSERS, different than above because one doesn't have paired MGX  

mbxPID <- mbx_mgxnum %>% distinct(Participant.ID.x)
mbxPID2 <- mbxPID %>% mutate(evernever154 = if_else(Participant.ID.x %in% userIDlist2$Participant.ID,"yes","no")) %>% rename(Participant.ID=Participant.ID.x) 
table(mbxPID2$evernever154) ######## THIS GIVES USERS VERSUS NON USERS####

#merge with the metadata file 
users_never_metadata <- inner_join(metadata2,mbxPID2,by="Participant.ID") #has all feature types here, therefore has some duplicate SampleIDs 
table1_metadata2 <- users_never_metadata %>% arrange(SampleID)  

#now update smoking and BMI data, fill in missing with patient's data at timepoint one, otherwise set to median (43% missing smoking, 15% BMI) 
table1_metadata2$smoking.status <- as.factor(table1_metadata2$smoking.status)
table1_metadata3 <- table1_metadata2 %>% mutate_at(37:38,as.numeric) #this gets rid of "" to NA
table1_metadata3$smoking.status[table1_metadata3$smoking.status == 1] <- NA

table1_metadata4 <- table1_metadata3 %>% 
		group_by(Participant.ID) %>% 
			fill(c(BMI:smoking.status),.direction="down")%>%
				fill(c(BMI:smoking.status),.direction="up")
##numeric key 
#4=never smoked
#2=current smoker
#3=former smoker 

table1_metadata4$smoking.status[table1_metadata4$smoking.status == 2] <- "ever"
table1_metadata4$smoking.status[table1_metadata4$smoking.status == 3] <- "ever"
table1_metadata4$smoking.status[table1_metadata4$smoking.status == 4] <- "never"

table1_metadata5 <- table1_metadata4 %>% arrange(SampleID) %>% select(BMI,smoking.status,SampleID,evernever154) %>%  distinct(SampleID,.keep_all=TRUE)

table1_metadata5$smoking.status[is.na(table1_metadata5$smoking.status)] <- "never" #to missing smoking data
table1_metadata5$BMI[is.na(table1_metadata5$BMI)] <- median(table1_metadata5$BMI,na.rm=TRUE) #median to missing BMI data

####now merge with my metadata file which has the other endpoints and metadata #### 
mymetadata <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t") 
mymeta2 <- inner_join(mymetadata, table1_metadata5, by="SampleID") 
dim(mymeta2) #1168, smaller because this is only among IBD patients who have metabolomics data  

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))] 
}


cols <- c("oral5asa_lc","diagnosis","biologics","site_name","Antibiotics","bowelsurg","dysbiosis","Participant.ID.x","steroids","pr5asa","bond5asa","Chemotherapy","hosp_2wks","smoking.status","evernever154")
mymeta2[cols] <- lapply(mymeta2[cols], factor) #converts these to a factor variable)
str(mymeta2)

mymeta2$dysbiosis <- relevel(mymeta2$dysbiosis, ref = "FALSE") # Set 'inactive' as reference
mymeta2 <- mymeta2 %>% mutate(IBDusers = if_else(oral5asa_lc == 1,2,if_else(oral5asa_lc==0 & as.numeric(diagnosis) == 3,1, if_else(oral5asa_lc==0 & as.numeric(diagnosis) == 2,1,0)))) 
mymeta2$IBDusers[mymeta2$IBDusers=="2"] <- "User"
mymeta2$IBDusers[mymeta2$IBDusers=="1"] <- "IBDNon-User"
mymeta2$IBDusers[mymeta2$IBDusers=="0"] <- "nonIBD"
mymeta2$IBDusers <- as.factor(mymeta2$IBDusers)
mymeta2$IBDusers <- relevel(mymeta2$IBDusers, ref = "nonIBD") # Set nonIBD as reference

#subset to IBD users only 
mymeta2_IBD <- mymeta2 %>% filter(diagnosis=="CD"|diagnosis=="UC") 

mymeta2_IBD$Antibiotics <- (as.numeric(mymeta2_IBD$Antibiotics)-1) 
mymeta2_IBD$bond5asa<- (as.numeric(mymeta2_IBD$bond5asa)-1) 
mymeta2_IBD$biologics <- (as.numeric(mymeta2_IBD$biologics)-1) 
mymeta2_IBD$steroids <- (as.numeric(mymeta2_IBD$steroids)-1) 
mymeta2_IBD$immunomod <- (as.numeric(mymeta2_IBD$immunomod)) 
mymeta2_IBD$dysbiosis <- (as.numeric(mymeta2_IBD$dysbiosis)-1) 
mymeta2_IBD$Chemotherapy <- (as.numeric(mymeta2_IBD$Chemotherapy)-1) 
mymeta2_IBD$hosp_2wks <- (as.numeric(mymeta2_IBD$hosp_2wks)-1) 
mymeta2_IBD$evernever154 <- (as.numeric(mymeta2_IBD$evernever154)-1) 


#Make table 1 with bonded 5-asa data 
###################################
table1data <- mymeta2_IBD

table1data2 <- table1data %>% 
	group_by(Participant.ID.x) %>% 
		summarize(mesaluse=max(evernever154),bio=max(biologics,na.rm=TRUE),immuno=max(immunomod,na.rm=TRUE), ster=max(steroids,na.rm=TRUE),age_dx=max(age_dx,na.rm=TRUE),age_consent=max(consent_age,na.rm=TRUE),dysbiosis=max(dysbiosis,na.rm=TRUE),abx=max(Antibiotics),
		dx=getmode(diagnosis), bond=max(bond5asa),chemo=max(Chemotherapy),hospital=max(as.numeric(hosp_2wks)),surgery=getmode(bowelsurg),calpro=max(fecalcal,na.rm=TRUE),sex=getmode(sex),race=getmode(race))

#replaces inf with NA 
table1data_clean <- do.call(data.frame,lapply(table1data2, function(x) replace(x, is.infinite(x),NA)))

table1data_clean <- table1data_clean %>% 
  mutate(peds = ifelse(age_consent < 16, 1, 0)) 

table1data_clean <- table1data_clean %>% 
  mutate(peds = ifelse(age_consent < 16, 1, 0)) 


myvars <- c("age_consent", "age_dx","dx","abx", "hospital", "immuno","ster","bio","dysbiosis","sex","race","peds","surgery","bond")
catvars<-c("abx", "hospital", "immuno","ster","bio","dysbiosis","sex","race","peds","dx","surgery","bond")

tab1 <- CreateTableOne(vars=myvars,data=table1data_clean,factorVars=catvars,strata="mesaluse",test=TRUE,includeNA=TRUE)
tab1mat <- print(tab1,formatOptions=list(big.mark=","),quote=FALSE, noSpaces=TRUE,printToggle=FALSE)
#write.csv(tab1mat,file="/home/rsm34/5asabb3/reviews/table1_bonded5asa.csv") #to be formatted in word 

##### find ids of bonded users to filter later #####

mymeta2_IBD %>% select(SampleID, evernever154,bond5asa,Participant.ID.x) %>% filter(bond5asa==0)

#### find concordance between self-report use and metabolite defined use #####
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt", header=T, sep="\t")
mymetadata <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t") 

mymetadata$SampleID <- as.character(mymetadata$SampleID)
cols <- c("oral5asa_lc","diagnosis","biologics","site_name","Antibiotics","bowelsurg","Participant.ID","steroids","pr5asa","bond5asa","sex","dysbiosis")
mymetadata[cols] <- lapply(mymetadata[cols], factor) #converts these to a factor variable)

raaj <- mymetadata %>% mutate(IBDusers2 = if_else(any5asa == 1,2,if_else(any5asa==0 & as.numeric(diagnosis) == 3,1, if_else(any5asa==0 & as.numeric(diagnosis) == 1,1,0)))) 
raaj$IBDusers2[raaj$IBDusers2=="2"] <- "User"
raaj$IBDusers2[raaj$IBDusers2=="1"] <- "IBDNon-User"
raaj$IBDusers2[raaj$IBDusers2=="0"] <- "nonIBD"
raaj$IBDusers2 <- as.factor(raaj$IBDusers2)
raaj$IBDusers2 <- relevel(raaj$IBDusers2, ref = "nonIBD") # Set nonIBD as reference

concord <- inner_join(metab_select,raaj, by="SampleID") %>% drop_na(any5asa) #removes missing drug data 

table(concord$IBDusers2,concord$binary154)                                    

#              0   1
#  nonIBD      105   0
#  IBDNon-User 115  45
#  User         23  57

#accuracy is TN + TP/ TN + TP + FN + FP
#(105+115+57)/(105+115+57+45+23) = .80289


pq <- ggplot(concord, aes(x=IBDusers2, y=X154.0502_3.83, color=binary154)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))

pq2 <- pq + scale_y_log10() + theme_bw() +labs(y= "Fecal 5-ASA level", x = "")+
	scale_x_discrete(labels=c("nonIBD" = "Non IBD", "IBDNon-User" = "IBD Non-user",
                              "User" = "5-ASA User")) + theme(text = element_text(size = 20))  + theme(legend.position="none")  

library(ggExtra)
pq3<-ggMarginal(pq2, groupColour = TRUE, margins="y",groupFill = TRUE)

ggsave(filename='/home/rsm34/5asabb3/reviews/suppfig1_rev.png', plot=pq3, width = 7, height = 6, dpi = 600) 


