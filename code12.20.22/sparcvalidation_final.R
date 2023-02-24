############################################# THESE ARE BASH COMMANDS TO BE DONE IN UNIX #############################################################################################
#steps to consolidate all relevant files 

#download file /n/holylfs05/LABS/nguyen_lab/Lab/data/molecular-data-production/clinical-data-by-project/Mehta_Cohort_2021-10-19.xlsx #n=300

#generate a list of IDs in excel
get rid of "_R1..." and concatenate "_genefamilies_relab.tsv" 

#then make sure unix compatible 
dos2unix /n/holyscratch01/nguyen_lab/raaj/sparcfilelist.txt

#then copy and move with script below  
cd to five paths: 

	/n/holyscratch01/nguyen_lab/raaj/sparc/sparc_cduc_diversigen_2021_subset_RAAJ/OUTPUT_subset/humann/relab/genes #n=103 (103)

	/n/holyscratch01/nguyen_lab/raaj/sparc/sparc_cduc_diversigen/humann/relab/genes #n=302  (105) 

	/n/holyscratch01/nguyen_lab/raaj/sparc/sparc_cduc_diversigen_2021/humann/relab/genes #n=221 (23) 

	/n/holyscratch01/nguyen_lab/raaj/sparc/original_sparc_runs/humann/relab/genes #n=402 (38) 

	/n/holyscratch01/nguyen_lab/raaj/sparc_cduc/humann/relab/genes #n=206 

file in $(cat /n/holyscratch01/nguyen_lab/raaj/sparcfilelist.txt); do cp "$file" /n/holyscratch01/nguyen_lab/raaj/temp/; done

#then collapse all into one file (load hutlab modules)  
humann_join_tables --input /n/holyscratch01/nguyen_lab/raaj/temp --output sparc300.txt

#shrink file 
# GREP not: 
$ grep -v "g_" sparc300.txt > UnirefPrecursorA.txt #keeps only the unirefs 
$ grep -v "|" UnirefPrecursorA.txt > UnirefPrecursorB.txt #gets rid of unclassified
$ sed 's/^#.//' UnirefPrecursorB.txt > UnirefPrecursorC.txt #gets rid of hashtag on first line 
$ wc -l UnirefPrecursorC.txt #2907276 (confirmed that this sums to 1) 

#now grep the relevant genes 
grep -w -F -f /n/holyscratch01/nguyen_lab/raaj/clinicauniref_match.txt /n/holyscratch01/nguyen_lab/raaj/UnirefPrecursorC.txt > clindnaselect.txt 

#now grep header
head -1 /n/holyscratch01/nguyen_lab/raaj/UnirefPrecursorC.txt > uniref.header_clin.txt

################################################################################################################################################################################

#now bring into R and process a bit more 

library(dplyr) 
library(readxl)

#read in files and combine (have more than just users)
unirefselect <- read.table("/n/holyscratch01/nguyen_lab/raaj/clindnaselect.txt",header=F,sep="\t") 
unirefheader <- read.table("/n/holyscratch01/nguyen_lab/raaj/uniref.header_clin.txt",header=F,sep="\t") 
unirefselectcomb <- rbind(unirefheader,unirefselect)
colnames(unirefselectcomb) <- unirefselectcomb[1,]
dnauniselect <- unirefselectcomb[2:nrow(unirefselectcomb),] 

#transpose and sort again 
dim(dnauniselect) 
dnauniselect_rotated = setNames(data.frame(t(dnauniselect[,-1])), dnauniselect[,1])
dnauniselect_rotated<- cbind (SampleID=row.names(dnauniselect_rotated),dnauniselect_rotated)
row.names(dnauniselect_rotated) <- NULL
dim(dnauniselect_rotated)

#now make binary (presence /absence) 
dnauniselect_rotated_bn <- dnauniselect_rotated[,2:13]
dnauniselect_rotated_bn[dnauniselect_rotated_bn >0] <-1
dnauniselect_rotated_bn2 <- cbind(dnauniselect_rotated$SampleID,dnauniselect_rotated_bn) %>% rename(SampleID = "dnauniselect_rotated$SampleID")



#merge with metadata 
ccfmeta <- read_excel("/n/holylfs05/LABS/nguyen_lab/Lab/data/molecular-data-production/clinical-data-by-project/Mehta_Cohort_2021-10-19.xlsx") %>% rename(patient_id = DEIDENTIFIED_MASTER_PATIENT_ID)
ccfmeta$SampleID <- gsub("R2.fastq.gz","Abundance-RPKs",ccfmeta$RAW.DATA.FILE.NAME_1)
ccfmeta$SampleID <- gsub("R1.fastq.gz","Abundance-RPKs",ccfmeta$SampleID)
mymeta3 <- inner_join(ccfmeta, dnauniselect_rotated_bn2, by="SampleID") 
dim(mymeta3) #692

library(tidyr)
library(stringr)
mymeta3$STEROIDS = mymeta3$STEROIDS %>% replace_na('none')
mymeta3[,c(1,4,7:9,13:24)] <- mymeta3 %>% select(1,4,7:9,13:24) %>% mutate_if(is.character,as.factor)
mymeta3$collyr <- substr(mymeta3$DATE_OF_CONSENT,1,4)
mymeta3$collyr<-as.numeric(mymeta3$collyr)
mymeta3$BIRTH_YEAR <- as.numeric(mymeta3$BIRTH_YEAR)
mymeta3$ageyr <- mymeta3$collyr - mymeta3$BIRTH_YEAR

## some filtration (eg get rid of current steroid users for this to be prospective, get rid of those who withdrew early), 
mymeta3 <- mymeta3 %>% filter(is.na(DATE_OF_CONSENT_WITHDRAWN)) %>% filter(!STEROIDS=="Current")  

#### now filter those non-indepedent samples (samples provided fewer than 30 days apart, b/c most steroid doses are 30 days long) 
mymeta3<-mymeta3 %>% arrange(patient_id, SAMPLE_COLLECTED_DATE) %>% group_by(patient_id) %>%
  mutate(diffDate = difftime(SAMPLE_COLLECTED_DATE, units="days",lag(SAMPLE_COLLECTED_DATE,1))) %>% ungroup()%>% filter(diffDate >30|is.na(diffDate))

mymeta3$unique <-paste0(mymeta3$patient_id,sep=".",mymeta3$RAW.DATA.FILE.NAME_1)


#censor those samples who provided stool after having the event (this was determined through a manual curation) 

mymeta3 <- mymeta3 %>% filter(!unique=="15115738.FB06103874_MQ2867_WGS_R1.fastq.gz") %>% filter(!unique=="15113540.FB05731806_MQ2867_WGS_R2.fastq.gz") %>% filter(!unique=="15113540.FB06104840_MQ2867_WGS_R2.fastq.gz") 

#15115738.FB06103874_MQ2867_WGS_R1.fastq.gz
#15113540.FB05731806_MQ2867_WGS_R2.fastq.gz
#15113540.FB05731806_MQ2867_WGS_R2.fastq.gz

############################################################################################### baseline characteristics ########################################################################

mean(mymeta3$ageyr)
#45.5
sd(mymeta3$ageyr)
#15

table(mymeta3$STEROIDS)
#  After    none
#     60     190

mymeta_uni <- mymeta3 %>% distinct(patient_id,.keep_all=TRUE)

table(mymeta_uni$DIAGNOSIS)
#Ulcerative Colitis    Crohn's Disease   IBD Unclassified
#               140                 67                  1

table(mymeta_uni$GENDER)
#Female   Male
#   140    113

mymeta3<- mymeta3 %>% mutate(steroidsany = if_else(STEROIDS=="After"|STEROIDS=="Current",1,0))
mymeta3$steroidsany <- as.factor(mymeta3$steroidsany) 
mymeta3$DIAGNOSIS <- relevel(mymeta3$DIAGNOSIS, ref="Ulcerative Colitis")

#check proportion present/absent 
DF<-mymeta3 %>% select(UniRef90_C7H1G6,UniRef90_R5CY66 ,UniRef90_R6CZ24 ,UniRef90_T5S060)
DF <- sapply(DF,as.numeric)
DF <- DF-1
colSums(DF)/nrow(DF)
#UniRef90_C7H1G6 UniRef90_R5CY66 UniRef90_R6CZ24 UniRef90_T5S060
#          0.440           0.204           0.072           0.188


############################################################################################### now performing analysis/testing  #######################################################################
#first create datasets and make the score variable 
df_steroids3 <- mymeta3 %>% drop_na(steroidsany,ageyr) 
dim(df_steroids3) #--> 250 samples 
df_steroids3 %>% distinct(patient_id) #from 208 individuals 
median(df_steroids3$diffDate,na.rm=T)#among those with repeated samples, median interval sampling is 133 days (4.4mos)

##########################
### OUTCOME: STEROIDS ####
##########################

### check the score like in the submission #####
##########################################################################

df_steroids_score <- df_steroids3 
df_steroids_score <- df_steroids_score%>% 
  mutate_at(vars(UniRef90_C7H1G6 , UniRef90_R5CY66 , UniRef90_R6CZ24 , UniRef90_T5S060), ~as.numeric(as.character(.))) 

df_steroids_score$score <- apply(df_steroids_score[,c(18,20,21,23)], 1, sum) 

#look at 2 categories (0-2,3+)  **** without adjustment for repeated measures ***** 
df_steroids_score<- df_steroids_score  %>% mutate(score3 = ifelse(score > 2, 1, 0)) 
score_lm3 <- glm(steroidsany~ ageyr + DIAGNOSIS + as.factor(score3)+GENDER, data=df_steroids_score, family=binomial) #gives by level 
summary(score_lm3)

exp(0.9801)
exp(0.9801-1.96* 0.5035) 
exp(0.9801+1.96* 0.5035)
#2.66 (0.99-7.15)

##### suspect our power is somewhat low ####
#ageyr                      -0.01651    0.01043  -1.582   0.1136
#DIAGNOSISCrohn's Disease    0.39641    0.31785   1.247   0.2123


######### gee (included in final manuscript)
#https://data.library.virginia.edu/getting-started-with-generalized-estimating-equations/ good for small number of small clusters 
library(gee)
try_gee <- gee(steroidsany~  as.factor(score3)+ageyr +GENDER,
               data = df_steroids_score, 
               id = as.factor(patient_id), 
               family = binomial,
               corstr = "independence")
summary(try_gee)
exp(1.0185)
exp(1.0185-1.96*0.5036)
exp(1.0185+1.96*0.5036)

#2.77 (1.03-7.43)
####

######### sensitivity analysis #1: I exclude only to baseline sample  ################
test<-df_steroids_score %>% group_by(patient_id) %>% filter(n() > 1) %>% select(patient_id,score3,SAMPLE_COLLECTED_DATE,steroidsany,unique)

data.frame(test)

test2 <- test[c(2,4,5,7,9,11,12,14,16,18,19,20,23,25,27,29,31,33,35,36,38,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,68,70,72,74,76,77),]

`%notin%` <- Negate(`%in%`)

df_steroids_score2 <-subset(df_steroids_score, unique %notin% test2$unique)


score_lm3_v2 <- glm(steroidsany~ ageyr + DIAGNOSIS + as.factor(score3)+GENDER, data=df_steroids_score2, family=binomial) #gives by level 
summary(score_lm3_v2)

exp(0.80-1.96*0.5616)
exp(0.80)
exp(0.80+1.96*0.5616)

#2.22 (0.74-6.69) ### still significant in meta-analysis, definitely underpowered

######### sensitivity analysis #2: I try mixed effects models (not as good given small number of small clusters)  ################
#https://data.library.virginia.edu/getting-started-with-generalized-estimating-equations/

#pt 6: glmm (null) 
#try mixed effects models for a few repeated measures --> null 
library(GLMMadaptive)
summary(mixed_model(fixed=steroidsany~ as.factor(score3), random= ~1|patient_id, data=df_steroids_score, family=binomial(link="logit"))) #gives binary 

#p=0.28


############################################################### META-ANALYSIS (done in R studio) ########################################################

library(metafor)

#### create dataframe to store these values for STEROIDS ### 

dat = NULL
dat = as.data.frame(matrix(-9,2,3))
colnames(dat) = c("study","yi","sei")
dat[1,1]<-"HMP2"
dat[2,1]<-"SPARC IBD"
dat[1,2]<-1.353
dat[2,2]<-1.0185
dat[1,3]<-0.6809
dat[2,3]<-0.5036
 
### fit random-effects model
res <- rma(yi, sei=sei, data=dat)

### forest plot with extra annotations
forest(res, atransf=exp, at=log(c(.05, .25, 1, 4,16)), xlim=c(-6,10),
       cex=1, header="Study", mlab="",slab=c("IBDMDB","SPARC IBD"))

k <- nrow(dat)
psize <- weights(res)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))
points(dat$yi, k:1, pch=15, cex=psize*1.15,col="white")
points(dat$yi, k:1, pch=22, cex=psize, col="black",bg="red")
 
### add the summary polygon
addpoly(res, row=-1, mlab="", efac=4, col="pink", border="black")

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-6, 3.5, pos=4, cex=0.9, bquote(paste("RE Model (Q = ",
     .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
     ", p = ", .(formatC(res$QEp, digits=3, format="f")), "; ", I^2, " = ",
     .(formatC(res$I2, digits=1, format="f")), "%)")))

#label overall 
text(-6, -1, pos=4, cex=1, bquote(paste("Overall")))


### lablel x axis (doesn't work) 
text(1, -2.5, pos=4,cex=0.75,bquote("Risk Ratio (log scale)"))


################################################################## check median follow up time for text using newly supplied file ###########################################

followuptime <- read.table("/n/home02/rmehta/plexus/validation/Mehta_Cohort_2021-10-19-updated_072022.txt",header=T,sep="\t") %>% rename(patient_id = DEIDENTIFIED_MASTER_PATIENT_ID) %>% select(patient_id,CORTICOSTEROIDS_ECRF,CORTICOSTEROIDS_EMR,RAW.DATA.FILE.NAME_1)
followuptime$patient_id <- as.character(followuptime$patient_id)

#generate unique ids to merge on 

ccfmeta$unique <-paste0(ccfmeta$patient_id,sep=".",ccfmeta$RAW.DATA.FILE.NAME_1)
followuptime$unique <-paste0(followuptime$patient_id,sep=".",followuptime$RAW.DATA.FILE.NAME_1)

ccfmeta2 <- inner_join(ccfmeta,followuptime, by="unique")


#first filter those rows that have steroid users after time of dx 

afterusers <- ccfmeta2 %>% filter(grepl('After', CORTICOSTEROIDS_EMR)) #not all after users have time data available 

library(tidyr)
median(extract_numeric(afterusers$CORTICOSTEROIDS_EMR ))
IQR(extract_numeric(afterusers$CORTICOSTEROIDS_EMR ))

####duplicates (to identify those to subset above)

duplicates <-mymeta3 %>% group_by(patient_id) %>% filter(n() > 1) %>% select(patient_id, SAMPLE_COLLECTED_DATE,score3,steroidsany,unique)


################################################################## Ordination plot ###########################################

followuptime <- read.table("/n/home02/rmehta/plexus/validation/Mehta_Cohort_2021-10-19-updated_072022.txt",header=T,sep="\t") %>% rename(patient_id = DEIDENTIFIED_MASTER_PATIENT_ID) %>% select(patient_id,CORTICOSTEROIDS_ECRF,CORTICOSTEROIDS_EMR,RAW.DATA.FILE.NAME_1)
followuptime$patient_id <- as.character(followuptime$patient_id)

#generate unique ids to merge on 

ccfmeta$unique <-paste0(ccfmeta$patient_id,sep=".",ccfmeta$RAW.DATA.FILE.NAME_1)
followuptime$unique <-paste0(followuptime$patient_id,sep=".",followuptime$RAW.DATA.FILE.NAME_1)

ccfmeta2 <- inner_join(ccfmeta,followuptime, by="unique")


#first filter those rows that have steroid users after time of dx 

afterusers <- ccfmeta2 %>% filter(grepl('After', CORTICOSTEROIDS_EMR)) #not all after users have time data available 

library(tidyr)
median(extract_numeric(afterusers$CORTICOSTEROIDS_EMR )) 
IQR(extract_numeric(afterusers$CORTICOSTEROIDS_EMR ))

####duplicates (to identify those to subset above)

duplicates <-mymeta3 %>% group_by(patient_id) %>% filter(n() > 1) %>% select(patient_id, SAMPLE_COLLECTED_DATE,score3,steroidsany,unique)






