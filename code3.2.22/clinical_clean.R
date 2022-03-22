######################################################################################################
# R programs for creating "clinical" plots					#
# ============================================					#
# Purposes:  odds ratios,   						#
#              distribution of P/A genes, host genetics 				#
######################################################################################################


############################################################### CODE ###################################################

library(dplyr)
library(purrr)
library(table1) 
library(tidyr) 
library(nlme) 
library(ggplot2)

#metadata prep -- extract BMI, smoking.status, etc for MV models
################################################################

## step 1: find those who are "Ever users" of 5-ASA based on their metabolites
#read in metabolite prep file (more condensed version with relevant metabs) #saves time in the future
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt", header=T, sep="\t") %>% select(c("SampleID","binary154","X196.0609_2.81")) 
metab154 <- metab_select %>% filter(binary154 == 1) #178 samples with use #we know there are only 80 IBD Pts who gave metabolites

#read in untouched metadata file & select columns of interest #NB i used col numbers here(try not to) bc var names are extremely long eg
####### "Soft.drinks..tea.or.coffee.with.sugar..corn.syrup..maple.syrup..cane.sugar..etc."

metadata <- read.csv("/home/rsm34/5asa/hmp2_metadata_load.csv",stringsAsFactors=FALSE)
metadata2 <- metadata %>% select(2,3,5,26:28,70:95,98,431,478,231,233,487) %>% rename(SampleID=External.ID) %>% filter(diagnosis=="CD"|diagnosis=="UC") 

#merge datasets and get list of unique IDs for users 
userIDlist <- inner_join(metab154,metadata2,by="SampleID") %>% select(Participant.ID, binary154,SampleID,BMI,smoking.status) 
userIDlist2 <- userIDlist %>% distinct(Participant.ID) #yields 46 IBD people ever using (not an enormous sample, not neccesarily with paired MGX data)  

#merge with the metadata file to find data across all time points for ever users
users_metadata <- inner_join(userIDlist2,metadata2,by="Participant.ID") #has all feature types here, therefore has some duplicate SampleIDs, but need all rows to fill missing below 
users_metadata2 <- users_metadata %>% arrange(SampleID) 

#now update smoking and BMI data, fill in missing with patient's data at timepoint one, otherwise set to median (43% missing smoking, 15% BMI) 
users_metadata2$smoking.status <- as.factor(users_metadata2$smoking.status)
users_metadata2 <- users_metadata2 %>% mutate_at(37:38,as.numeric) #this gets rid of "" to NA
users_metadata2$smoking.status[users_metadata2$smoking.status == 1] <- NA

#yes i skipped from 2->4
users_metadata4 <- users_metadata2 %>% 
		group_by(Participant.ID) %>% 
			fill(c(BMI:smoking.status),.direction="down")%>%
				fill(c(BMI:smoking.status),.direction="up")
#recode smoking 
##numeric key 
#4=never smoked
#2=current smoker
#3=former smoker 

users_metadata4$smoking.status[users_metadata4$smoking.status == 2] <- "ever"
users_metadata4$smoking.status[users_metadata4$smoking.status == 3] <- "ever"
users_metadata4$smoking.status[users_metadata4$smoking.status == 4] <- "never"

#get file with BMI and smoking status by each time point (fairly static) 
users_metadata5 <- users_metadata4 %>% arrange(SampleID) %>% select(BMI,smoking.status,SampleID) %>%  distinct(SampleID,.keep_all=TRUE)

users_metadata5$smoking.status[is.na(users_metadata5$smoking.status)] <- "never" #42% missing smoking data
users_metadata5$BMI[is.na(users_metadata5$BMI)] <- median(users_metadata5$BMI,na.rm=TRUE) #16% missing BMI data

####now merge with my metadata file which has the other endpoints and metadata #### 
mymetadata <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t") 
mymeta2 <- inner_join(mymetadata, users_metadata5, by="SampleID") 
dim(mymeta2) #692 x 32 , smaller because this is only among IBD patients who have ever used mesalamine containing products 

#some fine tuning of existing variables 
mymeta2$race[mymeta2$race == "American Indian or Alaska Native"] <- "nonwhite"
mymeta2$race[mymeta2$race == "Other"] <- "nonwhite"
mymeta2$race[mymeta2$race == "Black or African American"] <- "nonwhite"
mymeta2$race[is.na(mymeta2$race)] <- "nonwhite" #4% missing race data

mymeta2$sex[is.na(mymeta2$sex)] <- "Female" #4% missing sex data

cols <- c("smoking.status","diagnosis","bowelsurg","sex","race")
mymeta2[cols] <- lapply(mymeta2[cols], factor) #converts these to a factor variable)

mymeta2$smoking.status <- relevel(mymeta2$smoking.status, ref = "never") 
mymeta2$diagnosis <- relevel(mymeta2$diagnosis, ref = "UC") # Set UC as reference
mymeta2$race <- relevel(mymeta2$race, ref = "White") 

mymeta2$bowelsurg[mymeta2$bowelsurg==""] <- NA

mymeta2<- mymeta2 %>%
     mutate(BMIcat = case_when(BMI < 25 ~ 0,
                                BMI >= 25 ~1) )

mymeta2$BMIcat <- as.factor(mymeta2$BMIcat)
mymeta2$Participant.ID.y <- as.factor(mymeta2$Participant.ID.y) 
 
#####overview of endpoint######
table(mymeta2$steroids) #this is smaller than total b/c of missing steroid data (self-reported) 
#  0   1
#486 123

################################################# now prepare exposure data: extract and merge MGX data ##########################################################

#extract DNA uniref90s (had to download bb3 uniref90 files, see FAS), extract the 12 candidates from the overall file of users 
######################################
#good to save workspace before leaving  
#done in UNIX -- Uniref DNA data in the shell b/c file is unwieldly
#$head -1 /home/rsm34/5asabb3/archive/DNA_uniref/DNAunirefproptable_clinc.txt > /home/rsm34/5asabb3/archive/DNA_uniref/uniref.header_clin.txt 
#$grep -w -F -f /home/rsm34/5asabb3/submission/keyoutput/clinicauniref_match.txt /home/rsm34/5asabb3/archive/DNA_uniref/DNAunirefproptable_clinc.txt > /home/rsm34/5asabb3/submission/other/clindnaselect.txt 

#read in files and combine (have more than just users)
unirefselect <- read.table("/home/rsm34/5asabb3/submission/other/clindnaselect.txt",header=F,sep="\t") 
unirefheader <- read.table("/home/rsm34/5asabb3/archive/DNA_uniref/uniref.header_clin.txt",header=F,sep="\t") 
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
dnauniselect_rotated_bn <- dnauniselect_rotated[,2:14]
dnauniselect_rotated_bn[dnauniselect_rotated_bn >0] <-1
dnauniselect_rotated_bn2 <- cbind(dnauniselect_rotated$SampleID,dnauniselect_rotated_bn) %>% rename(SampleID = "dnauniselect_rotated$SampleID")

#merge with metadata 
mymeta3 <- inner_join(mymeta2, dnauniselect_rotated_bn2, by="SampleID") 
dim(mymeta3) #692

mymeta3[,34:46] <- mymeta3 %>% select(34:46) %>% mutate_if(is.character,as.factor)


############################################################################################### now performing analysis/testing  #######################################################################
#first create datasets 
df_steroids <- mymeta3 %>% drop_na(steroids,consent_age) 
df_steroids$steroids <- as.factor(df_steroids$steroids)
dim(df_steroids) #--> 609 samples 

##########################
### OUTCOME: STEROIDS ####
##########################

#univariate analysis (sensitivity)  
##########################################################################
fu <-glm(steroids~ UniRef90_C7H1G6 , data=df_steroids, family=binomial) #stronger
hu <-glm(steroids~ UniRef90_R5CY66 , data=df_steroids, family=binomial) #stronger
iu <-glm(steroids~ UniRef90_R6CZ24 , data=df_steroids, family=binomial) #stronger
ku <-glm(steroids~ UniRef90_T5S060 , data=df_steroids, family=binomial) #interesting this was null 

#MV analysis (without host genetics) #NB: (!) means neg ass'n # means null 
##########################################################################
ff <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_C7H1G6 , data=df_steroids, family=binomial) 
hh <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_R5CY66 , data=df_steroids, family=binomial) 
ii <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_R6CZ24 , data=df_steroids, family=binomial)
kk <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_T5S060 , data=df_steroids, family=binomial) 
#aa <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_A0A076IMG8, data=df_steroids, family=binomial) 
#bb <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_A0A076IRG9 , data=df_steroids, family=binomial) 
#cc <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_A0A076IXU4 , data=df_steroids, family=binomial) 
#dd <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_A0A0P0M2W7 , data=df_steroids, family=binomial) 
#(!)ee <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_A0A1C7H380 , data=df_steroids, family=binomial) 
#gg <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_I9R3P9 , data=df_steroids, family=binomial) 
#jj <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_R6TIX3 , data=df_steroids, family=binomial) 
#ll <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_U2QX51, data=df_steroids, family=binomial) 
#ecoli <-glm(steroids~ consent_age + sex + smoking.status + diagnosis + UniRef90_P77567, data=df_steroids, family=binomial) #null here, although not in the NAT group 


### check if these are additive eg if there's a score (without host genetics) #####
##########################################################################

df_steroids_score <- df_steroids 
df_steroids_score <- df_steroids_score%>% 
  mutate_at(vars(UniRef90_C7H1G6 , UniRef90_R5CY66 , UniRef90_R6CZ24 , UniRef90_T5S060), ~as.numeric(as.character(.))) 

df_steroids_score$score <- apply(df_steroids_score[,c(39,42,43,45)], 1, sum) 

#pt 1: look at 5 categories (0,1,2,3,4)
score_lm <- glm(steroids~ consent_age + sex + smoking.status + diagnosis + as.factor(score), data=df_steroids_score, family=binomial) #gives by level 
summary(glm(steroids~ consent_age + sex + smoking.status + diagnosis + score, data=df_steroids_score, family=binomial)) #gives trend

#pt 2: look at 5 categories (0,1,2,3+)
df_steroids_score<- df_steroids_score  %>% mutate(score2 = ifelse(score > 2, 3, score)) 
score_lm2 <- glm(steroids~ consent_age + sex + smoking.status + diagnosis + as.factor(score2), data=df_steroids_score, family=binomial) #gives by level 

#pt 3: look at 2 categories (0-2,3+)
df_steroids_score<- df_steroids_score  %>% mutate(score3 = ifelse(score > 2, 1, 0)) 
score_lm3 <- glm(steroids~ consent_age + sex + smoking.status + diagnosis + as.factor(score3), data=df_steroids_score, family=binomial) #gives by level 


#### now add host acetylation  ######
##########################################################################
####### NB exome sequencing data processed at broad servers for privacy reasons ########
snpselect <- read.table("/home/rsm34/5asabb3/genomics/snpselect.txt",header=T,sep="\t")
snpselect2 <- inner_join(snpselect,metadata2 ,by="SampleID") %>% select(Participant.ID,SampleID,acetylhost,snpsum,rs1041983,rs1801280,data_type,diagnosis)
snpselect3 <- snpselect2 %>% distinct(Participant.ID,.keep_all=TRUE) %>% rename(Participant.ID.x=Participant.ID)

#find out if ever steroids 
#find out if ever acetyltransferase microbe
df_steroids_ever <- df_steroids_score %>% mutate(bugacetyl = ifelse(score < 1, 0, 1)) 

df_steroids_ever2 <- df_steroids_ever %>% 
	group_by(Participant.ID.x) %>% 
		summarize(steroid_ever=max(as.numeric(steroids)-1),bugacetyl_ever=max(score3),consent_age=max(consent_age),sex=unique(sex),smoking.status=unique(smoking.status),diagnosis=unique(diagnosis))

df_steroids_ever2 <- df_steroids_ever2 %>% distinct(Participant.ID.x,.keep_all=TRUE) #do this b/c there are some m/f coding errors 

df_steroids_ever3 <- inner_join(df_steroids_ever2,snpselect3 ,by="Participant.ID.x")
df_steroids_ever3$acetylhost <- as.factor(df_steroids_ever3$acetylhost  ) 


#### makes proportions pie charts ####################
##########################################################################
### graph the results below #

#get proportions of those with 0-2 vs 3+ (who I will call fast vs slow) check and "fast vs slow acetylators"
table(df_steroids_score$score3) 
prop_bugs <- c(68,541)
pie(prop_bugs, labels = c("11% (n=68), Fast","89% (n=541), Slow"))

table(df_steroids_ever3$acetylhost) #14 fast,25 slow
df_steroids_ever3$acetylhost <- relevel(df_steroids_ever3$acetylhost, ref = "slow") 
Prop <- c(14,25)
pie(Prop , labels = c("38%, Fast","62%, Slow"))

#ever steroids vs single exposure (sensitivity) + baseline microbe prediction ==> these are NULL, maybe discouraging, but see mixed effects model below  
##########################################################################################################################################################
summary(glm(steroid_ever~ acetylhost + consent_age + sex + smoking.status + diagnosis.x, data=df_steroids_ever3, family=binomial)) # this is null  n=39
# > acetylhostfast      1.612e+00  1.056e+00   1.527   0.1267
summary(glm(steroid_ever~ bugacetyl_ever+ consent_age + sex + smoking.status + diagnosis.x, data=df_steroids_ever3, family=binomial)) #this is null n=39
# > bugacetyl_ever        0.53890    0.84808   0.635   0.5251

df_steroids_ever_bl <- df_steroids_ever %>% 
	group_by(Participant.ID.x) %>% 
  mutate(first_visit = min(visit_num)) %>% 
  ungroup() %>% 
  filter(visit_num==first_visit) %>% select(Participant.ID.x,SampleID, visit_num,sex,consent_age,smoking.status,diagnosis,score3,score,score2) %>% 
	distinct(Participant.ID.x,.keep_all=TRUE) %>% mutate(score_bl1 = ifelse(score < 1, 0, 1))  %>% mutate(score_bl2 = ifelse(score < 2, 0, 1)) 

everster <- df_steroids_ever3 %>% select(Participant.ID.x,steroid_ever)

baseline<- inner_join(everster ,df_steroids_ever_bl,by="Participant.ID.x")

table(baseline$steroid_ever,baseline$score_bl2)
#     0  1
#  0 19  7
#  1  9  4
table(baseline$steroid_ever,baseline$score_bl1)
#     0  1
#  0 10 16
#  1  4  9

summary(glm(steroid_ever~ score_bl1+ consent_age + sex + smoking.status + diagnosis, data=baseline, family=binomial)) #this is null n=39
summary(glm(steroid_ever~ score_bl2 + consent_age + sex + smoking.status + diagnosis, data=baseline, family=binomial)) #this is null n=39
#score_bl2             0.51143    1.00321   0.510   0.6102
summary(glm(steroid_ever~ score + consent_age + sex + smoking.status + diagnosis, data=baseline, family=binomial)) #this is null n=39


##########################################################################
#### now look at multiple time points -- MAIN RESULTS ####################
##########################################################################


df_steroids_ever4 <- inner_join(df_steroids_ever,snpselect3 ,by="Participant.ID.x")
df_steroids_ever4$acetylhost <- as.factor(df_steroids_ever4$acetylhost  ) 
df_steroids_ever4$acetylhost <- relevel(df_steroids_ever4$acetylhost, ref = "slow") 

#mutual adjustment (sensitivity) 
bbbb<-glm(steroids~ consent_age + sex + smoking.status +diagnosis.x + UniRef90_C7H1G6 + UniRef90_R5CY66 + UniRef90_R6CZ24 + UniRef90_T5S060 +acetylhost, data=df_steroids_ever4, family=binomial) 

#one by one for panel A
#####################
bbbb1<-glm(steroids~ consent_age + sex + smoking.status +diagnosis.x + UniRef90_C7H1G6 +acetylhost, data=df_steroids_ever4, family=binomial) 
bbbb2<-glm(steroids~ consent_age + sex + smoking.status +diagnosis.x + UniRef90_R5CY66 +acetylhost, data=df_steroids_ever4, family=binomial) 
bbbb3<-glm(steroids~ consent_age + sex + smoking.status +diagnosis.x + UniRef90_R6CZ24 +acetylhost, data=df_steroids_ever4, family=binomial) 
bbbb4<-glm(steroids~ consent_age + sex + smoking.status +diagnosis.x + UniRef90_T5S060 +acetylhost, data=df_steroids_ever4, family=binomial) 
ecoli2<-glm(steroids~ consent_age + sex + smoking.status +diagnosis.x + UniRef90_P77567 +acetylhost, data=df_steroids_ever4, family=binomial) 

df_steroids_ever4<- df_steroids_ever4  %>% mutate(score_bin = ifelse(score > 2, 1, 0)) 

#for panel C to compare host vs microbial genetics 
###################################################
bbbb5 <- glm(steroids~ consent_age + sex + smoking.status + diagnosis.x + as.factor(score)+acetylhost, data=df_steroids_ever4, family=binomial) #gives by level 
bbbb6 <- glm(steroids~ consent_age + sex + smoking.status + diagnosis.x + as.factor(score_bin)+acetylhost, data=df_steroids_ever4, family=binomial) #gives binary 

############ sensitivity analysis, seeing if mixed effects (Id as fixed effects) ######## 
### NB: get warning:https://stackoverflow.com/questions/38997371/glmer-warning-variance-covariance-matrix-is-not-positive-definite-or-cont
# "This warning suggests that your standard error estimates might be less accurate."
# "The estimates (point estimates, SEs, Z values ...) have changed a little bit, but not very much, so at the very least that should reassure you."

library(lme4)
bbbb7 <- glmer(steroids~ consent_age + sex + smoking.status + diagnosis.x + as.factor(score_bin)+acetylhost+ (1|Participant.ID.x), data=df_steroids_ever4, family=binomial) #gives binary 
exp(1.353)


summary(bbbb7)
#                        Estimate Std. Error z value Pr(>|z|)
#(Intercept)           -3.500e+00  1.969e+00  -1.777   0.0755 .
#consent_age           -1.539e-01  9.534e-02  -1.614   0.1065
#sexMale               -7.744e-02  6.337e-01  -0.122   0.9027
#smoking.statusever    -1.436e+02  8.455e+06   0.000   1.0000
#diagnosis.xCD          2.889e+00  1.628e+00   1.774   0.0761 .
#as.factor(score_bin)1  1.353e+00  6.809e-01   1.987   0.0469 *
#acetylhostfast         3.490e+00  1.743e+00   2.003   0.0452 *

#i think these are less reliable because of this warning "Hessian is numerically singular: parameters are not uniquely determined"
mixed1 <- glmer(steroids~ consent_age + sex + smoking.status + diagnosis.x + as.factor(UniRef90_C7H1G6)+acetylhost+ (1|Participant.ID.x), data=df_steroids_ever4, family=binomial) #borderline 
#as.factor(UniRef90_C7H1G6)1  7.250e-01  4.355e-01   1.665   0.0960 .
#acetylhostfast               3.576e+00  1.772e+00   2.018   0.0436 *

mixed2 <- glmer(steroids~ consent_age + sex + smoking.status + diagnosis.x + as.factor(UniRef90_R5CY66)+acetylhost+ (1|Participant.ID.x), data=df_steroids_ever4, family=binomial) #null
#as.factor(UniRef90_R5CY66)1  2.077e-01  4.460e-01   0.466   0.6414
#acetylhostfast               3.546e+00  1.798e+00   1.973   0.0485 *

mixed3 <- glmer(steroids~ consent_age + sex + smoking.status + diagnosis.x + as.factor(UniRef90_R6CZ24 )+acetylhost+ (1|Participant.ID.x), data=df_steroids_ever4, family=binomial) #sig 
#as.factor(UniRef90_R6CZ24)1  1.615e+00  6.221e-01   2.597  0.00941 **
#acetylhostfast               3.564e+00  1.858e+00   1.918  0.05511 .

mixed4 <- glmer(steroids~ consent_age + sex + smoking.status + diagnosis.x + as.factor(UniRef90_T5S060 )+acetylhost+ (1|Participant.ID.x), data=df_steroids_ever4, family=binomial) #null 
#as.factor(UniRef90_T5S060)1  6.045e-01  6.638e-01   0.911   0.3624
#acetylhostfast               3.552e+00  1.773e+00   2.004   0.0451 *


##############################################
#### EXPORT MAIN RESULTS ####################
##############################################

#### create dataframe to store these values for STEROIDS (pull age and CD from mutually adjusted models below) ### 

ORacetyl_table = NULL
ORacetyl_table = as.data.frame(matrix(-9,14,4))
colnames(ORacetyl_table) = c("Uniref","RiskRatio","LowerLimit","UpperLimit")
ORacetyl_table[1,1]<-"referent"
ORacetyl_table[1,2:4]<- c(1,1,1)
ORacetyl_table[8,2:4]<- c(1,1,1)
ORacetyl_table[2,1] <- "C7H1G6"
ORacetyl_table[3,1] <- "R5CY66"
ORacetyl_table[4,1] <- "R6CZ24"
ORacetyl_table[5,1] <- "T5S060"
ORacetyl_table[6,1] <- "age, y"
ORacetyl_table[7,1] <- "CD (vs UC)"
ORacetyl_table[8,1] <- "Zero (n=269)"
ORacetyl_table[9,1] <- "One (n=158)"
ORacetyl_table[10,1] <- "Two (n=114)"
ORacetyl_table[11,1] <- "Three (n=49)"
ORacetyl_table[12,1] <- "Four (n=19)"
ORacetyl_table[13,1] <- "Human NAT2 (Fast [0-1 SNP] v Slow [2+ SNP])"
ORacetyl_table[14,1] <- "Bacterial AT (Fast [3+] v Slow [0-2])"

#have to turn coefficients to ORs with 95% ci 
ORacetyl_table[2,2] <- exp(summary(bbbb1)$coefficients[6,1])
ORacetyl_table[2,3] <- exp(summary(bbbb1)$coefficients[6,1]-(1.96*summary(bbbb1)$coefficients[6,2]))
ORacetyl_table[2,4] <- exp(summary(bbbb1)$coefficients[6,1]+(1.96*summary(bbbb1)$coefficients[6,2]))

ORacetyl_table[3,2] <- exp(summary(bbbb2)$coefficients[6,1])
ORacetyl_table[3,3] <- exp(summary(bbbb2)$coefficients[6,1]-(1.96*summary(bbbb2)$coefficients[6,2]))
ORacetyl_table[3,4] <- exp(summary(bbbb2)$coefficients[6,1]+(1.96*summary(bbbb2)$coefficients[6,2]))

ORacetyl_table[4,2] <- exp(summary(bbbb3)$coefficients[6,1])
ORacetyl_table[4,3] <- exp(summary(bbbb3)$coefficients[6,1]-(1.96*summary(bbbb3)$coefficients[6,2]))
ORacetyl_table[4,4] <- exp(summary(bbbb3)$coefficients[6,1]+(1.96*summary(bbbb3)$coefficients[6,2]))

ORacetyl_table[5,2] <- exp(summary(bbbb4)$coefficients[6,1])
ORacetyl_table[5,3] <- exp(summary(bbbb4)$coefficients[6,1]-(1.96*summary(bbbb4)$coefficients[6,2]))
ORacetyl_table[5,4] <- exp(summary(bbbb4)$coefficients[6,1]+(1.96*summary(bbbb4)$coefficients[6,2]))

ORacetyl_table[6,2] <- exp(summary(bbbb)$coefficients[2,1])
ORacetyl_table[6,3] <- exp(summary(bbbb)$coefficients[2,1]-(1.96*summary(bbbb)$coefficients[2,2]))
ORacetyl_table[6,4] <- exp(summary(bbbb)$coefficients[2,1]+(1.96*summary(bbbb)$coefficients[2,2]))

ORacetyl_table[7,2] <- exp(summary(bbbb)$coefficients[5,1])
ORacetyl_table[7,3] <- exp(summary(bbbb)$coefficients[5,1]-(1.96*summary(bbbb)$coefficients[5,2]))
ORacetyl_table[7,4] <- exp(summary(bbbb)$coefficients[5,1]+(1.96*summary(bbbb)$coefficients[5,2]))
##
ORacetyl_table[9,2] <- exp(summary(bbbb5)$coefficients[6,1])
ORacetyl_table[9,3] <- exp(summary(bbbb5)$coefficients[6,1]-(1.96*summary(bbbb5)$coefficients[6,2]))
ORacetyl_table[9,4] <- exp(summary(bbbb5)$coefficients[6,1]+(1.96*summary(bbbb5)$coefficients[6,2]))

ORacetyl_table[10,2] <- exp(summary(bbbb5)$coefficients[7,1])
ORacetyl_table[10,3] <- exp(summary(bbbb5)$coefficients[7,1]-(1.96*summary(bbbb5)$coefficients[7,2]))
ORacetyl_table[10,4] <- exp(summary(bbbb5)$coefficients[7,1]+(1.96*summary(bbbb5)$coefficients[7,2]))


ORacetyl_table[11,2] <- exp(summary(bbbb5)$coefficients[8,1])
ORacetyl_table[11,3] <- exp(summary(bbbb5)$coefficients[8,1]-(1.96*summary(bbbb5)$coefficients[8,2]))
ORacetyl_table[11,4] <- exp(summary(bbbb5)$coefficients[8,1]+(1.96*summary(bbbb5)$coefficients[8,2]))


ORacetyl_table[12,2] <- exp(summary(bbbb5)$coefficients[9,1])
ORacetyl_table[12,3] <- exp(summary(bbbb5)$coefficients[9,1]-(1.96*summary(bbbb5)$coefficients[9,2]))
ORacetyl_table[12,4] <- exp(summary(bbbb5)$coefficients[9,1]+(1.96*summary(bbbb5)$coefficients[9,2]))
##
ORacetyl_table[13,2] <- exp(summary(bbbb6)$coefficients[7,1])
ORacetyl_table[13,3] <- exp(summary(bbbb6)$coefficients[7,1]-(1.96*summary(bbbb6)$coefficients[7,2]))
ORacetyl_table[13,4] <- exp(summary(bbbb6)$coefficients[7,1]+(1.96*summary(bbbb6)$coefficients[7,2]))

ORacetyl_table[14,2] <- exp(summary(bbbb6)$coefficients[6,1])
ORacetyl_table[14,3] <- exp(summary(bbbb6)$coefficients[6,1]-(1.96*summary(bbbb6)$coefficients[6,2]))
ORacetyl_table[14,4] <- exp(summary(bbbb6)$coefficients[6,1]+(1.96*summary(bbbb6)$coefficients[6,2]))

exp(summary(ecoli2)$coefficients[6,1]) #this is significant, but without host genetics is not 
exp(summary(ecoli2)$coefficients[6,1]-(1.96*summary(ecoli2)$coefficients[6,2]))
exp(summary(ecoli2)$coefficients[6,1]+(1.96*summary(ecoli2)$coefficients[6,2]))


#summary(bbbb) --> mutually adjusted (sensitivity)  
#UniRef90_C7H1G6      0.68612    0.28426   2.414   0.0158 *
#UniRef90_R5CY66      0.71975    0.29885   2.408   0.0160 *
#UniRef90_R6CZ24      0.35502    0.34484   1.030   0.3032
#UniRef90_T5S060      0.86607    0.37495   2.310   0.0209 *

### write these out for R studio plotting (see below for plotting code) 
write.table(ORacetyl_table, file="ORacetyl_table.txt", sep="\t", quote=FALSE,row.names=FALSE)


#### get proportions for users ####
##############################################

prop.table(table(df_steroids$UniRef90_C7H1G6))
prop.table(table(df_steroids$UniRef90_R6CZ24))
prop.table(table(df_steroids$UniRef90_R5CY66))
prop.table(table(df_steroids$UniRef90_T5S060))

#now input from other files (nonuserIBD) and make a small bar chart (won't use nonIBD)

library(ggplot2)
library(RColorBrewer)

group_names <- c(
  "Non-IBD (n=384)", "5-ASA non-users (n=382 samples)", "5-ASA users (n=609 samples)")
gene_pct <- data.frame(
  Group = factor(group_names, levels = group_names),
  C7H1G6 = c(67, 52, 48),
  R6CZ24 = c(29, 20, 17),
  R5CY66 = c(44, 14, 21),
  T5S060 = c(12, 09, 14)
)

#users 
#> prop.table(table(df_steroids$UniRef90_C7H1G6))
#        0         1
#0.5221675 0.4778325
#> prop.table(table(df_steroids$UniRef90_R6CZ24))
#        0         1
#0.8259442 0.1740558
#> prop.table(table(df_steroids$UniRef90_R5CY66))
#        0         1
#0.7931034 0.2068966
#> prop.table(table(df_steroids$UniRef90_T5S060))
#        0         1
#0.8587849 0.1412151

#IBD nonusers
#> prop.table(table(df_steroids$UniRef90_C7H1G6))
#        0         1
#0.4768519 0.5231481
#> prop.table(table(df_steroids$UniRef90_R6CZ24))
#        0         1
#0.7962963 0.2037037
#> prop.table(table(df_steroids$UniRef90_R5CY66))
#        0         1
#0.8564815 0.1435185
#> prop.table(table(df_steroids$UniRef90_T5S060))
#         0          1
#0.90740741 0.09259259

genepct_gather <- gene_pct %>% gather(gene,value,2:5) %>% filter(Group=="5-ASA non-users (n=382 samples)"  | Group=="5-ASA users (n=609 samples)")

### write these out for R studio plotting (see below for plotting code) 
write.table(genepct_gather , file="genepct_gather.txt", sep="\t", quote=FALSE,row.names=FALSE)


#################################################### GRAPHICS CODE ########################################################

library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(cowplot)

RR_data <- ORacetyl_table
RR_data$Uniref <- as.character(RR_data$Uniref)
RR_data[13,1] <- "Human NAT2 \n(Fast [0-1 SNP] \nv Slow [2+ SNP])"
RR_data[14,1] <- "Bacterial AT \n(Fast [3+] \nv Slow [0-2])"


blues <- brewer.pal(n=9,name="PuBu")
dkblue <- blues[8]
ltblue <- blues[6]

### PANEL A ######
RR_data1 <- RR_data[1:5,]


RR_data1$Uniref <- factor(RR_data1$Uniref , levels=c("referent", "C7H1G6", "R5CY66","R6CZ24","T5S060"))


pq = ggplot(data=RR_data1, aes(x = Uniref))+ geom_pointrange(aes(y = RiskRatio, ymin = LowerLimit, ymax = UpperLimit),size=2,fill="#DC0000FF",shape=22,colour = "grey47")+
    geom_hline(yintercept =1, linetype=2)+
    xlab('Uniref90 Name')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))  +
    scale_y_continuous(name = "Odds Ratio (95% Confidence Interval)",expand = expansion(mult = c(0.25, 0))) + 
    geom_text(aes(y = RiskRatio, label = round(RiskRatio,2)),hjust=0,vjust=1,nudge_x=0.2,size=5) +
    #scale_x_continuous(expand = expansion(mult = c(0.05, .1)),breaks=c(0,1,2,3,4,5),labels=RR_data$Group) + 
    theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.text=element_text(size=12), axis.title.x = element_text(size=13), axis.title.y = element_text(size=13))+theme_cowplot(14)+coord_trans(y = "log10")+ coord_cartesian(clip = 'off') 
pq

ggsave(filename='/home/rsm34/5asabb3/submission/figures/fig5a.png', plot=pq , width = 4.5, height = 5, dpi = 600) 


#4.5 x 6
### SCORE PANEL SUPP? ######
RR_data2 <- RR_data[8:12,]

RR_data2$Uniref <- factor(RR_data2$Uniref , levels=c("Zero (n=269)", "One (n=158)", "Two (n=114)","Three (n=49)","Four (n=19)"))


pq3 = ggplot(data=RR_data2, aes(x = Uniref))+ geom_pointrange(aes(y = RiskRatio, ymin = LowerLimit, ymax = UpperLimit),size=1.5,color="black",fill="firebrick3",shape=22)+
    geom_hline(yintercept =1, linetype=2)+
    xlab('Number of genes')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))  +
    scale_y_continuous(name = "Odds Ratio (95% Confidence Interval)",expand = expansion(mult = c(0.25, 0))) + 
    geom_text(aes(y = RiskRatio, label = round(RiskRatio,2)),hjust=0,vjust=1,nudge_x=0.2,size=4) +
    #scale_x_continuous(expand = expansion(mult = c(0.05, .1)),breaks=c(0,1,2,3,4,5),labels=RR_data$Group) + 
    theme(axis.text.x = element_text(size=11),axis.text.y = element_text(size=11),axis.text=element_text(size=11), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))+coord_cartesian(ylim=c(0, 23))
pq3



### PANEL C --> will probably nix ######


RR_data3 <- RR_data[c(13:14),]

pq3 = ggplot(data=RR_data3, aes(x = Uniref))+ geom_pointrange(aes(y = RiskRatio, ymin = LowerLimit, ymax = UpperLimit),size=1.5,color="black",fill="firebrick",shape=22)+
    geom_hline(yintercept =1, linetype=2)+
    xlab('')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))  +
    scale_y_continuous(name = "Odds Ratio (95% Confidence Interval)",expand = expansion(mult = c(0.25, 0))) + 
    geom_text(aes(y = RiskRatio, label = round(RiskRatio,2)),hjust=0,vjust=1,nudge_x=0.2,size=4) +
    #scale_x_continuous(expand = expansion(mult = c(0.05, .1)),breaks=c(0,1,2,3,4,5),labels=RR_data$Group) + 
    theme(axis.text.x = element_text(size=11),axis.text.y = element_text(size=11),axis.text=element_text(size=11), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))+coord_cartesian(ylim=c(0, 17.6))
pq3

#possibly to accompany panel C
#68 "fast" 541 "slow" as above 
prop_bugs <- c(68,541)
pie(prop_bugs, labels = c("11% (n=68),\n Fast","89% (n=541),\n Slow"))

#14 fast,25 slow
Prop <- c(14,25)
pie(Prop , labels = c("36% (n=14), Fast","64% (n=25), Slow"))


################################################# PANEL B: get proportions of presence /absence genes ######################


pq4 <- ggplot(data=genepct_gather, aes(x=gene,y=value,fill=factor(Group))) +
  geom_bar(position="dodge",stat="identity",width=0.7) + theme_bw() +
 scale_fill_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8)) +
  theme(legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
   text = element_text(size = 14),axis.text.x = element_text(angle = 45,hjust = 1),axis.title.x=element_blank()) + ylab("% prevalence")+theme_cowplot(14)+theme(legend.title = element_blank())+ xlab('Uniref90 Name')
pq4 

ggsave(filename='/home/rsm34/5asabb3/submission/figures/fig5b.png', plot=pq4 , width = 7, height = 5, dpi = 600) 


