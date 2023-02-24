######################################################################################################
# R programs for creating "metabolite" plots								#
# ============================================								#
# Purposes:  new use boxplots,    									#
#              networks, identificaiton of unannotated compounds, r2 parsing				#
######################################################################################################

#references: 
#https://yardstick.tidymodels.org/reference/roc_curve.html
#https://davetang.org/muse/2017/03/16/matrix-to-adjacency-list-in-r/
#https://www.r-graph-gallery.com/257-input-formats-for-network-charts.html
#http://mr.schochastics.net/netVizR.html check this out 
#https://kateto.net/networks-r-igraph
#https://rdrr.io/cran/KMDA/man/spearman.group.html
##https://data.library.virginia.edu/getting-started-with-factor-analysis/


############################################################### CODE ###################################################

library(dplyr)
library(purrr) 
#library(table1) 
library(tidyr)
library(magrittr)
#library(sm) 
library(ggplot2)
library(corrplot)
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
library(ggpubr)
library(viridis)

##### STEP 1: Identify changes in metabolites before and after adminstration ####

#read in "light file"  
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt",header=T,sep="\t") 


#reviewer fig 
revselect <- metab_select %>% filter(binary154==1)
rev1fig <- ggplot(revselect, aes(x=log(X154.0502_3.83+1), y=log(X196.0609_2.81+1))) + 
    geom_point(
        color="black",
        fill="#69b3a2",
        shape=22,
        alpha=0.5,
        size=3,
        stroke = 1
        )+theme_bw()+xlab("5-ASA levels")+ylab("N-Acetyl 5-ASA levels")

ggsave(filename='/home/rsm34/5asabb3/reviews/rev1fig.png', plot=rev1fig , width = 5, height = 5, dpi = 600) 


#read in metadata file 
meta_data <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t")
meta_data$SampleID <- as.character(meta_data$SampleID)
cols <- c("oral5asa_lc","diagnosis","biologics","site_name","Antibiotics","bowelsurg","Participant.ID","steroids","pr5asa","bond5asa","sex","dysbiosis")
meta_data[cols] <- lapply(meta_data[cols], factor) #converts these to a factor variable)
meta_data <- meta_data %>% select(c("SampleID","diagnosis","Participant.ID","dysbiosis","site_name","visit_num","oral5asa_lc","sex","consent_age","Antibiotics")) %>% filter(diagnosis == "CD" | diagnosis == "UC")  

mesal_metab<- inner_join(metab_select,meta_data,by="SampleID") %>% drop_na(diagnosis)
   	#write.table(mesal_metab, "newusersID.txt", sep="\t", quote=FALSE,row.names=FALSE) #this file helps us find new users shown below(based on binary154)


#I then selected participants who were new users of 5-asa based on metabolite, verified 12.17, can see below  
newids <- c("PSM7J19B", "PSM7J19J", "PSM6XBSE", "PSM6XBVM", "PSM6XBSK", "PSM6XBUG", "PSM6XBRK", "PSM7J1CU", "MSMA26AZ", "MSMB4LZ4",
"MSM6J2IG", "MSM6J2Q3", "MSM5LLDI", "MSM5LLDS", "HSMA33OJ", "HSMA33M8", "HSMA33OZ", "HSMA33MI", "HSM5MD73", "HSM6XRS8", "HSM6XRRV", 
"HSM6XRVO", "CSM79HRG", "CSM7KOO9", "CSM5MCXL", "CSM67UDN")


#subset
newiddata <- mesal_metab %>% filter(SampleID %in% newids)
newiddata %>% arrange(Participant.ID,binary154) %>% select(binary154,Participant.ID,visit_num,SampleID,X154.0502_3.83) #this shows the list n=13
#   binary154 Participant.ID visit_num SampleID X154.0502_3.83
#1          0          C3004        13 CSM5MCXL         283552
#2          1          C3004        20 CSM67UDN      507071413
#3          0          C3031        11 CSM79HRG         172164
#4          1          C3031        14 CSM7KOO9      129818420
#5          0          H4014        13 HSM6XRRV         290011
#6          1          H4014        19 HSM6XRVO     1530601385
#7          0          H4015         8 HSM5MD73         168868
#8          1          H4015        13 HSM6XRS8     3121179635
#9          0          H4035        23 HSMA33OZ        2352694
#10         1          H4035        29 HSMA33MI       82747705
#11         0          H4040        18 HSMA33OJ         390204
#12         1          H4040        29 HSMA33M8     1468731677
#13         0          M2008         4 MSM5LLDI        1773492
#14         1          M2008         9 MSM5LLDS      643195341
#15         0          M2028        13 MSM6J2IG         557418
#16         1          M2028        19 MSM6J2Q3      109827300
#17         0          M2071        12 MSMA26AZ        1036302
#18         1          M2071        16 MSMB4LZ4      323254933
#19         0          P6009         4 PSM6XBRK         169904
#20         1          P6009        25 PSM7J1CU      456525486
#21         0          P6010         8 PSM6XBSK         254310
#22         1          P6010        12 PSM6XBUG      879384985
#23         0          P6012         4 PSM6XBSE        1059704
#24         1          P6012         7 PSM6XBVM       12783515
#25         0          P6016         6 PSM7J19B         268885
#26         1          P6016        11 PSM7J19J     1171158465


#find out how many weeks on average between pre and post mesalamine (for text)
weeks <- read.csv("/home/rsm34/5asa/hmp2_metadata_load.csv",stringsAsFactors=FALSE) %>% select("External.ID","week_num") %>% rename(SampleID=External.ID)

weeks2 <- inner_join(newiddata,weeks,by="SampleID") %>% distinct(SampleID,.keep_all=TRUE)
weeks_pre <- weeks2 %>% filter(binary154==0) %>% arrange(Participant.ID) 
weeks_post <- weeks2 %>% filter(binary154==1) %>% arrange(Participant.ID) 
mean(weeks_post$week_num-weeks_pre$week_num) #13.1 
sd(weeks_post$week_num-weeks_pre$week_num) #8.7
weeks_post$diff <-  weeks_post$week_num-weeks_pre$week_num
weeks_diff <- weeks_post %>% select(Participant.ID,diff)
#read in FULL METABOLOMICS file after pre-processing script made it 
#########################
metab_rotated2 <- read.table("/home/rsm34/5asa/metab_data_nozero2.txt",header=T,sep="\t")

###### testing ##########

metabs_newids <- inner_join(newiddata,metab_rotated2,by="SampleID") %>% arrange(Participant.ID)%>% arrange(binary154)

#paired wilcoxon for metabolite data  
df = NULL 
df <- metabs_newids[,15:ncol(metabs_newids)]

kmetabs2 = NULL
kmetabs2 = as.data.frame(matrix(-9,ncol(df),2))
colnames(kmetabs2) = c("feature","p")

for(i in c(1:ncol(df)))
	{	b<-wilcox.test(df[,i] ~ metabs_newids$binary154, paired=TRUE)
 		kmetabs2[i,1]= colnames(df)[i]
 		kmetabs2[i,2]=b$p.value
	}
	kmetabs2$FDR <- p.adjust(kmetabs2$p,method="fdr")
	kmetabs2 <- kmetabs2 %>% arrange((p))
   	#write.table(kmetabs2 , "newusersMetabs_all.txt", sep="\t", quote=FALSE,row.names=FALSE)
	newusersMetabs <- kmetabs2
newusersmetab <- newusersMetabs


#### sensitivity analysis for reviewer to show duration doesn't matter  ####
library(nlme)
head(names(metabs_newids))
metabs_newids <- metabs_newids %>% select(nicotinuric.acid,X2.aminoadipate,nicotinate ,X196.0609_2.81.x,binary154, Participant.ID)
metabs_newids2<-inner_join(metabs_newids,weeks_diff,by="Participant.ID")
b <- lme(nicotinuric.acid ~ binary154 + diff, random=~1|Participant.ID, data=metabs_newids2,na.action=na.omit,control=lmeControl(returnObject=TRUE))
a <- lme(X2.aminoadipate ~ binary154 + diff, random=~1|Participant.ID, data=metabs_newids2,na.action=na.omit,control=lmeControl(returnObject=TRUE))
c <- lme(nicotinate ~ binary154 + diff, random=~1|Participant.ID, data=metabs_newids2,na.action=na.omit,control=lmeControl(returnObject=TRUE))
d <- lme(X196.0609_2.81.x~ binary154 + diff, random=~1|Participant.ID, data=metabs_newids2,na.action=na.omit,control=lmeControl(returnObject=TRUE))


summary(b)
summary(a)
summary(c)
summary(d)


## for future plotting
metabs_plotting2 <- metabs_newids %>% select(matches("SampleID|Participant.ID|nicotinuric.acid|X196.0609_2.81.x|binary154|X1.methylnicotinamide|nicotinate|X154.0502_3.83.x|X312.073_6.98|X458.1946_1.76|X189.0771_3.98|X2.aminoadipate|X373.1254_0.71|X242.0458_3"))
names(metabs_plotting2) <- gsub(".x","",names(metabs_plotting2))
newusers_metabs <- metabs_plotting2 %>% arrange(binary154,Participant.ID)
#write.table(newusers_metabs,"newusers_metabs.txt",sep="\t",quote=FALSE,row.names=FALSE)


################################### 1) Fig 2B: PLOTTING FOR NEW USERS ###############################################
###### write a file out for the new users plot and for the uniref analyses ########

mesallabels <- c("pre", "post")


P2 <- ggpaired(newusers_metabs, x = "binary154", y = "nicotinuric.acid",
         color = "binary154", line.color = "gray", line.size = 0.4,
         palette = c("#3C5488FF", "#DC0000FF"),short.panel.labs = FALSE, xlab="5-ASA use",ylab="rel.abund.",title = "nicotinuric.acid",point.size = 1.5,legend="none") +theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(labels=mesallabels)

P3 <- ggpaired(newusers_metabs, x = "binary154", y = "X2.aminoadipate",
         color = "binary154", line.color = "gray", line.size = 0.4,
         palette = c("#3C5488FF", "#DC0000FF"),short.panel.labs = FALSE, xlab="5-ASA use",ylab="rel.abund.",title = "X2.aminoadipate",point.size = 1.5,legend="none") +theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(labels=mesallabels)

pairedplotfxn <- function(metabolite,caption) {
	ggpaired(newusers_metabs, x = "binary154", y = metabolite,
         color = "binary154", line.color = "gray", line.size = 0.4,
         palette = c("#3C5488FF", "#DC0000FF"),short.panel.labs = FALSE, xlab="5-ASA use",ylab="abundance",title = caption,point.size = 1.5,legend="none",font.x=12,font.y=12,font.tickslab=12,font.main=12) +theme(plot.title = element_text(hjust = 0.5))+
scale_x_discrete(labels=mesallabels)
}

library(gridExtra)
pp1<-pairedplotfxn("X154.0502_3.83","5-ASA")
pp2<-pairedplotfxn("X2.aminoadipate","2-aminoadipate")
pp3<-pairedplotfxn("X196.0609_2.81","N-Acetyl 5-ASA")

pp4<-pairedplotfxn("X312.073_6.98","312.073\n module:hippurate")
pp5<-pairedplotfxn("X458.1946_1.76","458.1946\n module:Glycerol 3-phosphate")
pp6<-pairedplotfxn("X189.0771_3.98","189.0771\nmodule:Lithocholate")

pp7<-pairedplotfxn("nicotinuric.acid","Nicotinuric Acid")
pp8<-pairedplotfxn("X1.methylnicotinamide","N1-Methyl-nicotinamide")
pp9<-pairedplotfxn("nicotinate","Nicotinate")


fig2b <-grid.arrange(pp1,pp3,pp2,pp7,pp9,pp8,nrow=2) #5.8x6.9 pdf 

ggsave(filename='/home/rsm34/5asabb3/submission/figures/fig2b.png', plot=fig2b, width = 6.9, height = 5.8, dpi = 600) 

##########FIG S5: now get a time series for these people #### 
#subset to IBD users only 
meta_data_IBD <- meta_data %>% filter(diagnosis=="CD"|diagnosis=="UC") 

raaj<- inner_join(meta_data_IBD,metab_rotated2,by="SampleID") %>% filter(SampleID %in% newids) %>% select(X154.0502_3.83,SampleID,Participant.ID)

newusercommonID<-unique(raaj$Participant.ID)

timegraph_nico<- inner_join(meta_data_IBD,metab_rotated2,by="SampleID") %>% filter(Participant.ID %in% newusercommonID)%>% select(matches("nicotinuric.acid|X154.0502_3.83|X312.073_6.98|X458.1946_1.76|X189.0771_3.98|X196.0609_2.81|Participant.ID|visit_num"))

timegraph_nico <- timegraph_nico %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) 

timegraph_H4035 <- timegraph_nico %>% filter(Participant.ID=="H4035") #these are selected based on excel selection from metadata files showing 1) stoppage 2)starting or 3) both
timegraph_P6016 <- timegraph_nico %>% filter(Participant.ID=="P6016")
timegraph_P6010 <- timegraph_nico %>% filter(Participant.ID=="P6010")

timegraph_H4035g <- gather(timegraph_H4035, metab, value, nicotinuric.acid:X154.0502_3.83) %>% mutate(value = log10(value+1))
timegraph_P6016g <- gather(timegraph_P6016, metab, value, nicotinuric.acid:X154.0502_3.83) %>% mutate(value = log10(value+1))
timegraph_P6010g <- gather(timegraph_P6010, metab, value, nicotinuric.acid:X154.0502_3.83) %>% mutate(value = log10(value+1))

timegraph_P6010g$metab <- as.character(timegraph_P6010g$metab)

timegraph_P6010g$metab[timegraph_P6010g$metab=="X458.1946_1.76"] <- "458.1946"
timegraph_P6010g$metab[timegraph_P6010g$metab=="X189.0771_3.98"] <- "189.0771"
timegraph_P6010g$metab[timegraph_P6010g$metab=="X312.073_6.98"] <- "312.073"
timegraph_P6010g$metab[timegraph_P6010g$metab=="X196.0609_2.81"] <- "N-ac-5-ASA"
timegraph_P6010g$metab[timegraph_P6010g$metab=="X154.0502_3.83"] <- "5-ASA"
timegraph_P6010g$metab[timegraph_P6010g$metab=="nicotinuric.acid"] <- "Nicotinuric acid"


# Plot x 3 
figs5 <-ggplot(timegraph_P6010g, aes(x=visit_num, y=value, fill=metab)) + 
    geom_area(alpha=0.6 , size=.5, colour="white") +
    scale_fill_viridis(discrete = T) + theme_bw() +
    ggtitle("P6010") +
  scale_x_continuous(breaks = c(0,5,8,12,19,26,30),expand=c(0,0)) + expand_limits(x=0)  


figs5 <- ggplot(
  timegraph_P6010g,
  aes(visit_num, value, group = metab, color = factor(metab))
  ) +
  geom_line() +
  scale_color_viridis_d() +
  labs(x = "Visit Number", y = "value, log()") +
  theme(legend.position = "top") + theme_bw()+ expand_limits(x=0)  +
  scale_x_continuous(breaks = c(0,5,8,12,19,26,30)) +  labs(color="Metabolite")

ggsave(filename='/home/rsm34/5asabb3/submission/figures/suppfig5.png', plot=figs5, width = 7.5, height = 5, dpi = 600) 


############## #STEP 2 ###################
### subset list and merge with metadata to get IBD only  ###############################
#new users 
newusersmetab2 <- newusersmetab %>% filter(FDR < 0.25) 
newusersmetab2$feature <- gsub(".y","",newusersmetab2$feature)
dim(newusersmetab2) #2,306

#subset results based on this 
###############################
selectvars <- newusersmetab2 %>% select(feature)    
selectvars <- rbind(c("SampleID"),c("X154.0502_3.83"),selectvars) 

metab_selected <- metab_rotated2 %>% select(any_of(selectvars$feature)) %>% mutate(binary196 = if_else(X196.0609_2.81 > 10000000,1,0)) %>% relocate(binary196, .after = SampleID) %>% 
	mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% relocate(binary154, .after = binary196)

#prepare metadata (will read it in again to be sure) 
meta_data <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t")
cols <- c("oral5asa_lc","diagnosis","biologics","site_name","Antibiotics","bowelsurg","dysbiosis","Participant.ID","steroids","pr5asa","bond5asa")
meta_data[cols] <- lapply(meta_data[cols], factor) #converts these to a factor variable)

meta_data$dysbiosis <- relevel(meta_data$dysbiosis, ref = "FALSE") # Set 'inactive' as reference
meta_data$diagnosis <- relevel(meta_data$diagnosis, ref = "nonIBD") # Set nonIBD as reference
meta_data <- meta_data %>% mutate(IBDusers = if_else(oral5asa_lc == 1,2,if_else(oral5asa_lc==0 & as.numeric(diagnosis) == 3,1, if_else(oral5asa_lc==0 & as.numeric(diagnosis) == 2,1,0)))) 
meta_data$IBDusers[meta_data$IBDusers=="2"] <- "User"
meta_data$IBDusers[meta_data$IBDusers=="1"] <- "IBDNon-User"
meta_data$IBDusers[meta_data$IBDusers=="0"] <- "nonIBD"
meta_data$IBDusers <- as.factor(meta_data$IBDusers)
meta_data$IBDusers <- relevel(meta_data$IBDusers, ref = "nonIBD") # Set nonIBD as reference

#subset to IBD users only 
meta_data_IBD <- meta_data %>% filter(diagnosis=="CD"|diagnosis=="UC") 

mesal_metab<- inner_join(meta_data_IBD,metab_selected,by="SampleID") %>% drop_na(diagnosis)

mesal_metab$binary196 <- as.factor(mesal_metab$binary196)
mesal_metab$binary154 <- as.factor(mesal_metab$binary154)

dim(mesal_metab) 
#283 x 2336 

#find out how many are non-zero####
medianzero <- mesal_metab  %>% group_by(binary154) %>% summarise(across(X2.aminoadipate:X299.1042_1.73, median, na.rm= TRUE)) 
medianzero2 <- as.data.frame(t(medianzero)) %>% arrange(V1,V2)
colnames(medianzero2) <- c("nonuser","user")
medianzero2$nonuser <- as.numeric(medianzero2$nonuser)
medianzero2$user <- as.numeric(medianzero2$user)
(sum(medianzero2$nonuser==0)-1)/(nrow(medianzero2)-1) #24.7%
medianzero2 <- medianzero2 %>% arrange(user) 
medianzero2$foldchange <- medianzero2$user/(medianzero2$nonuser+1)
medianzero2$id <- row.names(medianzero2)
medianzero2 %>% filter(id == "X373.1254_0.71")#73
medianzero2 %>% filter(id == "X242.0458_3")#37661

medianzero2$feature <- rownames(medianzero2) 
medianzero2 <- medianzero2 %>% relocate(feature , .before = nonuser)
rownames(medianzero2) = NULL 
medianzero2 <- medianzero2 %>% arrange(user) 
(sum(medianzero2$foldchange>1000))/(nrow(medianzero2)-1) #28.9% (number with fold change >1000) 
(sum(medianzero2$user>1000000))/(nrow(medianzero2)-1) #8.0% (relative abundance close >1000000) 



#################### STEP 3: GENERATE ROC CURVES WITH AUC ESTIMATION AS A WAY TO SUBSET ##################################################


metab_gather <- gather(mesal_metab, metab, value, X2.aminoadipate:X299.1042_1.73)

metabAUCs <- metab_gather %>%
  group_by(metab) %>%
  roc_auc(binary154, value, event_level = ("second")) %>% 
	arrange(desc(.estimate))

#write.table(metabAUCs,"metabAUCs_all.txt",sep="\t",quote=FALSE,row.names=FALSE)


############## first subset to bottom 20% (these are anticorrelated with 5-ASA use) #########
#subset to bottom 20% 
dim(metabAUCs) #2303
metabAUC20 <- metabAUCs %>% filter(.estimate<0.20) #250
dim(metabAUC20)

candidates_low <- c("X154.0502_3.83","X175.0614_3.96","nicotinate","nicotinuric.acid","X288.0487_2.25","X200.0922_0.95","X257.1136_3.75","SampleID","X337.141_1.36","X189.0771_3.98","X349.1772_5.83")
metab_low <- metab_rotated2 %>% select(any_of(candidates_low)) %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% relocate(binary154, .after = SampleID)
metab_low2 <- inner_join(meta_data_IBD,metab_low,by="SampleID") %>% drop_na(diagnosis) %>% select(1,31:41) %>% relocate(binary154, .after = SampleID)
metab_low3 <- gather(metab_low2, metab, value, X154.0502_3.83:X349.1772_5.83)

#graph a few to see what they look like (supplement?) 
#random<-metab_low3 %>%
# ggplot(aes(x=as.factor(binary154),y=value+1)) +
#  geom_boxplot() + 
#  theme_bw() +
#  scale_fill_viridis_d() +
#  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
#  facet_wrap(~ metab) + scale_y_log10()

#ggsave(filename='/home/rsm34/5asabb3/submission/figures/random.png', plot=random, width = 7.5, height = 5, dpi = 600) 


#example of a feature that is poorly correlated with drug levels (for supplement)? 
#metab_low3 %>% filter(metab == "X288.0487_2.25") %>% 
#	 ggplot(aes(x=as.factor(binary154),y=value+1)) +
#  geom_boxplot() + 
#  theme_bw() +
#  scale_fill_viridis_d() +
#  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + scale_y_log10()



######### NOW LOOK AT TOP 95% and further characterize ##########
#subset to 95% after density plot , see what the list is, n = 683 
dim(metabAUCs)
metabAUC95 <- metabAUCs %>% filter(.estimate>0.95)
dim(metabAUC95)#683

#merge with original file to find out how many are hilic pos, neg etc 
hilic  <- read.table("/home/rsm34/5asa/HMP2_metabolomics_edit.txt", header=T, sep="\t", check.names=FALSE) %>% select("Method","Metabolite") %>%
		rename(metab=Metabolite)
hilic$metab <- ifelse(substring(hilic$metab, 1,1) %in% c(0:9), paste0("X", hilic$metab), hilic$metab)
hilic$metab <- gsub(" ", ".", hilic$metab, fixed = TRUE)
metabAUC95 <- inner_join(metabAUC95, hilic,by="metab") 

#get most common metabolomic methods    
hilictable <- sort(table(metabAUC95$Method),decreasing=TRUE)[1:4]
hilictable
#  C18-neg HILIC-pos HILIC-neg    C8-pos
#      320       221       132        10

#find metabolites with large fold change and with zero values among users ## could be used in the text --> 45% negative. 
metabAUC95 <- metabAUC95 %>% rename(feature=metab)
AUClist<- metabAUC95 %>% select(feature)
AUClist2 <- rbind(c("SampleID"),c("binary154"),AUClist)
metab_683 <- mesal_metab %>% select(any_of(AUClist2$feature))
medianzero <- metab_683  %>% group_by(binary154) %>% summarise(across(X107.0364_0.71:X316.2282_13.53, median, na.rm= TRUE)) 
medianzero2 <- as.data.frame(t(medianzero)) %>% arrange(V1,V2)
colnames(medianzero2) <- c("nonuser","user")
medianzero2$nonuser <- as.numeric(medianzero2$nonuser)
medianzero2$user <- as.numeric(medianzero2$user)
(sum(medianzero2$nonuser==0)-1)/(nrow(medianzero2)-1) #45.2% (n-1 b/c binary154 is included, number of median=0) 
medianzero2 <- medianzero2 %>% arrange(user) 
medianzero2$foldchange <- medianzero2$user/(medianzero2$nonuser+1)
head(medianzero2,289)
tail(medianzero2,10)
medianzero2$feature <- rownames(medianzero2) 
medianzero2 <- medianzero2 %>% relocate(feature , .before = nonuser)
rownames(medianzero2) = NULL 
#which(medianzero2$feature=="X196.0609_2.81") #683 
medianzero2 <- medianzero2 %>% arrange(user) 
(sum(medianzero2$foldchange>1000))/(nrow(medianzero2)-1) #(number with fold change >1000) 
(sum(medianzero2$user>10000000))/(nrow(medianzero2)-1) #(relative abundance close >10000000) 

#################### STEP 4: CALCULATE THE DELTA M/Z to SEE WHAT BIOTRANSFORMATIONS THERE ARE ##################################################

#extract name
metabAUC95$mz <- gsub("X","",metabAUC95$feature)
metabAUC95$mz <- gsub("_.*","",metabAUC95$mz) 
metabAUC95$mz <- as.numeric(metabAUC95$mz) #this rounds automatically, which is annoying, but in some ways it is helpful, warning for nicotinuric acid  
metabAUC95$deltamz <- metabAUC95$mz - 154.0502 #parent compound
metabAUC95 <- metabAUC95 %>% arrange(deltamz)
metabAUC95$deltamz_round <- round(metabAUC95$deltamz,0)
metabAUC95$Method <- as.factor(metabAUC95$Method)

#write this table out to make a density plot (see below)
#write.table(metabAUC95,"metab_deltamz.txt",sep="\t",quote=FALSE,row.names=FALSE)
  


#### graph a few candidates for examples of biotransformations/host shifts #####
##### these were selected based on a lit search and/or annotations ##### 

candidates <- c("X154.0502_3.83","X316.0664_4.26","X224.0922_2.24","X196.0609_2.81","X210.0764_2.52","nicotinate","nicotinuric.acid","SampleID","X1.methylnicotinamide", "X373.1254_0.71","X242.0458_3")

metab_biotransf <- metab_rotated2 %>% select(any_of(candidates)) %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% relocate(binary154, .after = SampleID)

metab_biotransf2 <- inner_join(meta_data_IBD,metab_biotransf,by="SampleID") %>% drop_na(diagnosis) %>% select(1,31:41) %>% relocate(binary154, .after = SampleID) 

metab_biotransf3 <- gather(metab_biotransf2, metab, value, X154.0502_3.83:X242.0458_3)
metab_biotransf3$round <- substr(metab_biotransf3$metab,1,4)
metab_biotransf3$cohort <- "HMP2"

#test graph 
metab_biotransf3 %>%
 ggplot(aes(x=as.factor(binary154),y=value+1)) +
  geom_boxplot() + 
  theme_bw() +
  scale_fill_viridis_d() +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
  facet_wrap(~ metab) + scale_y_log10()

#write this table out to make side by side plots with PRISM (see biotransf_clean.R )
write.table(metab_biotransf3,"metab_biotransf3.txt",sep="\t",quote=FALSE,row.names=FALSE)
 
############################### STEP 5: CHECK FOR CORRELATIONS WITHIN METABOLITES, SEE HOW THEY CLUSTER, ARE THEY ALL RELATED? HOW MANY GROUPS ARE THERE? ##################

############################################################### correlation network ###############################################
#https://davetang.org/muse/2017/03/16/matrix-to-adjacency-list-in-r/
#https://www.r-graph-gallery.com/257-input-formats-for-network-charts.html
#http://mr.schochastics.net/netVizR.html check this out 
#https://kateto.net/networks-r-igraph

ibd154 <- metab_biotransf2 %>% select(SampleID,X154.0502_3.83)

metab_683_v2 <- inner_join(metab_683,ibd154,by="SampleID")  
metab_683_alone <- metab_683_v2 %>% filter(binary154==1) %>% select(-SampleID) %>% select(-binary154) 
metab_683_log <- log(metab_683_alone+1) 
my_cor_matrix <- cor(metab_683_log) 
my_cor_matrix[1:6, 1:6]
dim(my_cor_matrix)
my_cor_matrix[upper.tri(my_cor_matrix)] <- 42

my_cor_df <- melt(my_cor_matrix) %>% arrange(value)

my_cor_df <- filter(my_cor_df, value != 42) %>% filter(Var1 != Var2)
dim(my_cor_df)
head(my_cor_df)

summary(my_cor_df$value)

cor.test(metab_683_alone$X196.0609_2.81,metab_683_alone$X154.0502_3.83,method="spearman") # 0.5002581, P 4e-9


#this outputs a table of LCMS method for annotation in network analysis 
mergemethod <- metabAUC95 %>% select(feature,Method) %>% rename(Var1 = feature)
de<-data.frame("X154.0502_3.83","HILIC-pos")
names(de) <- c("Var1","Method") 
mergemethod2 <- rbind(mergemethod,de) 
write.table(mergemethod2  ,"mergemethod2.txt",sep="\t",quote=FALSE,row.names=FALSE) 

my_cor_df2 <- inner_join(my_cor_df,mergemethod2, by="Var1") 

#this one is just right 
my_adj_list85 <- my_cor_df2 %>% filter(abs(value) > 0.85)
names(my_adj_list85) <- c('from', 'to', 'weight','method')
dim(my_adj_list85)

#will output 0.85 version for graphing below 
#write.table(my_adj_list80 ,"my_adj_list80.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(my_adj_list90 ,"my_adj_list90.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(my_adj_list85 ,"my_adj_list85.txt",sep="\t",quote=FALSE,row.names=FALSE)
#graph in R studio, see below because iGraph is finicky in the computing cluster 

####### Grouping Based on Spearman Correlation Coefficients ############ --> similar to what elinav did in nature paper
###### find metabolite groups --> https://rdrr.io/cran/KMDA/man/spearman.group.html #####
library(KMDA)

metablist <- metabAUC95 %>% select(feature) 
grouping_list <- rbind(c("X154.0502_3.83"),metablist)
grouping_pre <- metab_rotated2 %>% select(any_of(grouping_list$feature))%>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% filter(binary154==1)%>% select(-binary154) 

cor(grouping_pre$X196.0609_2.81,grouping_pre$X154.0502_3.83,method="spearman") #checks correlation betwn drug and feature 

grouping_trans <- t(grouping_pre)
grouping_tran_mat <- as.matrix(grouping_trans)

groups <-spearman.group(grouping_tran_mat,0.85)
names <- row.names(grouping_trans)
grouped <- data.frame(names, groups) 
grouped %>% filter(groups==1) ## looks at a group 

#finds out where some of my friends are, all in different groups from parent 
grouped %>% filter(names=="X196.0609_2.81") #6
grouped %>% filter(names=="X316.0664_4.26") #461
grouped %>% filter(names=="X224.0922_2.24") #229
grouped %>% filter(names=="nicotinuric.acid") #148

#finds total number of groups
grouped %>% distinct(groups) #351

#finds singletons (subtract from total)
test<-grouped$groups 
test2 <- test[!(duplicated(test)|duplicated(test, fromLast=TRUE))]
length(test2) #266

#get means of each metabolite, to then select based on max in each cluster 
grouping_trans2 <- as.data.frame(grouping_trans)
grouping_trans2$mean <- (rowMeans(grouping_trans2))
grouping_trans2$names <- row.names(grouping_trans2)
row.names(grouping_trans2) = NULL 
groupedmeans <- grouping_trans2 %>% select(names, mean)

#now merge with hilic, means, and then select the max -- this is for halla (takes the max from each group for correlations)   
grouped3 <- inner_join(grouped, groupedmeans ,by="names") 
subset_groups <- grouped3 %>% group_by(groups) %>% slice_max(mean, n = 1) %>% arrange(mean)
write.table(subset_groups ,"groupedmetabnames_halla.txt",sep="\t",quote=FALSE,row.names=FALSE)


############################################# STEP 6: SEE HOW MUCH VARIANCE IS EXPLAINED BY MICROBES, HOST, DRUG ETC ########################################

###### make graphs between 154 and other metabolites among users #######

metab_5Ausers <- metab_rotated2 %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% filter(binary154==1) %>% select(-binary154) 

AUClist3 <- rbind(c("SampleID"),c("X154.0502_3.83"),c("binary154"),AUClist)
metab_users_95AUC <- metab_5Ausers %>% select(any_of(AUClist3$feature))
metab_users_95AUC[,2:684] <- log(metab_users_95AUC[,2:684]+1)
write.table(metab_users_95AUC ,"metab_halla_input_pre.txt",sep="\t",quote=FALSE,row.names=FALSE) #write table for halla --> see hallaprep_clean.R 


metab_regression <- inner_join(metab_users_95AUC,meta_data_IBD,by="SampleID") 

######read in the SPECIES file after pre-processing script made it (1632 --> n=1464)
species_rotated <- read.table("/home/rsm34/5asabb3/archive/taxa/metaphlan_rotated.txt",header=T,sep="\t")
species_rotated<-species_rotated[rowSums(species_rotated[,2:ncol(species_rotated)])!=0,,drop=FALSE] 	#there are a few rows that sum to zero (also a few that sum to <0.5

#subset to Ids above 
idlist <- metab_regression %>% select(SampleID) 
species_rotated2 <- inner_join(species_rotated,idlist,by="SampleID") 
row.names(species_rotated2) <- species_rotated2$SampleID
species_data <- species_rotated2 %>% select(matches("s__"))

species_data2 <- asin(sqrt(species_data*0.01))
data.bray=vegdist(species_data2) #Bray-Curtis distance
data.b.pcoa=cmdscale(data.bray,k=(nrow(species_data2)-1),eig=TRUE) #ordination

#cmdscale(data.bray,k=10,eig=TRUE)$GOF #ordination --> [1] 0.5854432 0.6361382 gives sense of variance explained 

str(data.b.pcoa)
species_pcoa = data.frame(PC1 = data.b.pcoa$points[,1], PC2 = data.b.pcoa$points[,2],PC3 = data.b.pcoa$points[,3],PC4 = data.b.pcoa$points[,4],PC5 = data.b.pcoa$points[,5],PC6 = data.b.pcoa$points[,6],PC7 = data.b.pcoa$points[,7],PC8 = data.b.pcoa$points[,8],
PC9 = data.b.pcoa$points[,9],PC10 = data.b.pcoa$points[,10])
species_pcoa$SampleID <- row.names(species_pcoa)
row.names(species_pcoa)=NULL 

###### add in diet data ######

#dietary data 
##############
#read in file 
dietarypre <- read.csv("/home/rsm34/5asa/hmp2_metadata_load.csv",stringsAsFactors=FALSE)%>% filter(data_type=="metabolomics")

#select columns of interest
dietary <- dietarypre %>% select(2,26:28,70:95,98,478) %>% rename(SampleID=External.ID) %>% select(1,7:27) 


dietary2 <- data.frame(lapply(dietary, function(x) {
               gsub("No, I did not consume these products in the last 7 days", "0", x)
        }))

dietary2 <- data.frame(lapply(dietary2, function(x) {
               gsub("Within the past 4 to 7 days", "0.18", x)
        }))

dietary2 <- data.frame(lapply(dietary2, function(x) {
               gsub("Within the past 2 to 3 days", "0.4", x)
        }))

dietary2 <- data.frame(lapply(dietary2, function(x) {
               gsub("Yesterday, 1 to 2 times", "0.66", x)
        }))

dietary2 <- data.frame(lapply(dietary2, function(x) {
               gsub("Yesterday, 3 or more times", "1.5", x)
        }))

dietary2 <- dietary2 %>% rename(soft_drinks=Soft.drinks..tea.or.coffee.with.sugar..corn.syrup..maple.syrup..cane.sugar..etc.,
                   diet_drinks=Diet.soft.drinks..tea.or.coffee.with.sugar..Stevia..Equal..Splenda.etc.,
                   fruit_juice=Fruit.juice..orange..apple..cranberry..prune.etc..,
                   alcohol=Alcohol..beer..brandy..spirits..hard.liquor..wine..aperitif..etc..,
		  yogurt=Yogurt.or.other.foods.containing.active.bacterial.cultures..kefir..sauerkraut.,
                   dairy=Dairy..milk..cream..ice.cream..cheese..cream.cheese.,
		  fruits=Fruits..no.juice...Apples..raisins..bananas..oranges..strawberries..blueberries,
		  veggies=Vegetables..salad..tomatoes..onions..greens..carrots..peppers..green.beans..etc.,
		  legumes=Beans..tofu..soy..soy.burgers..lentils..Mexican.beans..lima.beans.etc.,
		wholegrains=Whole.grains..wheat..oats..brown.rice..rye..quinoa..wheat.bread..wheat.pasta.,
		refgrains = Starch..white.rice..bread..pizza..potatoes..yams..cereals..pancakes..etc..,
		procmeat=Processed.meat..other.red.or.white.meat.such.as.lunch.meat..ham..salami..bologna,
		redmeat=Red.meat..beef..hamburger..pork..lamb.,
		whitemeat=White.meat..chicken..turkey..etc..,
		shellfish=Shellfish..shrimp..lobster..scallops..etc..,
		fish=Fish..fish.nuggets..breaded.fish..fish.cakes..salmon..tuna..etc..,
		sweets=Sweets..pies..jam..chocolate..cake..cookies..etc..) %>% select(-c(9,12)) 

dietary2[,c(2:20)] <- as.data.frame(sapply(dietary2[,c(2:20)] ,as.numeric)) #convert to numeric 

dietary3 <- dietary2[,c(2:20)]

dietary3 <- dietary3 %>% mutate_all(~replace(., is.na(.), 0)) #make sure NA is = 0 for diet data. 
rownames(dietary3) <- dietary2$SampleID

#now will try to make dietary data into two patterns  
#https://data.library.virginia.edu/getting-started-with-factor-analysis/
fa <- factanal( ~., 2, data = dietary3 ,rotation = "varimax")

prudent <- fa$loadings[,1] %*% t(dietary3 )
prudent <- t((prudent))
colnames(prudent)<-"prudscore"
prudent <- as.data.frame(prudent)
prudent$SampleID <- rownames(prudent)

prudent$prudscore[prudent$prudscore==0] <- median(prudent$prudscore) ##this sets missing to median value 

########### merge ####### 
metab_regression2<- inner_join(metab_regression,species_pcoa, by="SampleID") 
metab_regression3<- inner_join(metab_regression2,prudent, by="SampleID") 
#delete?
#write.table(metab_regression3 ,"metab_regression3.txt",sep="\t",quote=FALSE,row.names=FALSE)#write table to run this as a batch file 

#### test 1 - n ac 5-asa ######
library(relaimpo) 

#first generate linear model to get R2 which in this example is 0.50 
linmod3_5 <- lm(X196.0609_2.81 ~ X154.0502_3.83 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + diagnosis + consent_age+Antibiotics+prudscore,data=metab_regression3)

#now we estimate the groups contribution (drug (0.60), microbiome (0.27), host (0.12))
metrics <- calc.relimp(linmod3_5,type=("lmg"))
aaa<- data.frame(metrics@lmg)
sum(aaa[2:11,])/sum(aaa) #microbiome contribution 
sum(aaa[1,])/sum(aaa) #drug contribution 
sum(aaa[12:16,])/sum(aaa) #host contribution 


### test 2 - nicotinuric acid ####### --> no host contribution 
linmod4 <- lm(nicotinuric.acid ~ X154.0502_3.83 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + diagnosis + consent_age+Antibiotics,data=metab_regression2) 
metrics2 <- calc.relimp(linmod4,type=("lmg"))

bbb<- data.frame(metrics2@lmg)
sum(bbb[2:11,])/sum(bbb) #microbiome contributio
sum(bbb[1,])/sum(bbb) #drug contributio
sum(bbb[12:16,])/sum(bbb) #host contribution --> 0% 


linmod5 <- lm(X288.0487_2.25 ~ X154.0502_3.83 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + diagnosis + consent_age+Antibiotics,data=metab_regression2) 
metrics4 <- calc.relimp(linmod5,type=("lmg"))
ccc<- data.frame(metrics4@lmg)
sum(ccc[2:11,])/sum(ccc) #microbiome contributi
sum(ccc[1,])/sum(ccc) #drug contributi


#### now run glms for each of the metabolites, and then check the proportion of variance explained according to drug, microbiome, and host factors (dx, age, sex, diet, abx) --- takes a long time (hours) ####
	### see relaimpo.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(relaimpo) 

metab_regression3 <- read.table("/home/rsm34/5asabb3/intermediate_files/metab_regression3.txt",header=T,sep="\t")

kmetabs_glm = NULL
kmetabs_glm= as.data.frame(matrix(-9,690,5))
colnames(kmetabs_glm) = c("feature","totalR2","drugR2","microbiomeR2","hostR2")

for(i in c(3:684)) #### ran for a subset on 4.7
	{	
  		print(paste("Running entity:", i, "now"))

		b<-lm(metab_regression3[,i]~metab_regression3$X154.0502_3.83 + metab_regression3$PC1 + metab_regression3$PC2+ metab_regression3$PC3+ metab_regression3$PC4+ metab_regression3$PC5+ metab_regression3$PC6+ metab_regression3$PC7+ metab_regression3$PC8+ metab_regression3$PC9+ metab_regression3$PC10+
                      metab_regression3$sex +metab_regression3$diagnosis +metab_regression3$consent_age+metab_regression3$Antibiotics+metab_regression3$prudscore)
 		kmetabs_glm[i,1]= colnames(metab_regression3)[i]
 		kmetabs_glm[i,2]=summary(b)$r.squared
		metrics <- calc.relimp(b,type=("lmg"))
		aaa<- data.frame(metrics@lmg)
 		kmetabs_glm[i,3]=sum(aaa[1,])/sum(aaa)
 		kmetabs_glm[i,4]=sum(aaa[2:11,])/sum(aaa)
 		kmetabs_glm[i,5]=sum(aaa[12:16,])/sum(aaa)
	}
	kmetabs_glm <- kmetabs_glm %>% arrange(drugR2)
write.table(kmetabs_glm ,"kmetabs_relaimpo.txt",sep="\t",quote=FALSE,row.names=FALSE)




###################################################################################################### PLOTTING #################################################################################################

### in Rstudio 


###################################### Supp Fig 4) NETWORK ANALYSIS ##########################################

##### plot these using iGraph   #####

my_adj_list85 <- read.table("my_adj_list85.txt",header=T,sep="\t")

mergemethod <- read.table("mergemethod2.txt",header=T,sep="\t")

library(igraph)
library(RColorBrewer)
library(viridis)
library(dplyr)

 
# create igraph S3 object
links <-my_adj_list85[,1:3]

idlist <- my_adj_list85%>%  
  select(from, to) %>% 
  t %>% c %>% unique
idlist <- as.data.frame(idlist) 
names(idlist) <- c("Var1")

nodes<-inner_join(idlist,mergemethod) 

nodes$Method <- as.character(nodes$Method)

nodes[which(nodes$Var1=="X196.0609_2.81"), 2] <- "N_Ac_5ASA"
nodes[which(nodes$Var1=="X154.0502_3.83"), 2] <- "mesalamine"


net2 <- graph_from_data_frame(d=links, vertices=nodes, directed=F)

# Create a vector of color

coul <- brewer.pal(6, "Set1") 
coul2 <- c( "#DC0000FF","#3C5488FF","#00A087FF","#7E6148FF","#FF7F00", "#FFFF33")
coul3 <-(viridis_pal(option="B")(6))
my_color <- coul2[as.numeric(as.factor(V(net2)$Method))]

# store original margins
#orig_mar <- par()$mar 
 
# set new margins to limit whitespace in plot
par(mar=rep(.1, 4))


node.size<-replicate(383, 5)
node.size2<-setNames(node.size,nodes$Var1)
qq<-which(names(node.size2) == "X154.0502_3.83")
node.size2[qq] <- 15

set.seed(1)

plot(net2,  edge.width = E(net2)$weight,vertex.size=node.size2,vertex.color=my_color,layout=layout.fruchterman.reingold, vertex.label = ifelse(degree(net2) > 88, "center", NA),vertex.label.color="black",vertex.label.font=c(2),vertex.label.cex=c(1.5))
legend("bottomright", 
       legend=paste( levels(as.factor(V(net2)$Method)), " ", sep=""), 
       col = coul2 , 
       bty = "n", pch=20 , pt.cex = 2, cex = 1,
       text.col="black" , horiz = F)

#7x6 pdf, edited a little in inkscape to bring orange node to front #

############################ SUPP FIG 3A:  CONTRIBUTION TO VARIANCE #####################################
library(tidyr)
library(ggplot2)
library(dplyr)

kmetabs_glm <- read.table("kmetabs_relaimpo.txt",header=T,sep="\t")




##### OLD########

kmetabs_glm <- kmetabs_glm %>% filter(totalR2 >-2)	%>% select(-2)

kmetabs_gather <- kmetabs_glm %>% gather(Class, Value, -feature) 

# sum the abundance for each class, across all IDs, & sort the result
sort.class <- kmetabs_gather %>% 
  count(Class, wt = Value) %>%
  arrange(desc(n)) %>%
  pull(Class)

feature.order <- kmetabs_gather%>%
  filter(Class == sort.class[1]) %>%
  arrange(desc(Value)) %>%
  pull(feature)

pp<-kmetabs_gather %>%
  mutate(feature = factor(feature, levels = feature.order)) %>%
  mutate(Class = factor(Class, levels = rev(sort.class))) %>%
  ggplot(aes(x = feature, y = Value, fill = Class)) +
  geom_col(width = 1) #geom_col is equivalent to geom_bar(stat = "identity")
pp + ylab("% of variance explained") + scale_fill_manual(values = c("#3C5488FF","#00A087FF","#E64B35FF"),labels = c("Host", "Microbiome", "5-ASA"),name = "Variable")+ xlab("Metabolite") +
 theme( axis.text.x = element_blank(), axis.ticks = element_blank())

#5.5 x 3 

#### NEW 3.16.22 --> has unexplained variance #####

kmetabs_glm <- kmetabs_glm %>% mutate(drugR2_adj = totalR2 * drugR2) %>% mutate(microbiomeR2_adj = totalR2 * microbiomeR2)%>% mutate(hostR2_adj = totalR2 * hostR2) %>% mutate(unexplained = (1-totalR2)) %>% filter(totalR2 >-8) 

kmetabs_glm <- kmetabs_glm %>% select(-c(2,3,4,5))

kmetabs_gather <- kmetabs_glm %>% gather(Class, Value, -feature) 

# sum the abundance for each class, across all IDs, & sort the result
sort.class <- kmetabs_gather %>% 
  count(Class, wt = Value) %>%
  arrange(desc(n)) %>%
  pull(Class)

feature.order <- kmetabs_gather%>%
  filter(Class == sort.class[1]) %>%
  arrange((Value)) %>%
  pull(feature)

pp<-kmetabs_gather %>%
  mutate(feature = factor(feature, levels = feature.order)) %>%
  mutate(Class = factor(Class, levels = rev(sort.class))) %>%
  ggplot(aes(x = feature, y = Value, fill = Class)) +
  geom_col(width = 1) #geom_col is equivalent to geom_bar(stat = "identity")
pp + ylab("% of variance explained") + scale_fill_manual(values = c("#3C5488FF","#00A087FF","#E64B35FF","grey"),labels = c("Host", "Microbiome", "5-ASA","Unexplained"),name = "Variable")+ xlab("Metabolite") +
 theme( axis.text.x = element_blank(), axis.ticks = element_blank())

#5.5 x 3 

############ BONUS STUFF (UNUSED) ########

#time charts

myanim<-p+   geom_point() +   transition_reveal(visit_num)

library("magick")

animate(myanim, duration = 10, fps=10, renderer = magick_renderer())

anim_save("metab6010smooth.gif", animation = last_animation())


# SUPP: delta M/z 
###############
metabmz <- read.table("metab_deltamz.txt",header=T,sep="\t")
ggplot(metabmz, aes(x=deltamz_round)) + geom_density(alpha=0.3,fill="yellow") +theme_cowplot() +xlab("Delta (m/z)") + theme(text = element_text(size=20))


