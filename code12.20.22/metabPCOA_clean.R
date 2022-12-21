######################################################################################################
# R programs for creating "metabolite" PCOA					#
# ============================================					#
# Purposes:  ordination
######################################################################################################

#references: 
#https://bitbucket.org/biobakery/physalia_2018/wiki/Lab%20%237:%20Metagenomic%20Visualization

############################################################### CODE ###################################################

library(dplyr)
library(tidyr)
library(nlme)
library(vegan) 
library(ggplot2)
library(cowplot) 

#prepare metadata 
meta_data <- read.table("/home/rsm34/5asabb3/archive/metadata/meta_data.txt",header=T,sep="\t")
cols <- c("oral5asa_lc","diagnosis","biologics","site_name","Antibiotics","bowelsurg","dysbiosis","Participant.ID","steroids","pr5asa","bond5asa","immunomod","sex")
meta_data[cols] <- lapply(meta_data[cols], factor) #converts these to a factor variable)

meta_data$dysbiosis <- relevel(meta_data$dysbiosis, ref = "FALSE") # Set 'inactive' as reference
meta_data$diagnosis <- relevel(meta_data$diagnosis, ref = "nonIBD") # Set nonIBD as reference
meta_data <- meta_data %>% mutate(IBDusers = if_else(oral5asa_lc == 1,2,if_else(oral5asa_lc==0 & as.numeric(diagnosis) == 3,1, if_else(oral5asa_lc==0 & as.numeric(diagnosis) == 2,1,0)))) 
meta_data$IBDusers[meta_data$IBDusers=="2"] <- "User"
meta_data$IBDusers[meta_data$IBDusers=="1"] <- "IBDNon-User"
meta_data$IBDusers[meta_data$IBDusers=="0"] <- "nonIBD"
meta_data$IBDusers <- as.factor(meta_data$IBDusers)
meta_data$IBDusers <- relevel(meta_data$IBDusers, ref = "nonIBD") # Set nonIBD as reference


meta_data_IBD <- meta_data %>% filter(diagnosis=="CD"|diagnosis=="UC")
meta_data_IBD$diagnosis <- droplevels(meta_data_IBD$diagnosis)
 

#read in METABOLOMICS file after pre-processing script made it 
metab_rotated2 <- read.table("/home/rsm34/5asabb3/archive/metabolite/metab_data_nozero2.txt",header=T,sep="\t")
metab_IBD<- inner_join(meta_data_IBD,metab_rotated2,by="SampleID") %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% relocate(binary154, .after = IBDusers)

metab_data_logt = metab_IBD
dim(metab_data_logt) #283x81898

metab_data_logt[,32:ncol(metab_data_logt)] <- log((metab_data_logt[,32:ncol(metab_data_logt)] + 1))
yy <-round(.1*nrow(metab_data_logt))
mesal_metab_filt <- metab_data_logt[,32:ncol(metab_data_logt)]
mesal_metab_filt <- mesal_metab_filt[,colSums(mesal_metab_filt != 0) >= yy] #this filters features with more than 10% missing. 
meta_data_metab <-metab_data_logt[,c(1:31)] #get ready to bind 
mesal_metab_filt_comb <- cbind(meta_data_metab,mesal_metab_filt)

metab_IBD_nlme <- mesal_metab_filt_comb 
dim(metab_IBD_nlme)
#dim 283 x 74330 

########## now look at differences between UC/5ASA+ and UC/5ASA-; look at CD/5ASA+ vs UC/5ASA+ : PCOA  ##############################

metab_bray <- mesal_metab_filt_comb
metab_bray <- metab_bray %>% mutate(mesalID = if_else(binary154==0 & as.numeric(diagnosis)==1,1,if_else(binary154==1 & as.numeric(diagnosis) == 1,2, if_else(binary154==0 & as.numeric(diagnosis) == 2,103,if_else(binary154==1 & as.numeric(diagnosis) == 2,104,0))))) 
metab_bray$mesalID2 <- paste0(metab_bray$SampleID,"_",metab_bray$mesalID)
metab_bray <- metab_bray %>% relocate(mesalID2 , .after = IBDusers) %>% select(-mesalID)

row.names(metab_bray) <- metab_bray$mesalID2
check2 <- vegdist(metab_bray[,34:ncol(metab_bray)], method="bray")

library(tidyr) 
data.b.pcoa=cmdscale(check2,k=(nrow(metab_bray)-1),eig=TRUE) #ordination on bray curtis log transformed metabolites 
str(data.b.pcoa)
pcoa = data.frame(PC1 = data.b.pcoa$points[,1], PC2 = data.b.pcoa$points[,2])
pcoa$diagnosis = metab_bray$mesalID2

#for axes labels
cmdscale(check2,k=1,eig=TRUE)$GOF #0.1047772 
cmdscale(check2,k=2,eig=TRUE)$GOF #0.1867013 - 0.1047772 = 0.0819241


#write.table(pcoa, "pcoa5asametab.txt", sep="\t", quote=FALSE,row.names=FALSE)


####### adonis ########

adonis(check2  ~ metab_bray$diagnosis) #2.2%, p <0.001 
adonis(check2  ~ as.factor(metab_bray$binary154)) #6.8%, p<0.001

metab_bray_filt <- metab_bray %>% drop_na(dysbiosis)
check3 <- vegdist(metab_bray_filt[,34:ncol(metab_bray_filt)], method="bray")

adonis(check3  ~ metab_bray_filt$dysbiosis) #2.6%, p<0.001 


################################################################################# PLOTTING ########################################################################################################
###in rstudio
#pcoa <- read.table("pcoa5asametab.txt",header=T,sep="\t")
pcoa$diagnosis2 <- gsub(".*_","",pcoa$diagnosis) #assigns to groups 
 
pcoa <- pcoa %>% mutate(mesalstatus = if_else(diagnosis2==104,1,if_else(diagnosis2==2,1,0)))
 
pcoa$diagnosis2[pcoa$diagnosis2=="1"] <- "CD 5-ASA(-)"
pcoa$diagnosis2[pcoa$diagnosis2=="2"] <- "CD 5-ASA(+)"
pcoa$diagnosis2[pcoa$diagnosis2=="103"] <- "UC 5-ASA(-)"
pcoa$diagnosis2[pcoa$diagnosis2=="104"] <- "UC 5-ASA(+)"

#5-ASA plot
my_col_med = c('CD 5-ASA(-)'='#3C5488FF','UC 5-ASA(-)'='#8491B4FF','UC 5-ASA(+)'='#F39B7FFF','CD 5-ASA(+)'='#DC0000FF')
p1 = ggplot(pcoa, aes(x=PC1, y=PC2,fill = diagnosis2)) + geom_point(shape=21,size=2.5,color="black") + scale_fill_manual(values=my_col_med)   + theme_cowplot(12) + labs(fill="5-ASA use by Dx") + ylab("PCoA2, 8.2%") + xlab("PCoA1, 10.4%") 
#p1

#add ellipses --> pdf 6x5
p3 = ggplot(pcoa,aes(x=PC1, y=PC2)) + geom_point(shape=21,aes(fill=diagnosis2),size=2.5,color="black") + scale_fill_manual(values=my_col_med) + 
  stat_ellipse(aes(x=PC1, y=PC2,color=as.factor(mesalstatus)),type = "norm",linetype=2) + 
   theme_cowplot(16) + labs(fill="5-ASA use by Dx") + ylab("PCoA2, 8.2%") + xlab("PCoA1, 10.4%") + scale_colour_manual(values=c("black","black"),guide=FALSE) 
#p3

ggsave(filename='/home/rsm34/5asabb3/submission/figures/fig2a.png', plot=p3, width = 7.5, height = 6, dpi = 600) 
