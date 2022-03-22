############################################################### CODE ###################################################

library(dplyr)
library(tidyr)
library(nlme)
library(ggplot2)
library(stringr)
library(caret)


###### metadata ###################
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
dim(meta_data)

#subset to IBD users only 

IBDids <- meta_data %>% filter(diagnosis=="CD"|diagnosis=="UC") %>% select(SampleID) 

####### MTX uniref data ############

rnauni_rotated_metab <- read.table("/home/rsm34/5asabb3/archive/RNA_uniref/rnauni_subset_bb3.txt",header=T,sep="\t") %>% mutate(binary196 = if_else(X196.0609_2.81 > 10000000,1,0)) %>% relocate(binary196, .after = UNMAPPED)
dim(rnauni_rotated_metab)#213x 56546 

##### sen/spec function ######
##########################
rnauni_rotated_metab2 <- inner_join(IBDids,rnauni_rotated_metab)
mtx_uniref_pa <- rnauni_rotated_metab2[,5:ncol(rnauni_rotated_metab2)] 
mtx_uniref_pa[mtx_uniref_pa>0] <-1

df <- mtx_uniref_pa

kmetabs2 = NULL
kmetabs2 = as.data.frame(matrix(-9,ncol(df),3))
colnames(kmetabs2) = c("metab","sens","spec")

for(i in c(1:ncol(df))){	
		kmetabs2[i,1]= colnames(df)[i]
 		kmetabs2[i,2]= (sum(df[,i]*df[,c("binary196")])/sum(df[,c("binary196")]))
		kmetabs2[i,3]= (sum(df[,i]+df[,c("binary196")]==0)/sum(df[,c("binary196")]==0))
	}
	kmetabs2$sum <- (kmetabs2$sens + kmetabs2$spec)
	kmetabs2 <- kmetabs2 %>% arrange(sum)

write.table(kmetabs2, "sen_spec_196_ibdmtxbb3_all.txt", sep="\t", quote=FALSE,row.names=FALSE) #### --> use this for base plotting purposes 

### verify with caret package, makes sense ###
df$UniRef90_D4KA95 <- factor(df$UniRef90_D4KA95, levels=rev(levels(df$UniRef90_D4KA95)))
df$binary196 <- factor(df$binary196,levels=rev(levels(df$binary196)))
xtab <- table(df$UniRef90_D4KA95, df$binary196)
sensitivity(df$UniRef90_D4KA95, df$binary196) 

#filter 
senspec196_2 <- kmetabs2 %>% filter(sens > 0.5) %>% filter(spec>0.5) %>% rename(feature=metab)
dim(senspec196_2) #2015, some of which could be promising, reviewed on uniprot 
senspec196_2$uniprot <- str_replace(senspec196_2$feature,"UniRef90_","")

write.table(senspec196_2, "sen_spec_196_ibdmtxbb3_select.txt", sep="\t", quote=FALSE,row.names=FALSE)


######## sensitivity spec plotting #######
library(ggrepel) 

#rest in Rstudio with txt files above
senspec196 <- read.table("/home/rsm34/5asabb3/submission/keyoutput/sen_spec_196_ibdmtxbb3_all.txt",header=T, sep="\t") %>% filter(metab!="binary196")

sen_spec_names <- c("UniRef90_R6TIX3","UniRef90_C7H1G6","UniRef90_A0A1C7H380","UniRef90_A0A076IMG8","UniRef90_A0A076IRG9","UniRef90_A0A076IXU4","UniRef90_C7H1G6","UniRef90_I9R3P9","UniRef90_A0A0P0M2W7","UniRef90_U2QX51","UniRef90_R6CZ24","UniRef90_R5CY66","UniRef90_T5S060")
senspecmap <- senspec196[senspec196$metab %in% sen_spec_names,]

senspecmap$labelsraaj <- str_replace(senspecmap$metab,"UniRef90_","")


p3 <- ggplot(senspec196, aes(1-spec,sens)) + geom_point(alpha=0.4, colour = "#8491B4FF",size=2,shape=16) + theme_bw(base_size=14) + geom_point(data = senspecmap , aes(1-spec,sens), fill = "red",shape=21,colour="white",stroke=1) + geom_segment(aes(x = 0, y = 0.5, xend = 0.5, yend = 0.5),linetype="dashed") + 
  geom_segment(aes(x = 0.5, y = 0.5, xend = 0.5, yend = 1.0),linetype="dashed") +
	  xlab("1-Specificity (TN/N)") + ylab("Sensitivity (TP/P)")

p3 + geom_segment(aes(x = 0, y = 0.5, xend = 0.5, yend = 0.5),linetype="dashed")


fig3b<-p3 + geom_label_repel(data = senspecmap ,aes(label=labelsraaj),box.padding = 0.1) #5X5 PDF 

ggsave(filename='/home/rsm34/5asabb3/submission/figures/fig3b.png', plot=fig3b , width = 5, height = 5, dpi = 600) 

 





