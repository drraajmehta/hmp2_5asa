
############################################################### CODE (re-reviewed 1/22/21) ###################################################

library(dplyr)
library(nlme)
library(stringr)
library(tidyr)
library(ggplot2)


##### STEP 1: HAVE TO GET DNA AND RNA FILES WITH THE SAME PARTICIPANTS AND WITH THE SAME FEATURES (IDENTICAL TABLE FORMAT) SO THAT I CAN ADJUST FOR DNA COPY ################################

#check overlap between rna and dna 
###################################
#$ cat /home/rsm34/5asabb3/archive/RNA_uniref/MTXunirefproptable2.txt | cut -f1 > mtxnamessubset.txt
#$ cat /home/rsm34/5asabb3/archive/DNA_uniref/DNAunirefproptable2.txt | cut -f1 > dnanamessubset.txt  #they are labeled subset since they are based off filtered prop tables,  which don't include all features  

mtxnames <- read.table("/home/rsm34/5asabb3/archive/mtxnamessubset.txt",header=T,sep="\t")
dnanames <- read.table("/home/rsm34/5asabb3/archive/dnanamessubset.txt",header=T,sep="\t") #file made with the above in unix

dim(mtxnames) (far fewer mtx with bb3, reviewed with Eric and Curtis) 
#[1] 61256
dim(dnanames)
#[1] 379036

matchingnames <- inner_join(mtxnames,dnanames,by="ID")
dim(matchingnames) 
#56542

#write.table(matchingnames, "matchedMTXDNAuninames_bb3.txt", sep="\t", quote=FALSE,row.names=FALSE,col.names=FALSE)

############## build DNA file ##########################
chunk1 <- read.table("/home/rsm34/5asabb3/archive/DNA_uniref/DNAunirefproptable2.txt", header=T, sep="\t")
dnauniselect <- inner_join(matchingnames, chunk1)

#transpose and sort again 
dim(dnauniselect) #56542x 214
dnauniselect_rotated = setNames(data.frame(t(dnauniselect[,-1])), dnauniselect[,1])
dnauniselect_rotated<- cbind (SampleID=row.names(dnauniselect_rotated),dnauniselect_rotated)
row.names(dnauniselect_rotated) <- NULL
dim(dnauniselect_rotated) #213*56543

############## build RNA file ########################## (start with output from uniref_proptable_mtx.R)
chunk1 <- read.table("/home/rsm34/5asabb3/archive/RNA_uniref/MTXunirefproptable2.txt", header=T, sep="\t")
rnauniselect <- inner_join(matchingnames, chunk1)

dim(rnauniselect) #56542x214 
rnauniselect_rotated = setNames(data.frame(t(rnauniselect[,-1])), rnauniselect[,1])
rnauniselect_rotated<- cbind (SampleID=row.names(rnauniselect_rotated),rnauniselect_rotated)
row.names(rnauniselect_rotated) <- NULL
dim(rnauniselect_rotated) #213x56543 

######## now clean to make sure they fit with metab file ###########
#read in binary154 made elsewhere 
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt",header=T,sep="\t") 

#merge with metabs file to see which participants get selected 
metab_ids <- metab_select[,c(1,2,5)]
dnauni_ids <- dnauniselect_rotated %>% select(SampleID)
rnauni_ids <- rnauniselect_rotated %>% select(SampleID)
gene_ids <- inner_join(dnauni_ids,rnauni_ids,by="SampleID")
gene_metabids <- inner_join(metab_ids,gene_ids,by="SampleID") 

#> dim(gene_metabids)
#[1] 213   3

dnauni_rotated_metab<- inner_join(gene_metabids,dnauniselect_rotated,by="SampleID") 
rnauni_rotated_metab<- inner_join(gene_metabids,rnauniselect_rotated,by="SampleID") 


#sort column variables 
dnauni_rotated_metab <- dnauni_rotated_metab %>% select("SampleID","binary154","X196.0609_2.81","UNMAPPED",sort(colnames(.)))
rnauni_rotated_metab <- rnauni_rotated_metab %>% select("SampleID","binary154","X196.0609_2.81","UNMAPPED",sort(colnames(.)))

#transform these 
rnauni_rotated_metab[,5:ncol(rnauni_rotated_metab)] <- asin(sqrt(rnauni_rotated_metab[,5:ncol(rnauni_rotated_metab)]))
dnauni_rotated_metab[,5:ncol(dnauni_rotated_metab)] <- asin(sqrt(dnauni_rotated_metab[,5:ncol(dnauni_rotated_metab)]))

write.table(dnauni_rotated_metab, file="dnauni_subset_bb3.txt", sep="\t", quote=FALSE,row.names=FALSE)
write.table(rnauni_rotated_metab, file="rnauni_subset_bb3.txt", sep="\t", quote=FALSE,row.names=FALSE)

### can start here if running interactively 1/22 checked ### 
#dnauni_rotated_metab <- read.table("/home/rsm34/5asabb3/archive/DNA_uniref/dnauni_subset_bb3.txt",header=T,sep="\t")
#rnauni_rotated_metab <- read.table("/home/rsm34/5asabb3/archive/RNA_uniref/rnauni_subset_bb3.txt",header=T,sep="\t")

dnauni_rotated_metab <- dnauni_rotated_metab %>% mutate(binary196 = if_else(X196.0609_2.81 > 10000000,1,0)) 
rnauni_rotated_metab <- rnauni_rotated_metab %>% mutate(binary196 = if_else(X196.0609_2.81 > 10000000,1,0))

dnauni_rotated_metab <- dnauni_rotated_metab %>% relocate(binary196, .after = binary154)
rnauni_rotated_metab <- rnauni_rotated_metab %>% relocate(binary196, .after = binary154)

dim(rnauni_rotated_metab)
dim(dnauni_rotated_metab)


#213x 56546 for both 



#build metadata file
#####################
meta_data <- read.table("/home/rsm34/5asa/meta_data.txt",header=T,sep="\t")
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

dnauni_ids <- dnauni_rotated_metab[,1:3]

meta_data_uni <- inner_join(dnauni_ids,meta_data)

dim(meta_data_uni)

##### nlme function ######
##########################

#subset to IBD users only --> this was already taken into account when i pulled the bb3 files

IBDids <- meta_data_uni %>% filter(diagnosis=="CD"|diagnosis=="UC") %>% select(SampleID) 
dnauni_rotated_metab2 <- inner_join(IBDids,dnauni_rotated_metab)
rnauni_rotated_metab2 <- inner_join(IBDids,rnauni_rotated_metab)
meta_data_uni2 <- inner_join(IBDids,meta_data_uni)

dim(dnauni_rotated_metab2)
dim(rnauni_rotated_metab2)
dim(meta_data_uni2) 


	feat.names <- colnames(dnauni_rotated_metab2)[6:ncol(dnauni_rotated_metab2)]
	no.feat <- length(feat.names)


		#create a file to hold the coefficients
		diag_sig = NULL
		diag_sig = as.data.frame(matrix(-9,1,3))
		colnames(diag_sig) = c("feature","beta","p")

		# loop over feat names
		for(i in feat.names){ 

  			# print status
  			print(paste("Running entity:", i, "which is", which(feat.names==i), "out of", no.feat))

  			# create temporary data matrix and model formula
			a1 = which(colnames(dnauni_rotated_metab2)==i)
			a2 = which(colnames(rnauni_rotated_metab2)==i)
			a3 = as.data.frame(dnauni_rotated_metab2[,a1])
			colnames(a3) = "dna"
			a4 = as.data.frame(rnauni_rotated_metab2[,a2])
			colnames(a4) = "rna"
			a5 = as.data.frame(cbind(a3,a4))

  			tmp <- meta_data_uni2[, c("binary154","Antibiotics","consent_age","dysbiosis","Participant.ID","diagnosis")]
			tmp2 <- cbind(a5,tmp) 
  			fml <- as.formula( paste( "rna", "~", paste(c("binary154","Antibiotics","consent_age","dysbiosis","diagnosis","dna"), collapse="+") ) )

  			# assign fit to list by name
  			b <- lme(fml, random=~1|Participant.ID, data=tmp2,na.action=na.omit,control=lmeControl(returnObject=TRUE))
 			diag_sig[i,1]= feat.names[i]
 			diag_sig[i,2]=(b$coefficients)$fixed[2]
 			diag_sig[i,3]=summary(b)$tTable[2,5]
			diag_sig$names=rownames(diag_sig)
		}
	diag_sig$FDR <- p.adjust(diag_sig$p,method="fdr")
	diag_sig <- diag_sig %>% arrange((p)) 
	write.table(diag_sig, file="dnarna_ibd_nlmebb3.txt", sep="\t", quote=FALSE,row.names=FALSE)


####getting a list of names to upload to uniprot ######

dnarna_nlmebb3 <- read.table("dnarna_ibd_nlmebb3.txt",header=T,sep="\t")
dnarna_nlmebb3$uniprot <- str_replace(dnarna_nlmebb3$feature,"UniRef90_","")
write.table(dnarna_nlmebb3, file="dnarna_uniprotbb3.txt", sep="\t", quote=FALSE,row.names=FALSE)

############################################################## make a few plots based on Uniprot search #############################################################################################################

#candidates
#UniRef90_R6TIX3 - Acetyl-CoA acetyltransferase
#UniRef90_C7H1G6 - GNAT domain 

#80% homology
#UniRef90_R6CZ24
#UniRef90_R5CY66
#UniRef90_T5S060
#UniRef90_U2QX51

###### density plots #######

######################
#uniref density #
#####################
unirefs_plotting <- rnauni_rotated_metab2 %>% select(matches("UniRef90_R6TIX3|UniRef90_C7H1G6|binary154|UniRef90_R6CZ24|UniRef90_R5CY66|UniRef90_T5S060|UniRef90_U2QX51|UniRef90_A0A3E5DVA8"))

u1 <- unirefs_plotting %>% select(matches("UniRef90_R6TIX3|binary154")) %>% filter(UniRef90_R6TIX3 >0) %>% mutate(median=median(UniRef90_R6TIX3)) %>% mutate(scaled = log10(UniRef90_R6TIX3/median)) %>% mutate(feature_name = "UniRef90_R6TIX3") %>% rename(abundance=UniRef90_R6TIX3)
u2 <- unirefs_plotting %>% select(matches("UniRef90_C7H1G6|binary154")) %>% filter(UniRef90_C7H1G6 >0) %>% mutate(median=median(UniRef90_C7H1G6)) %>% mutate(scaled = log10(UniRef90_C7H1G6/median)) %>% mutate(feature_name = "UniRef90_C7H1G6") %>% rename(abundance=UniRef90_C7H1G6)
u3 <- unirefs_plotting %>% select(matches("UniRef90_R6CZ24|binary154")) %>% filter(UniRef90_R6CZ24 >0) %>% mutate(median=median(UniRef90_R6CZ24)) %>% mutate(scaled = log10(UniRef90_R6CZ24/median)) %>% mutate(feature_name = "UniRef90_R6CZ24") %>% rename(abundance=UniRef90_R6CZ24)
u4 <- unirefs_plotting %>% select(matches("UniRef90_R5CY66|binary154")) %>% filter(UniRef90_R5CY66 >0) %>% mutate(median=median(UniRef90_R5CY66)) %>% mutate(scaled = log10(UniRef90_R5CY66/median)) %>% mutate(feature_name = "UniRef90_R5CY66") %>% rename(abundance=UniRef90_R5CY66)
u5 <- unirefs_plotting %>% select(matches("UniRef90_T5S060|binary154")) %>% filter(UniRef90_T5S060 >0) %>% mutate(median=median(UniRef90_T5S060)) %>% mutate(scaled = log10(UniRef90_T5S060/median)) %>% mutate(feature_name = "UniRef90_T5S060") %>% rename(abundance=UniRef90_T5S060)
u6 <- unirefs_plotting %>% select(matches("UniRef90_U2QX51|binary154")) %>% filter(UniRef90_U2QX51>0) %>% mutate(median=median(UniRef90_U2QX51)) %>% mutate(scaled = log10(UniRef90_U2QX51/median)) %>% mutate(feature_name = "UniRef90_U2QX51") %>% rename(abundance=UniRef90_U2QX51)

uni_density_outfile <- do.call("rbind",list(u1,u2,u3,u4,u5,u6))
write.table(uni_density_outfile,"/home/rsm34/5asabb3/reviews/keyfiles/uni_dens_outfile_rev.txt",sep="\t",quote=FALSE,row.names=FALSE)

rev_dens <- unirefs_plotting %>% select(matches("UniRef90_A0A3E5DVA8|binary154")) %>% filter(UniRef90_A0A3E5DVA8>0) %>% mutate(median=median(UniRef90_A0A3E5DVA8)) %>% mutate(scaled = log10(UniRef90_A0A3E5DVA8/median)) %>% mutate(feature_name = "UniRef90_A0A3E5DVA8") %>% rename(abundance=UniRef90_A0A3E5DVA8)
write.table(rev_dens,"/home/rsm34/5asabb3/reviews/keyfiles/rev_dens.txt",sep="\t",quote=FALSE,row.names=FALSE)


#################### 
# adjoing barplots #
####################
########
u_a <- unirefs_plotting %>% select(matches("UniRef90_R6TIX3|binary154")) %>% group_by(binary154) %>% summarize(pctzero=sum(UniRef90_R6TIX3==0)/nrow(unirefs_plotting)*100) %>% mutate(feature_name = "UniRef90_R6TIX3")
u_b <- unirefs_plotting %>% select(matches("UniRef90_C7H1G6|binary154")) %>% group_by(binary154) %>% summarize(pctzero=sum(UniRef90_C7H1G6==0)/nrow(unirefs_plotting)*100) %>% mutate(feature_name = "UniRef90_C7H1G6")
u_c <- unirefs_plotting %>% select(matches("UniRef90_R6CZ24|binary154")) %>% group_by(binary154) %>% summarize(pctzero=sum(UniRef90_R6CZ24==0)/nrow(unirefs_plotting)*100) %>% mutate(feature_name = "UniRef90_R6CZ24")
u_d <- unirefs_plotting %>% select(matches("UniRef90_R5CY66|binary154")) %>% group_by(binary154) %>% summarize(pctzero=sum(UniRef90_R5CY66==0)/nrow(unirefs_plotting)*100) %>% mutate(feature_name = "UniRef90_R5CY66")
u_e <- unirefs_plotting %>% select(matches("UniRef90_T5S060|binary154")) %>% group_by(binary154) %>% summarize(pctzero=sum(UniRef90_T5S060==0)/nrow(unirefs_plotting)*100) %>% mutate(feature_name = "UniRef90_T5S060")
u_f <- unirefs_plotting %>% select(matches("UniRef90_U2QX51|binary154")) %>% group_by(binary154) %>% summarize(pctzero=sum(UniRef90_U2QX51==0)/nrow(unirefs_plotting)*100) %>% mutate(feature_name = "UniRef90_U2QX51")

uni_bar_outfile <- do.call("rbind",list(u_a,u_b,u_c,u_d,u_e,u_f))
write.table(uni_bar_outfile,"/home/rsm34/5asabb3/reviews/keyfiles/uni_bar_outfile_rev.txt",sep="\t",quote=FALSE,row.names=FALSE)

rev_bar <- unirefs_plotting %>% select(matches("UniRef90_A0A3E5DVA8|binary154")) %>% group_by(binary154) %>% summarize(pctzero=sum(UniRef90_A0A3E5DVA8==0)/nrow(unirefs_plotting)*100) %>% mutate(feature_name = "UniRef90_A0A3E5DVA8")
write.table(rev_bar,"/home/rsm34/5asabb3/reviews/keyfiles/rev_bar.txt",sep="\t",quote=FALSE,row.names=FALSE)

#### enrichment nominal p values calcuated with these proportions and 119 non-users and 94 users

################# other ##############
rnauni_rotated_metab3 <- rnauni_rotated_metab2 %>% filter(binary154==1)
acetyl196 <- rnauni_rotated_metab3 %>% select(matches("UniRef90_R6TIX3|UniRef90_C7H1G6|X196.0609_2.81|UniRef90_A0A1C7H380|UniRef90_A0A076IMG8|UniRef90_A0A076IRG9|UniRef90_A0A076IXU4|UniRef90_C7H1G6|UniRef90_I9R3P9|UniRef90_A0A0P0M2W7|UniRef90_U2QX51|UniRef90_R6CZ24|UniRef90_R5CY66|UniRef90_T5S060")) %>% mutate(Nacetyl5ASA=log(X196.0609_2.81+1))

kmetabs2 = NULL
kmetabs2 = as.data.frame(matrix(-9,ncol(acetyl196 ),4))
colnames(kmetabs2) = c("feature","beta","p","R2")

for(i in c(2:13))
	{	b<-lm(acetyl196$Nacetyl5ASA~acetyl196[,i])
 		kmetabs2[i,1]= colnames(acetyl196 )[i]
 		kmetabs2[i,2]=summary(b)$coefficients[2,1] 
 		kmetabs2[i,3]=summary(b)$coefficients[2,4]
 		kmetabs2[i,4]=summary(b)$r.squared
	}
	kmetabs2$FDR <- p.adjust(kmetabs2$p,method="fdr")
	kmetabs2 <- kmetabs2 %>% arrange(desc(R2)) #interestingly, some of these have a negative Beta... 

gather196sig <- acetyl196 %>% gather(gene,arcsqt_abundance,2:13)

write.table(gather196sig ,"acetyl196plot_outfile.txt",sep="\t",quote=FALSE,row.names=FALSE)

########################################################################################### PLOTTING ##############################################################################################


################### PANEL A: DENSITY + BAR PLOTS #################### 

library(cowplot)

#density
uni_dens_outfile <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/uni_dens_outfile_rev.txt",header=T,sep="\t")
uni_dens_outfile$binary154 <- as.factor(uni_dens_outfile$binary154)
uni_dens_outfile$feature_name = factor(uni_dens_outfile$feature_name, levels = c("UniRef90_C7H1G6", "UniRef90_R6TIX3", "UniRef90_R6CZ24", "UniRef90_R5CY66","UniRef90_T5S060","UniRef90_U2QX51"))
my_col_med = c("0"="#3C5488FF","1"="#DC0000FF") 
unirefdens <-ggplot(uni_dens_outfile, aes(scaled, fill = binary154,color=binary154)) +
  geom_density(alpha = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+facet_wrap(~feature_name,dir="v",ncol=1) + theme_cowplot(14) + theme(strip.background = element_blank(),legend.position="top") + geom_rug() + 
scale_fill_manual(values=my_col_med,labels=c("Non-users","Users"),name = "5-ASA use") + scale_color_manual(values=my_col_med,guide = FALSE) #4x12 pdf 

ggsave(filename='/home/rsm34/5asabb3/reviews/keyfiles/fig3a_dens.png', plot=unirefdens, width = 3.6, height = 10.1, dpi = 600) 


#Bars
uni_bar_outfile <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/uni_bar_outfile_rev.txt",header=T,sep="\t")
uni_bar_outfile$binary154 <- as.factor(uni_bar_outfile$binary154)

uni_bar_outfile$feature_name = factor(uni_bar_outfile$feature_name, levels = c("UniRef90_C7H1G6", "UniRef90_R6TIX3", "UniRef90_R6CZ24", "UniRef90_R5CY66","UniRef90_T5S060","UniRef90_U2QX51"))

unirefbar <-ggplot(uni_bar_outfile, aes(x=binary154, y=pctzero, fill=binary154)) +
  geom_bar(stat="identity",width=0.6,alpha=0.5)+ theme_cowplot(14) +facet_wrap(~feature_name,dir="v",ncol=1) + 
theme(strip.background = element_blank(),strip.text.x = element_blank(),aspect.ratio = 2/1,axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values=my_col_med,labels=c("Non-users","Users"),guide = FALSE) + scale_y_continuous(name="% zero",breaks=c(0,25,50))
 
ggsave(filename='/home/rsm34/5asabb3/reviews/keyfiles/fig3a_bar.png', plot=unirefbar, width = 1.5, height = 10, dpi = 600) 


### rev fig ####

#density
rev_dens_outfile <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/rev_dens.txt",header=T,sep="\t")
rev_dens_outfile$binary154 <- as.factor(rev_dens_outfile$binary154)
my_col_med = c("0"="#3C5488FF","1"="#DC0000FF") 
unirefdens <-ggplot(uni_dens_outfile, aes(scaled, fill = binary154,color=binary154)) +
  geom_density(alpha = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme_cowplot(14) + theme(strip.background = element_blank(),legend.position="top") + geom_rug() + 
scale_fill_manual(values=my_col_med,labels=c("Non-users","Users"),name = "5-ASA use") + scale_color_manual(values=my_col_med,guide = FALSE) #4x12 pdf 

ggsave(filename='/home/rsm34/5asabb3/reviews/keyfiles/rev_dens.png', plot=unirefdens, width = 3.6, height = 2, dpi = 600) 


#Bars
rev_bar_outfile <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/rev_bar.txt",header=T,sep="\t")
rev_bar_outfile$binary154 <- as.factor(rev_bar_outfile$binary154)

rev_bar_outfile$feature_name = factor(rev_bar_outfile$feature_name, levels = c("UniRef90_C7H1G6", "UniRef90_R6TIX3", "UniRef90_R6CZ24", "UniRef90_R5CY66","UniRef90_T5S060","UniRef90_U2QX51"))

unirefbar <-ggplot(rev_bar_outfile, aes(x=binary154, y=pctzero, fill=binary154)) +
  geom_bar(stat="identity",width=0.6,alpha=0.5)+ theme_cowplot(14)  + 
theme(strip.background = element_blank(),strip.text.x = element_blank(),aspect.ratio = 2/1,axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values=my_col_med,labels=c("Non-users","Users"),guide = FALSE) + scale_y_continuous(name="% zero")
 
ggsave(filename='/home/rsm34/5asabb3/reviews/keyfiles/rev_bar.png', plot=unirefbar, width = 1.5, height = 2, dpi = 600) 








############# PANEL C: CORRELATIONS #######################

#make figure in R studio 
gather196sig  <- read.table("/home/rsm34/5asabb3/submission/keyoutput/acetyl196plot_outfile.txt",header=T,sep="\t")

gather196sig2  <- gather196sig %>% filter(gene == "UniRef90_R6TIX3"|gene=="UniRef90_U2QX51"|gene=="UniRef90_T5S060")

library(scales)
library(ggpmisc)

my.formula <- y~x 
unirefcorr <- ggplot(gather196sig2 , aes(x = arcsqt_abundance, y = Nacetyl5ASA)) +
    geom_point(col = "#3C5488FF",alpha=0.8) + theme_bw()+ geom_smooth(method = "lm",
        col = "#C42126", formula = my.formula,
        #se = FALSE,
        size = 1) + facet_wrap(~ gene,ncol=2,scales = "free_x")+ stat_poly_eq(aes(label = paste(..p.value.label..)),
                      label.x.npc = "right", label.y.npc = 0.15, formula = my.formula, 
                      parse = TRUE, size = 4) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
                      label.x.npc = "right", label.y.npc = 0.05, formula = my.formula, 
                      parse = TRUE, size = 4)+theme_cowplot(14) +
theme(strip.background = element_blank(),axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+xlab("Rel. abund. (asin sqrt)")+ylab("N-acetyl 5-ASA")+
scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) #pdf 4x4 

ggsave(filename='/home/rsm34/5asabb3/reviews/fig3c_rev.png', plot=unirefcorr , width = 4, height = 5, dpi = 600) 





