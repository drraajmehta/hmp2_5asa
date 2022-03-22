
############################################################### CODE ###################################################

library(dplyr)
library(nlme)
library(tidyr) 
library(stringr)

#read in the files 
###################################

#bring in metabolite data
########################

#prepare metabolite data
#read in METABOLOMICS file after pre-processing script made it 
metab_select <- read.table("/home/rsm34/5asabb3/archive/metabolite/metabs_export.txt", header=T, sep="\t") %>% select(c("SampleID","binary154","X196.0609_2.81")) #saves time in the future. 

#bring in species data
########################
#read in the SPECIES file after pre-processing script made it (1632 --> n=1464)
species_rotated <- read.table("/home/rsm34/5asabb3/archive/taxa/metaphlan_rotated.txt",header=T,sep="\t")
species_rotated<-species_rotated[rowSums(species_rotated[,2:ncol(species_rotated)])!=0,,drop=FALSE] 	#there are a few rows that sum to zero (also a few that sum to <0.5

#build metadata file
#####################
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


######## now merge ###########
metab_species <- inner_join(metab_select,species_rotated,by="SampleID") 
dim(metab_species)
metab_species_meta <- inner_join(meta_data_IBD,metab_species,by="SampleID")
dim(metab_species_meta)

metab_species_meta[,33:ncol(metab_species_meta)] <- asin(sqrt(metab_species_meta[,33:ncol(metab_species_meta)]*0.01)) #need to scale to 1 for transformation

##### nlme function ######
##########################

#subset to IBD users only: test binary154 (use vs non use)

spec_IBD <- metab_species_meta  

	feat.names <- colnames(spec_IBD)[33:ncol(spec_IBD)]
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
  			tmp <- spec_IBD[, c(i,"binary154","Antibiotics","consent_age","dysbiosis","site_name","Participant.ID","diagnosis")]
  			fml <- as.formula( paste( i, "~", paste(c("binary154","Antibiotics","consent_age","dysbiosis","diagnosis"), collapse="+") ) )

  			# assign fit to list by name
  			b <- lme(fml, random=~1|site_name/Participant.ID, data=tmp,na.action=na.omit,control=lmeControl(returnObject=TRUE))
 			diag_sig[i,1]= feat.names[i]
 			diag_sig[i,2]=(b$coefficients)$fixed[2]
 			diag_sig[i,3]=summary(b)$tTable[2,5]
			diag_sig$names=rownames(diag_sig)
		}
	diag_sig$FDR <- p.adjust(diag_sig$p,method="fdr")
	diag_sig <- diag_sig %>% arrange((p)) 
	row.names(diag_sig) <- NULL 
	diag_sig$names <- gsub(".*s__","",diag_sig$names)
	write.table(diag_sig, file="spec_ibd_nlme.txt", sep="\t", quote=FALSE,row.names=FALSE)


#try with 196 now 
	spec_IBD2 <- spec_IBD %>% filter(binary154==1) 
	spec_IBD2$X196.0609_2.81 <- log(spec_IBD2$X196.0609_2.81)
	feat.names <- colnames(spec_IBD2)[33:ncol(spec_IBD2)]
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
  			tmp <- spec_IBD2[, c(i,"X196.0609_2.81","Antibiotics","consent_age","dysbiosis","site_name","Participant.ID","diagnosis")]
  			fml <- as.formula( paste( i, "~", paste(c("X196.0609_2.81","Antibiotics","consent_age","dysbiosis","diagnosis"), collapse="+") ) )

  			# assign fit to list by name
  			b <- lme(fml, random=~1|site_name/Participant.ID, data=tmp,na.action=na.omit,control=lmeControl(returnObject=TRUE))
 			diag_sig[i,1]= feat.names[i]
 			diag_sig[i,2]=(b$coefficients)$fixed[2]
 			diag_sig[i,3]=summary(b)$tTable[2,5]
			diag_sig$names=rownames(diag_sig)
		}
	diag_sig$FDR <- p.adjust(diag_sig$p,method="fdr")
	diag_sig <- diag_sig %>% arrange((p)) 
	row.names(diag_sig) <- NULL 
	diag_sig$names <- gsub(".*s__","",diag_sig$names)
	write.table(diag_sig, file="spec_196_nlme.txt", sep="\t", quote=FALSE,row.names=FALSE)












#remants of experimentation 


#dnauni<- dnauni_rotated_metab[,1:20]
#rnauni<- rnauni_rotated_metab[,1:20]


a1=which(colnames(dnauni)=="UniRef90_A0A010YBU9")
a2=which(colnames(rnauni)=="UniRef90_A0A010YBU9")
a3=as.data.frame(dnauni[,5])
