module load gcc/6.2.0 R/4.0.1
library(dplyr)

########################### METADATA ######################################################
$ wget https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv --no-check-certificate   

########################### TAXONOMIC DATA ######################################################
$ wget https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/taxonomic_profiles_3.tsv.gz --no-check-certificate
$ gunzip taxonomic_profiles_3.tsv.gz
$ sed 's/Feature\\//' taxonomic_profiles_3.tsv > taxonomic_profiles.tsv #gets rid of hashtag on first line 

#filter to species only 
$ grep "|s__" taxonomic_profiles.tsv > speciesprecursor1.tsv
$ head -1 taxonomic_profiles.tsv > column.header
$ cat column.header speciesprecursor1.tsv > speciesprecursor2.tsv
$ grep -v "Vir" speciesprecursor2.tsv > speciesprecursor3.txt
$ grep -v "|t__" speciesprecursor3.txt > metaphlan_species.txt

# remove taxa with relative abundance below 1e-4 in at least 10 % samples --> R (stopped here)
metaphlan_species <- read.table("/home/rsm34/5asabb3/archive/metaphlan_species.txt", header=T, sep="\t")
metaphlan_species2 <- metaphlan_species[ rowSums(metaphlan_species >= 0.0001 ) >= 150, ]
dim(metaphlan_species2)
dim(metaphlan_species)
write.table(metaphlan_species2,  "metaphlan_species2.txt", sep = "\t", ,row.names=FALSE,quote=F )

$ cat metaphlan_species2.txt > speciesbiom.txt 

#rotates file (in R)
metaphlan_species <- read.table("/home/rsm34/5asabb3/archive/speciesbiom.txt", header=T, sep="\t")
metaphlan_rotated = setNames(data.frame(t(metaphlan_species[,-1])), metaphlan_species[,1])
metaphlan_rotated<- cbind (SampleID=row.names(metaphlan_rotated),metaphlan_rotated)
row.names(metaphlan_rotated) <- NULL
metaphlan_rotated <- metaphlan_rotated[ ,c(1,ncol(metaphlan_rotated),2:(ncol(metaphlan_rotated)-1))]
metaphlan_rotated$SampleID <- substr(metaphlan_rotated$SampleID,1,nchar(metaphlan_rotated$SampleID)-8)

write.table(metaphlan_rotated,"metaphlan_rotated.txt",sep="\t",row.names=FALSE,quote=F)


###### phylum level data ########

#To Be done in UNIX because loading large files overwhelms R -- the follow code provides phlyum data# 

$ grep "|p__" taxonomic_profiles.tsv > phlyumprecursor1.tsv
$ head -1 taxonomic_profiles.tsv > column.header
$ cat column.header phlyumprecursor1.tsv > phlyumprecursor2.tsv
$ grep -v "Vir" phlyumprecursor2.tsv > phlyumprecursor3.txt
$ grep -v "|c__" phlyumprecursor3.txt > metaphlan_phlyum.txt 

#R
#now we shorten the phlyum names for graphing purposes -- this is hardly elegant 
metaphlan_phlyum <- read.table("/home/rsm34/5asa/archive/metaphlan_phlyum.txt", header=T, sep="\t", check.names=FALSE, row.names=1)
longphlyum <- row.names(metaphlan_phlyum)
SampleID <-  gsub(".*p__","",longphlyum)
rownames(metaphlan_phlyum) <- NULL
phlyumfile <- cbind(SampleID,metaphlan_phlyum)
metaphlan_rotatedP = setNames(data.frame(t(phlyumfile [,-1])), phlyumfile [,1])
metaphlan_rotatedP<- cbind (SampleID=row.names(metaphlan_rotatedP),metaphlan_rotatedP)
row.names(metaphlan_rotatedP) <- NULL
metaphlan_rotatedP <- metaphlan_rotatedP[ ,c(1,ncol(metaphlan_rotatedP),2:(ncol(metaphlan_rotatedP)-1))] #reorders matrix
write.table(metaphlan_rotatedP,"metaphlan_rotated_phy.txt",sep="\t",row.names=FALSE,quote=F)

############################# BIOPSY 16S ########################################################
$ wget https://ibdmdb.org/tunnel/products/HMP2_Pilot/16S/1547/taxonomic_profiles.tsv.gz --no-check-certificate
$ gunzip taxonomic_profiles.tsv.gz
#make sure to delete # from first row 
$ sed 's/^#//' taxonomic_profiles.tsv > 16sbiopsy_pre1.txt #gets rid of hashtag on first line 

#make prop table file (i checked and the rows don't sum to 1)
      # in R 
      data_16s <- read.table("/home/rsm34/5asabb3/archive/16sbiopsy_pre1.txt", header=T, sep="\t")
      data_16s [,2:ncol(data_16s )] <- prop.table(as.matrix(data_16s [,2:ncol(data_16s )]),2)
      write.table(data_16s ,"data_16s_pre2.txt",sep="\t", row.names=FALSE,quote=F)

$ sed -i -e 's/Clade_name/ID/g' data_16s_pre2.txt

#### rotates files in R #####

data_16s <- read.table("/home/rsm34/5asa/archive/data_16s_pre2.txt", header=T, sep="\t", check.names=FALSE)
data_16s_rotated = setNames(data.frame(t(data_16s[,-1])), data_16s[,1])
data_16s_rotated<- cbind (SampleID=row.names(data_16s_rotated),data_16s_rotated)
row.names(data_16s_rotated) <- NULL
write.table(data_16s_rotated,"data_16s_rotated.txt",sep="\t",row.names=FALSE,quote=F)


########################### DNA PATHWAYS ######################################################  HAVE NOT LOADED B/C USING UNIREFS #####
$ wget https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/pathabundances.tsv.gz --no-check-certificate
$ gunzip pathabundances.tsv.gz

$ grep -v "g_" pathabundances.tsv > DNAPathPrecursor1.txt
$ grep -v "|" DNAPathPrecursor1.txt > DNAPathPrecursor2.txt
$ sed 's/:[^\t]*//' DNAPathPrecursor2.txt > DNAPathPrecursor3.txt #gets rid of english names after colon 

$ cat DNAPathPrecursor3.txt > DNAPathbiom.txt
$ sed -i -e 's/Pathway/SampleID/g' DNAPathbiom.txt
#may need to delete # from to cell 1,row1

#rotates file (in R)
DNApath <- read.table("/home/rsm34/5asa/DNAPathbiom.txt", header=T, sep="\t", check.names=FALSE)
DNApath_rotated = setNames(data.frame(t(DNApath[,-1])), DNApath[,1])
DNApath_rotated<- cbind (SampleID=row.names(DNApath_rotated),DNApath_rotated)
row.names(DNApath_rotated) <- NULL
DNApath_rotated <- DNApath_rotated[ ,c(1,ncol(DNApath_rotated),2:(ncol(DNApath_rotated)-1))]
write.table(DNApath_rotated,"DNApath_rotated.txt",sep="\t",row.names=FALSE,quote=F)

########################### RNA PATHWAYS ########################################################  HAVE NOT LOADED B/C USING UNIREFS #######
$ cd /home/rsm34/5asa/transcript 
$ wget https://ibdmdb.org/tunnel/products/HMP2/MTX/1750/pathabundances.tsv.gz --no-check-certificate
$ gunzip pathabundances.tsv.gz

$ grep -v "g_" pathabundances.tsv > RNAPathPrecursor1.txt
$ grep -v "|" RNAPathPrecursor1.txt > RNAPathPrecursor2.txt
$ sed 's/:[^\t]*//' RNAPathPrecursor2.txt > RNAPathPrecursor3.txt #gets rid of english names after colon 

$ cat RNAPathPrecursor3.txt > RNAPathbiom.txt
$ sed -i -e 's/Pathway/SampleID/g' RNAPathbiom.txt
#may need to delete # from to cell 1,row1

RNApath <- read.table("/home/rsm34/5asa/transcript/RNAPathbiom.txt", header=T, sep="\t", check.names=FALSE)
RNApath_rotated = setNames(data.frame(t(RNApath[,-1])), RNApath[,1])
RNApath_rotated<- cbind (SampleID=row.names(RNApath_rotated),RNApath_rotated)
row.names(RNApath_rotated) <- NULL
RNApath_rotated <- RNApath_rotated[ ,c(1,ncol(RNApath_rotated),2:(ncol(RNApath_rotated)-1))]
write.table(RNApath_rotated,"RNApath_rotated.txt",sep="\t",row.names=FALSE,quote=F)


########################### DNA: EC Data #######################################################  HAVE NOT LOADED B/C USING UNIREFS ########
$ wget https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/ecs.tsv.gz --no-check-certificate
$ gunzip ecs.tsv.gz

# GREP not: 
$ grep -v "g_" ecs.tsv > GenePrecursor1.txt #keeps only the ECs/pathways
$ grep -v "|" GenePrecursor1.txt > GenePrecursor2.txt #gets rid of unclassified
$ sed 's/:[^\t]*//' GenePrecursor2.txt > GenePrecursor3.txt #gets rid of english names after colon 

#make prop table file (i checked and the rows don't sum to 1)
      # in R (FOR DNA)
      dnabiom <- read.table("/home/rsm34/5asa/GenePrecursor3.txt", header=T, sep="\t", check.names=FALSE)
      dnabiom[,2:ncol(dnabiom)] <- prop.table(as.matrix(dnabiom[,2:ncol(dnabiom)]),2)
      write.table(dnabiom,"Genebiom.txt",sep="\t", row.names=FALSE,quote=F)
$ sed -i -e 's/Gene Family/ID/g' Genebiom.txt

#### rotates files in R #####

gene <- read.table("/home/rsm34/5asa/Genebiom.txt", header=T, sep="\t", check.names=FALSE)
gene_rotated = setNames(data.frame(t(gene[,-1])), gene[,1])
gene_rotated<- cbind (SampleID=row.names(gene_rotated),gene_rotated)
row.names(gene_rotated) <- NULL
write.table(gene_rotated,"gene_rotated.txt",sep="\t",row.names=FALSE,quote=F)

gene_rotated <- read.table("/home/rsm34/5asa/gene_rotated.txt",header=T,sep="\t")

gene_rotated$SampleID <- as.character(gene_rotated$SampleID)



########################### RNA: EC Data ######################################################  HAVE NOT LOADED B/C USING UNIREFS #########
$ wget https://ibdmdb.org/tunnel/products/HMP2/MTX/1750/ecs.tsv.gz --no-check-certificate
$ gunzip ecs.tsv.gz

# GREP not: 
$ grep -v "g_" ecs.tsv > RNAPrecursor1.txt #keeps only the ECs/pathways
$ grep -v "|" RNAPrecursor1.txt > RNAPrecursor2.txt #gets rid of unclassified
$ sed 's/:[^\t]*//' RNAPrecursor2.txt > RNAPrecursor3.txt #gets rid of english names after colon 
$ sed 's/^#.//' RNAPrecursor3.txt > RNAPrecursor4.txt #gets rid of hashtag on first line 


#make prop table file (i checked and the rows don't sum to 1)
      # in R (FOR DNA)
      rnabiom <- read.table("/home/rsm34/5asa/transcript/RNAPrecursor4.txt", header=T, sep="\t", check.names=FALSE)
      rnabiom[,2:ncol(rnabiom)] <- prop.table(as.matrix(rnabiom[,2:ncol(rnabiom)]),2)
      write.table(rnabiom,"rnabiom.txt",sep="\t", row.names=FALSE,quote=F)
$ sed -i -e 's/Gene Family/ID/g' rnabiom.txt

#### rotates files in R #####

rna <- read.table("/home/rsm34/5asa/transcript/rnabiom.txt", header=T, sep="\t", check.names=FALSE)
rna_rotated = setNames(data.frame(t(rna[,-1])), rna[,1])
rna_rotated<- cbind (SampleID=row.names(rna_rotated),rna_rotated)
row.names(rna_rotated) <- NULL
write.table(rna_rotated,"rna_rotated.txt",sep="\t",row.names=FALSE,quote=F)


############################# HOST BIOPSY TRANSCRIPTS ######################################################## NO UPDATE, TRANSFERRED FROM PRIOR FOLDER #### 
$ wget https://ibdmdb.org/tunnel/products/HMP2/HTX/1730/host_tx_counts.tsv.gz --no-check-certificate
$ gunzip host_tx_counts.tsv.gz

#### rotates files in R #####
host_mtx <- read.table("/home/rsm34/5asa/archive/host_tx_counts.tsv", header=T, sep="\t", check.names=FALSE)
host_mtx$ID <- row.names(host_mtx) 
host_mtx <- host_mtx %>% select(ID, everything())
row.names(host_mtx) = NULL 
write.table(host_mtx,"host_mtx.txt",sep="\t",row.names=FALSE, quote=F)
host_mtx_rotated = setNames(data.frame(t(host_mtx[,-1])), host_mtx[,1])
host_mtx_rotated<- cbind (SampleID=row.names(host_mtx_rotated),host_mtx_rotated)
row.names(host_mtx_rotated) <- NULL
write.table(host_mtx_rotated,"host_mtx_rotated.txt",sep="\t",row.names=FALSE,quote=F)


########################### UNIREF DNA ######################################################

##bb2 was used first in this project, with this file from Eric/Yancong 

$ gunzip /home/rsm34/5asa/archive/uniref/HMP2_MGX_genefamilies.tsv.gz

# GREP not: 
$ grep -v "g_" HMP2_MGX_genefamilies.tsv > UnirefPrecursor1.txt #keeps only the unirefs 
$ grep -v "|" UnirefPrecursor1.txt > UnirefPrecursor2.txt #gets rid of unclassified
$ sed 's/^#.//' UnirefPrecursor2.txt > UnirefPrecursor3.txt #gets rid of hashtag on first line
$ wc -l UnirefPrecursor3.txt #1767839 
$ awk -F $'\t' '{c = 0; for(i=1; i<=NF; i++) {c += $i == "0" ? 1 : 0}} c <= 1475' UnirefPrecursor3.txt > UnirefPrecursor4.txt #get rid of features with more than 90% zeroes  
#at this point there are 354232 rows  

$ sed -i -e 's/Gene Family/ID/g' UnirefPrecursor4.txt #changes name 


###bb3 for MTX/DNA analysis###

#ids from dna file (n=213) were extracted from bb2 file (overlap mbx, mtx, mgx) --> see submission/keyoutput/dnaids.txt
##then navigate to /n/home02/rmehta/hmp2bb3/ in fas where a list of these files was downloaded 
##then they were unzipped and filtered using "rm *level4*", "rm *pathabund*", "rm *log*" 
##humann_join_tables --input /n/home02/rmehta/hmp2bb3/mgx --output dna_merge_genefams.txt

# GREP not: 
$ grep -v "g_" dna_merge_genefams.txt > UnirefPrecursor1.txt #keeps only the unirefs 
$ grep -v "|" UnirefPrecursor1.txt > UnirefPrecursor2.txt #gets rid of unclassified
$ sed 's/^#.//' UnirefPrecursor2.txt > UnirefPrecursor3.txt #gets rid of hashtag on first line 
$ wc -l UnirefPrecursor3.txt

#at this point there are 1902185 rows  (genes, unfiltered)

$ sed -i -e 's/Gene Family/ID/g' UnirefPrecursor3.txt #changes name 

###bb3 for clinical DNA analysis###

#ids from clinical data were extracted --> see submission/keyoutput/dnaremainder.txt
##then navigate to /n/home02/rmehta/hmp2bb3/dnaremainder.txt in fas
##then they were unzipped and filtered using "rm *level4*", "rm *pathabund*", "rm *log*" 
##humann_join_tables --input /n/home02/rmehta/hmp2bb3/mgx --output dna_merge_genefams2.txt
##humann_join_tables --input /n/home02/rmehta/hmp2bb3/mgx/nonibd --output dna_merge_genefams_nonibd.txt
##humann_join_tables --input /n/home02/rmehta/hmp2bb3/mgx/no5asa --output dna_merge_genefams_no5asa.txt


# GREP not: 
$ grep -v "g_" dna_merge_genefams2.txt > UnirefPrecursorA.txt #keeps only the unirefs 
$ grep -v "|" UnirefPrecursorA.txt > UnirefPrecursorB.txt #gets rid of unclassified
$ sed 's/^#.//' UnirefPrecursorB.txt > UnirefPrecursorC.txt #gets rid of hashtag on first line 
$ wc -l UnirefPrecursorC.txt

#at this point there are 2305404 rows  (genes, unfiltered)

$ sed -i -e 's/Gene Family/ID/g' UnirefPrecursorC.txt #changes name 


#for no5asa
#ids were extracted --> see submission/keyoutput/IBD_no5asa_mgxids.txt
##then navigate to /n/home02/rmehta/hmp2bb3/mgx in fas

# GREP not: 
$ grep -v "g_" dna_merge_genefams_no5asa.txt > UnirefPrecursorA1.txt #keeps only the unirefs 
$ grep -v "|" UnirefPrecursorA1.txt > UnirefPrecursorB1.txt #gets rid of unclassified
$ sed 's/^#.//' UnirefPrecursorB1.txt > UnirefPrecursorC1.txt #gets rid of hashtag on first line 
$ wc -l UnirefPrecursorC1.txt

#at this point there are 2305404 rows  (genes, unfiltered)

$ sed -i -e 's/Gene Family/ID/g' UnirefPrecursorC1.txt #changes name 

#for noibd
#ids were extracted --> see submission/keyoutput/nonIBDmgxids.txt
# GREP not: 
$ grep -v "g_" dna_merge_genefams_nonibd.txt > UnirefPrecursorA2.txt #keeps only the unirefs 
$ grep -v "|" UnirefPrecursorA2.txt > UnirefPrecursorB2.txt #gets rid of unclassified
$ sed 's/^#.//' UnirefPrecursorB2.txt > UnirefPrecursorC2.txt #gets rid of hashtag on first line 
$ wc -l UnirefPrecursorC2.txt

#at this point there are 2305404 rows  (genes, unfiltered)

$ sed -i -e 's/Gene Family/ID/g' UnirefPrecursorC2.txt #changes name 



########################### UNIREF RNA ######################################################

$ gunzip /home/rsm34/5asa/archive/uniref/MTX/HMP2_MTX_genefamilies_relab.tsv.gz

# GREP not: 
$ grep -v "g_" HMP2_MTX_genefamilies_relab.tsv > UnirefPrecursor1.txt #keeps only the unirefs 
$ grep -v "|" UnirefPrecursor1.txt > UnirefPrecursor2.txt #gets rid of unclassified
$ sed 's/^#.//' UnirefPrecursor2.txt > UnirefPrecursor3.txt #gets rid of hashtag on first line 
$ wc -l UnirefPrecursor3.txt #941782
$ awk -F $'\t' '{c = 0; for(i=1; i<=NF; i++) {c += $i == "0" ? 1 : 0}} c <= 190' UnirefPrecursor3.txt > UnirefPre_MTX.txt #get rid of features with more than 90% zeroes   

$ sed -i -e 's/Gene Family/ID/g' UnirefPre_MTX.txt #changes name 
$ wc -l UnirefPre_MTX.txt #941782

###bb3 ###

#ids from dna file (n=213) were extracted from bb2 file (overlap mbx, mtx, mgx)
##then navigate to /n/home02/rmehta/hmp2bb3/rnaids.txt in fas, these were all grabbed from the internet 
##then they were unzipped and filtered using "rm *level4*", "rm *pathabund*", "rm *log*" 
## then run this: humann_join_tables --input /n/home02/rmehta/hmp2bb3/mtx --output rna_merge_genefams.txt

$ grep -v "g_" rna_merge_genefams.txt > UnirefPrecursor1.txt #keeps only the unirefs 
$ grep -v "|" UnirefPrecursor1.txt > UnirefPrecursor2.txt #gets rid of unclassified
$ sed 's/^#.//' UnirefPrecursor2.txt > UnirefPre_MTX.txt  #gets rid of hashtag on first line 
$ wc -l UnirefPre_MTX.txt

#at this point there are 683944 rows (genes, unfiltered)  

$ sed -i -e 's/Gene Family/ID/g' UnirefPre_MTX.txt #changes name 


############################# METABOLOMICS ########################################################


## metabolomics file is edited somewhat in excel because i moved the m/z ratios to the names for the unnamed compounds and got rid of prime values, was too lazy to do this in R , but verified 12/15/20 ####

metabolomics <- read.table("/home/rsm34/5asa/HMP2_metabolomics_edit.txt", header=T, sep="\t", check.names=FALSE)

#5asa 

metab_5asa <- metabolomics %>% select(-"Pooled QC sample CV",-"m/z",-"HMDB (*Representative ID)",-"Method",-"RT",-"Compound") %>%
		rename(SampleID=Metabolite) %>% 
			filter(SampleID == "carboxyibuprofen" | SampleID == "metronidazole" | SampleID == "salicylate" | SampleID == "sulfapyridine" | SampleID == "154.0502_3.83" |
			SampleID == "132.0454_1.66" | SampleID == "194.046_3.93" | SampleID == "154.0502_3.83" | SampleID == "196.0481_7.88" | SampleID == "196.0533_2.89" |
 SampleID == "196.0579_7.11" | SampleID == "196.0602_5.19" | SampleID == "196.0609_2.81" |
 SampleID == "196.0609_5.03" | SampleID == "196.061_1.91" | SampleID == "196.0612_1.69" | SampleID == "196.0615_5.22" |
 SampleID == "196.0631_4.36" | SampleID == "196.0633_4.86" | SampleID == "196.0647_3.68" | SampleID == "196.0721_7.08" |
 SampleID == "196.073_6.08" | SampleID == "196.0798_4.09" | SampleID == "196.081_3.52" | SampleID == "196.0827_4.94" |
 SampleID == "196.0833_4.8" | SampleID == "196.0837_4.57" | SampleID == "196.0837_5.06" | SampleID == "196.0956_7.06" |
 SampleID == "196.0959_7.05" | SampleID == "196.0968_8.17" | SampleID == "196.097_4.83" | SampleID == "196.0972_2.73" |
 SampleID == "196.0972_2.96" | SampleID == "196.0972_3.71" | SampleID == "196.0972_3.78" | SampleID == "196.0972_4.47" |
 SampleID == "196.0972_5.57" |  SampleID == "196.0974_5.86" | SampleID == "196.0974_7.14" | SampleID == "196.0974_7.25")

metab_rotated = setNames(data.frame(t(metab_5asa[,-1])), metab_5asa[,1])
metab_rotated<- cbind (SampleID=row.names(metab_rotated),metab_rotated)
row.names(metab_rotated) <- NULL
mesal_metab_rotated <- metab_rotated
mesal_metab_rotated$SampleID <- as.character(mesal_metab_rotated$SampleID)

write.table(mesal_metab_rotated,"mesal_metab_rotated2.txt",sep="\t",row.names=FALSE,quote=F)

#all metab rotated 

metab_all <- metabolomics %>% select(-"Pooled QC sample CV",-"m/z",-"HMDB (*Representative ID)",-"Method",-"RT",-"Compound") %>%
		rename(SampleID=Metabolite) 

metab_rotated_all = setNames(data.frame(t(metab_all[,-1])), metab_all[,1])
metab_rotated_all<- cbind (SampleID=row.names(metab_rotated_all),metab_rotated_all)
row.names(metab_rotated_all) <- NULL
metab_rotated_all$SampleID <- as.character(metab_rotated_all$SampleID)

write.table(metab_rotated_all,"metab_rotated_all.txt",sep="\t",row.names=FALSE,quote=F)

metab_rotated <- read.table("/home/rsm34/5asa/metab_rotated_all.txt",header=T,sep="\t")
metab_rotated2 <- metab_rotated %>% mutate_all(~replace(., is.na(.), 0)) #make sure NA is = 0. 
write.table(metab_rotated2 , "metab_data_nozero2.txt", sep="\t", quote=FALSE,row.names=FALSE

#now link with metabolomics data, after pre-processing script made it 
metab_rotated2 <- read.table("/home/rsm34/5asabb3/archive/metabolite/metab_data_nozero2.txt",header=T,sep="\t")
metab_IBD<- inner_join(meta_data_IBD,metab_rotated2,by="SampleID") %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% relocate(binary154, .after = IBDusers)

metab_154196<- metab_rotated2 %>% mutate(binary154 = if_else(X154.0502_3.83 > 10000000,1,0)) %>% relocate(binary154, .after = SampleID) %>% relocate(X154.0502_3.83, .after=binary154) %>% 
		mutate(binary196 = if_else(X196.0609_2.81 > 10000000,1,0)) %>% relocate(binary196, .after = X154.0502_3.83) %>% relocate(X196.0609_2.81, .after=binary196) %>% select(,1:5)
write.table(metab_154196, "metabs_export.txt", sep="\t", quote=FALSE,row.names=FALSE) ### to be used for other analyses as  short cut 

########################## MERGE ##################################

metab_metag <- left_join(metab_rotated,gene_rotated,by="SampleID")
metab_metag <- na.omit(metab_metag)
write.table(metab_metag,"metab_metag.txt",sep="\t",row.names=FALSE,quote=F)




 