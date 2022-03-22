############################################################### CODE ###################################################

library(dplyr)
library(nlme)
library(stringr) 

####### uniref data ############
chunk1 <- read.table("/home/rsm34/5asabb3/archive/DNA_uniref/UnirefPrecursorC.txt", header=T, sep="\t") 
actualColumnNames <- str_replace(names(chunk1),"_Abundance.CPM","")
colnames(chunk1) <- actualColumnNames
chunk1[,2:ncol(chunk1)] <- prop.table(as.matrix(chunk1[,2:ncol(chunk1)]),2)
dim(chunk1) #1902184     214
#write.table(chunk1, file="unirefproptable.txt", sep="\t", quote=FALSE,row.names=FALSE)
#chunk1 <- read.table("unirefproptable.txt",header=T, sep="\t")


#remove features with lt 1e-8 in 10% of samples 
chunk2 <- chunk1[ rowSums(chunk1 >= 0.00000001 ) >= 20, ] #leave this filter as is 
dim(chunk2) 
# 379036    214
write.table(chunk2, file="DNAunirefproptable_clinc.txt", sep="\t", quote=FALSE,row.names=FALSE)
