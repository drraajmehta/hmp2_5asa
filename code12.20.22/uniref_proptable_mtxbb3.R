
############################################################### CODE ################################################### 

library(dplyr)
library(purrr)
library(table1) 
library(tidyr)
library(magrittr)
library(sm)
library(ggplot2)
library(corrplot)
library(caret) 
library(vegan)
library(grid)
library(cowplot)
library(viridis) 
library(nlme)
library(stringr) 

####### uniref data ############
chunk1=NULL
chunk2=NULL 
chunk1 <- read.table("/home/rsm34/5asabb3/archive/RNA_uniref/UnirefPre_MTX.txt", header=T, sep="\t") 
actualColumnNames <- str_replace(names(chunk1),"_Abundance.CPM","")
colnames(chunk1) <- actualColumnNames
chunk1[,2:ncol(chunk1)] <- prop.table(as.matrix(chunk1[,2:ncol(chunk1)]),2)
dim(chunk1) 
#write.table(chunk1, file="MTXunirefproptable.txt", sep="\t", quote=FALSE,row.names=FALSE)
# [1] 683943    214

#chunk1 <- read.table("MTXunirefproptable.txt",header=T, sep="\t")

#remove features with lt 1e-8 in 10% of samples 
chunk2 <- chunk1[ rowSums(chunk1 >= 0.00000001 ) >= 20, ]
dim(chunk2) 
# [1] 61256   214
write.table(chunk2, file="MTXunirefproptable2.txt", sep="\t", quote=FALSE,row.names=FALSE)

