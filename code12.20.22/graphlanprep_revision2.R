library(dplyr)

#read in simmer data 
########################
#read in the SPECIES file after pre-processing script made it (1632 --> n=1464) 
simmerdata <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/DM1_enzyme_predictions.txt",header=T,sep="\t")

simmerdata$Lineage <- gsub(".*p__","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("c__","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("o__","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("f__","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("g__","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("s__","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("\\;",".",simmerdata$Lineage)
simmerdata$Lineage <- gsub("\\_A","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("\\_B","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("\\_C","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("\\_D","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("\\_I","",simmerdata$Lineage)
simmerdata$Lineage <- gsub("\\_G","",simmerdata$Lineage)
simmerdata$Lineage <- sub(" ", "_", simmerdata$Lineage)

simmerfile <- simmerdata %>% select(3) %>% arrange(Lineage) %>% rename(species = "Lineage") 
simmerfile <-subset(simmerfile , !(endsWith(species, '.')))

simmerfile <- simmerfile %>% distinct(species) #77-> 62 

#read in MGX data 
######################
$ get names of organisms that contribute MGX uniref90s
$cut -f 1 /home/rsm34/5asabb3/archive/DNA_uniref/dna_merge_genefams.txt > /home/rsm34/5asabb3/reviews/keyfiles/dnabugnames.txt
$grep "g_" /home/rsm34/5asabb3/reviews/keyfiles/dnabugnames.txt > /home/rsm34/5asabb3/reviews/keyfiles/dnabugnames2.txt

$ now map them to full taxonomic names from unfiltered metaphlan file 
$cut -f 1 /home/rsm34/5asabb3/archive/taxa/metaphlan_species_raw.txt > /home/rsm34/5asabb3/reviews/keyfiles/rawbugnames.txt 

#read into R and make graphlan compatible
library(dplyr)
dnabugs <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/dnabugnames2.txt", header=F, sep="\t")
colnames(dnabugs) <- c("uniref") 
dnabugs$species <- gsub(".*\\|","",dnabugs$uniref) 
dnabugs$uniref <- gsub("\\|.*","",dnabugs$uniref) 

#read in filtered uniref90 names (10e-8 etc filtration) 
filtuniref <- read.table("/home/rsm34/5asabb3/reviews/uniref90filtMGXnames.txt",header=T,sep="\t")
colnames(filtuniref) <- "uniref"

dnafiltbugs <- inner_join(dnabugs,filtuniref, by="uniref")

dnafiltbugs <- dnafiltbugs %>% distinct(species)

rawbugs <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/rawbugnames.txt",header=T,sep="\t") 
rawbugs$Sample <- gsub("\\|",".",rawbugs$Sample)
rawbugs$species <- gsub(".*f__","",rawbugs$Sample) 
rawbugs$species <- sub("^[^.]*.", "", rawbugs$species)

bugMGX <- inner_join(rawbugs,dnafiltbugs , by="species") #322 
bugMGX <- bugMGX %>% select(1) %>% rename(species=Sample)

bugMGX <- bugMGX %>% distinct(species)  #n=322


#read in MTX data
#################
$ get names of organisms that contribute MTX uniref90s
$cut -f 1 /home/rsm34/5asabb3/archive/RNA_uniref/rna_merge_genefams.txt > /home/rsm34/5asabb3/reviews/keyfiles/rnabugnames.txt
$grep "g_" /home/rsm34/5asabb3/reviews/keyfiles/rnabugnames.txt > /home/rsm34/5asabb3/reviews/keyfiles/rnabugnames2.txt

$ now map them to full taxonomic names from unfiltered metaphlan file 
$cut -f 1 /home/rsm34/5asabb3/archive/taxa/metaphlan_species_raw.txt > /home/rsm34/5asabb3/reviews/keyfiles/rawbugnames.txt 

#read into R and make graphlan compatible
library(dplyr)
rnabugs <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/rnabugnames2.txt", header=F, sep="\t")
colnames(rnabugs) <- c("uniref") 
rnabugs$species <- gsub(".*\\|","",rnabugs$uniref) 
rnabugs$uniref <- gsub("\\|.*","",rnabugs$uniref) 

#read in filtered uniref90 names (10e-8 etc filtration) 
filtunirefrna <- read.table("/home/rsm34/5asabb3/reviews/uniref90filtMTXnames.txt",header=T,sep="\t")
colnames(filtunirefrna) <- "uniref"

rnafiltbugs <- inner_join(rnabugs,filtunirefrna, by="uniref")

rnafiltbugs <- rnafiltbugs %>% distinct(species)

rawbugs <- read.table("/home/rsm34/5asabb3/reviews/keyfiles/rawbugnames.txt",header=T,sep="\t") 
rawbugs$Sample <- gsub("\\|",".",rawbugs$Sample)
rawbugs$species <- gsub(".*f__","",rawbugs$Sample) 
rawbugs$species <- sub("^[^.]*.", "", rawbugs$species)

bugMTX <- inner_join(rawbugs,rnafiltbugs, by="species") #330 
bugMTX <- bugMTX %>% select(1) %>% rename(species=Sample)

bugMTX <- bugMTX %>% distinct(species) #301 


##read in delomenie bugs####
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC99640/

delomenie <- c(
"k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Salmonella|s__Salmonella_enterica", 
"k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter|s__Citrobacter_amalonaticus",
"k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter|s__Citrobacter_farmeri",
"k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter|s__Citrobacter_freundii",
"k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter|s__Citrobacter_koseri",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae.Escherichia.Escherichia_coli",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae.Klebsiella.Klebsiella_pneumoniae",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae.Klebsiella.Klebsiella_oxytoca",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae.Klebsiella.Klebsiella_rhinoscleromatis",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Morganellaceae.Morganella.Morganella_morganii",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae.Proteus.Proteus_vulgaris",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Morganellaceae.Providencia.Providencia_alcalifaciens",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Yersiniaceae.Serratia.Serratia_marcescens",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae.Shigella.Shigella_flexneri",
"Proteobacteria.Gammaproteobacteria.Aeromonadales.Aeromonadaceae.Aeromonas.Aeromonas_hydrophila",
"k__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Helicobacteraceae|g__Helicobacter|s__Helicobacter_pylori",
"Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae.Plesiomonas.Plesiomonas_shigelloides",
"Proteobacteria.Gammaproteobacteria.Vibrionales.Vibrionaceae.Vibrio.Vibrio_cholerae",
"k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Alcaligenaceae|g__Achromobacter|s__Achromobacter_xylosoxidans",
"Proteobacteria.Gammaproteobacteria.Pseudomonadales.Pseudomonadaceae.Pseudomonas.Pseudomonas_aeruginosa")

delomenie <- data.frame(delomenie)
colnames(delomenie) <- c("species") 


### now concatenate ####

temp1<-rbind(bugMGX, bugMTX)
graphlanring <- rbind(temp1, delomenie) 

graphlanring$species <- gsub(".*p__","",graphlanring$species)
graphlanring$species <- gsub("c__","",graphlanring$species)
graphlanring$species <- gsub("o__","",graphlanring$species)
graphlanring$species <- gsub("f__","",graphlanring$species)
graphlanring$species <- gsub("g__","",graphlanring$species)
graphlanring$species <- gsub("s__","",graphlanring$species)
graphlanring$species <- gsub("\\;",".",graphlanring$species)
graphlanring$species <- gsub("\\|",".",graphlanring$species)


graphlanring$species <- gsub("Actinobacteriota","Actinobacteria",graphlanring$species) #have to then add 1 manually to prevent overlap with class of same name 

graphlanring$species  <- gsub("a.Actinobacteria","a.Actinobacteria1",graphlanring$species)

graphlaninput <- graphlanring %>% distinct(species) %>% arrange(species)

write.table(graphlaninput , file="graphlanrevise.txt", sep="\t", quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(graphlanring, file="graphlanring2.txt", sep="\t", quote=FALSE,row.names=FALSE,col.names=FALSE)

#####make ring file in excel with formatting #####

#then go to galaxy and upload files
#add in phyla annotation (white color)
#shrink nodes to size 4 (marker) 
#add in 6 bugs at size 10 (black marker)
#add in ring file 
#edit in adobe 








