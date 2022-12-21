library(gggenes)
library(ggplot2)

context <- read.table("genecontext_all3.txt",header=T,sep="\t")
context$molecule<- gsub("n(", "\n(", context$molecule,fixed=TRUE)


#NB the file above was made in excel, from combing through yancongs data 


ggplot(context , aes(xmin = start, xmax = end, y = molecule, fill = gene, 
                          forward = orientation)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position="none")


context %>% filter(molecule == "NMTT01000010 \n(GCA_002549855 \n(Faecalibacterium prausnitzii (A2-165))") %>% 
ggplot(aes(xmin = start, xmax = end, y = molecule, fill = gene, 
                          forward = orientation)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
 scale_fill_brewer(palette = "Set3") +
  theme_genes() 

context %>% filter(molecule == "DS995534.1 \n(GCA_000156075 \n(Phocaeicola dorei)") %>% 
ggplot(aes(xmin = start, xmax = end, y = molecule, fill = gene, 
                          forward = orientation)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
 scale_fill_brewer(palette = "Set3") +
  theme_genes() 

context %>% filter(molecule == "KI271734.1 \n(GCA_000469445 \n(Oscillibacter_sp_KLE_1745)") %>% 
ggplot(aes(xmin = start, xmax = end, y = molecule, fill = gene, 
                          forward = orientation)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
 scale_fill_brewer(palette = "Set3") +
  theme_genes() 



