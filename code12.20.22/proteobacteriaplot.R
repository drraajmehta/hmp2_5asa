##########

library(ggplot2)
library(dplyr)
library(tidyr)

#### NEW 3.16.22 --> has unexplained variance #####

phylum <- read.table("/home/rsm34/5asabb3/archive/taxa/metaphlan_rotated_phy.txt",header=T,sep="\t") %>% select(SampleID, Proteobacteria, Bacteroidetes, Firmicutes, Actinobacteria, Verrucomicrobia)


phylum_gather <- phylum %>% gather(Class, Value,-SampleID) 

##########################
# sum the abundance for each class, across all IDs, & sort the result
sort.class <- phylum_gather %>% 
  count(Class, wt = Value) %>%
  arrange(desc(n)) %>%
  pull(Class)

feature.order <- phylum_gather%>%
  filter(Class == sort.class[3]) %>%
  arrange((Value)) %>%
  pull(SampleID)

pp<-phylum_gather %>%
  mutate(feature = factor(SampleID, levels = feature.order)) %>%
  mutate(Class = factor(Class, levels = rev(sort.class))) %>%
  ggplot(aes(x = SampleID, y = Value, fill = Class)) +
  geom_col(width = 1) #geom_col is equivalent to geom_bar(stat = "identity")
pp 


p <- ggplot(df, aes(x=weight)) + geom_area(stat = "bin") +
  scale_fill_brewer(palette="Dark2") 

#####
proteoorder <-phylum %>% arrange(desc(Proteobacteria)) %>% select(SampleID)

proteoorder <- proteoorder[,1]

phylum_gather$SampleID <- factor(phylum_gather$SampleID , levels=proteoorder)

p <- ggplot(phylum_gather, aes(x=SampleID, y=Value, fill=Class)) + 
    geom_area()

pp<-phylum_gather %>%
  mutate(feature = factor(SampleID, levels = proteoorder )) %>%
  mutate(Class = factor(Class, levels = c("Proteobacteria","Bacteroidetes","Firmicutes","Actinobacteria","Verrucomicrobia"))) %>%
  ggplot(aes(x = SampleID, y = Value, fill = Class)) +
  geom_col(width = 1)+ ylab("Relative abundance") + scale_fill_manual(values = c("#00A087FF","#3C5488FF","#E64B35FF","purple","grey"),name = "Phylum")+ xlab("Sample") +
 theme( axis.text.x = element_blank(), axis.ticks = element_blank())


ggsave(filename='./protoplot.png', plot=pp , width = 6, height = 4, dpi = 600)
