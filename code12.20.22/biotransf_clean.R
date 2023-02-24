######################################################################################################
# R programs for creating "metabolite" plots for part 2				          #
# ============================================					                 #
# Purposes:  user v non user boxplots,    					#
#            							#
######################################################################################################

#references: 
#https://www.genome.jp/kegg-bin/show_pathway?map00760+C05380
#Am J Health-Syst Pharm—Vol 60 Jul 1, 2003 Suppl 2
#structures made in chemdraw https://chemdrawdirect.perkinelmer.cloud/js/sample/index.html#

library(ggplot2)
library(dplyr) 

################# CODE to be done in Rstudio ###############

biotranHMP2 <- read.table("/home/rsm34/5asabb3/submission/keyoutput/metab_biotransf3.txt",header=T,sep="\t") 

biotranPRISM <- read.table("/home/rsm34/5asabb3/submission/keyoutput/metab_biotransf2_prism.txt",header=T,sep="\t")
biotranPRISM$SampleID <- as.factor(biotranPRISM$SampleID)

biotran_both <- rbind(biotranHMP2,biotranPRISM)

biotran_both$binary154 <- as.factor(biotran_both$binary154)

#### FIGURE 2 ######
mesal154 <- biotran_both %>% filter(round=="X154")
mesal154fig<-mesal154 %>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=binary154)) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
  facet_wrap(~ cohort) + scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #4x3 pdf

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesal154fig.png', plot=mesal154fig, width = 4, height = 3, dpi = 600) 


mesal196 <- biotran_both %>% filter(round=="X196")
mesal196fig<-mesal196 %>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=binary154)) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
  facet_wrap(~ cohort) + scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #340x230 png

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesal196fig.png', plot=mesal196fig, width = 3.4, height = 2.3, dpi = 600) 


mesal210 <- biotran_both %>% filter(round=="X210")
mesal210fig<-mesal210 %>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=binary154)) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
  facet_wrap(~ cohort) + scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #340x230 png

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesal210fig.png', plot=mesal210fig, width = 3.4, height = 2.3, dpi = 600) 

mesal224 <- biotran_both %>% filter(round=="X224")
mesal224fig<-mesal224 %>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=binary154)) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
  facet_wrap(~ cohort) + scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #340x230 png

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesal224fig.png', plot=mesal224fig, width = 3.4, height = 2.3, dpi = 600) 

########SUPPLEMENT 1

mesalnicotinate <- biotran_both %>% filter(metab=="nicotinate")
mesalnicotinatefig<-mesalnicotinate %>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=binary154)) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
   scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #3x3 pdf

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesalnicotinatefig.png', plot=mesalnicotinatefig, width = 3, height = 3, dpi = 600) 


mesalnicotinur <- biotran_both %>% filter(metab=="nicotinuric.acid")
mesalnicotinurfig<-mesalnicotinur %>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=binary154)) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
   scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #3x3 pdf

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesalnicotinurfig.png', plot=mesalnicotinurfig, width = 3, height = 3, dpi = 600) 

mesalNAM <- biotran_both %>% filter(round=="X1.m")
mesalNAMfig<-mesalNAM %>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=binary154)) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
   scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #3x3 pdf

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesalNAMfig.png', plot=mesalNAMfig, width = 3, height = 3, dpi = 600)

########SUPPLEMENT 2

mesal242 <- biotranHMP2  %>% filter(round=="X242")
mesal242fig<-mesal242 %>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=as.factor(binary154))) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
   scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #3x3 pdf

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesal242fig.png', plot=mesal242fig, width = 3, height = 3, dpi = 600)

mesal373 <- biotranHMP2 %>% filter(metab=="X373.1254_0.71")
mesal373fig<-mesal373%>%
 ggplot(aes(x=as.factor(binary154),y=value+1,color=as.factor(binary154))) +
  geom_boxplot() + 
  theme_bw(14) +
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4) + 
   scale_y_log10() + xlab("")+theme(strip.background = element_blank(),legend.position = "none")+scale_color_manual(values=alpha(c("#3C5488FF","#DC0000FF"),0.8))+
scale_x_discrete(labels=c("non-user","user")) #3x3 pdf

ggsave(filename='/home/rsm34/5asabb3/submission/figures/mesal373fig.png', plot=mesal373fig, width = 3, height = 3, dpi = 600) 

