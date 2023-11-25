library("ggpubr")
library("ggplot2")
library("reshape")
library("ggprism")
library("rstatix")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("aplot")




group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")

bgc_class = c("RiPPs","NRPS","PKS","Terpene","Others")
bgc_colors =c ("#e78ac3","#8da0cb","#fc8d62","#66c2a5","#cccccc") 

family <- read.table("../data//Family/Alpha_div/All_samples.alpha_diversity.xls",header=T, sep = "\t",check.name = F,row.names=1)
phylum <- read.table("../data//Phylum/Alpha_div/All_samples.alpha_diversity.xls",header=T, sep = "\t",check.name = F,row.names=1)

family$Group <- rownames(family)
phylum$Group <- rownames(phylum)

phylum$Group <- factor(phylum$Group,levels=group_list)
family$Group <- factor(family$Group,levels=group_list)

bar1 <- ggbarplot(family, x = "Group", y = "Shannon.index", fill = "Group",position = position_dodge(0.75),alpha=0.7,
                  ggtheme = theme_bw() +
                    theme(axis.text.x = element_text(color="black",size=0,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                          axis.text.y = element_text(color="black",size=10,face="bold"),
                          axis.title.y=element_text(color="black",size=10),
                          axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          strip.text.x = element_text(size = 8, colour = "black", face="bold"), 
                          strip.background  = element_blank(),
                          legend.position = "none",
                          legend.text = element_text(size = 6, colour = "black"),
                          legend.title = element_text(size = 6),
                          legend.key.width = unit(0.3, 'cm'),
                          legend.key.size = unit(0, 'lines'),
                          panel.grid = element_blank(),
                          panel.background = element_blank()),
                  legend = "none",title = "",xlab = '', ylab = 'Family Shannon index',width = 0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))+
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)
bar1


pdf("../data//Family-Shannon-index.pdf",width=1.2,height=2,onefile = FALSE)
bar1
dev.off()
png("../data//Family-Shannon-index.png",width = 1.2,height = 2,units='in',res=600)
bar1
dev.off()
#phylum
bar2 <- ggbarplot(phylum, x = "Group", y = "Shannon.index", fill = "Group",position = position_dodge(0.75),alpha=0.7,
                  ggtheme = theme_bw() +
                    theme(axis.text.x = element_text(color="black",size=0,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                          axis.text.y = element_text(color="black",size=10,face="bold"),
                          axis.title.y=element_text(color="black",size=10),
                          axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          strip.text.x = element_text(size = 8, colour = "black", face="bold"), 
                          strip.background  = element_blank(),
                          legend.position = "none",
                          legend.text = element_text(size = 6, colour = "black"),
                          legend.title = element_text(size = 6),
                          legend.key.width = unit(0.3, 'cm'),
                          legend.key.size = unit(0, 'lines'),
                          panel.grid = element_blank(),
                          panel.background = element_blank()),
                  legend = "none",title = "",xlab = '', ylab = 'Phylum Shannon index',width = 0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))+
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)
bar2



pdf("../data//Phylum-Shannon-index.pdf",width=1.2,height=2,onefile = FALSE)
bar2
dev.off()
png("../data//Phylum-Shannon-index.png",width = 1.2,height = 2,units='in',res=600)
bar2
dev.off()

#number 

bar3 <- ggbarplot(family, x = "Group", y = "Family.number", fill = "Group",position = position_dodge(0.75),alpha=0.7,
                  ggtheme = theme_bw() +
                    theme(axis.text.x = element_text(color="black",size=0,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                          axis.text.y = element_text(color="black",size=10,face="bold"),
                          axis.title.y=element_text(color="black",size=10),
                          axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          strip.text.x = element_text(size = 8, colour = "black", face="bold"), 
                          strip.background  = element_blank(),
                          legend.position = "none",
                          legend.text = element_text(size = 6, colour = "black"),
                          legend.title = element_text(size = 6),
                          legend.key.width = unit(0.3, 'cm'),
                          legend.key.size = unit(0, 'lines'),
                          panel.grid = element_blank(),
                          panel.background = element_blank()),
                  legend = "none",title = "",xlab = '', ylab = '# Family',width = 0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))+
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)



pdf("../data//Family-number.pdf",width=1.2,height=2,onefile = FALSE)
bar3
dev.off()
png("../data//Family-number.png",width = 1.2,height = 2,units='in',res=600)
bar3
dev.off()
#phylum

phylum <- read.table("../data//upset/antisMash.phylum.exit.txt",header=T,sep = "	",check.names=F,row.names=1)
family <- read.table("../data//upset/antisMash.family.exit.txt",header=T,sep = "	",check.names=F,row.names=1)


Ldata2 <-reshape2::melt(family,variable.name = "Group", value.name = "total") 
plot <- Ldata2 %>%
  group_by(Group)%>%
  summarize(total = sum(total))

bar4 <- ggbarplot(plot, x = "Group", y = "total", fill = "Group",position = position_dodge(0.75),alpha=0.7,
                  ggtheme = theme_bw() +
                    theme(axis.text.x = element_text(color="black",size=0,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                          axis.text.y = element_text(color="black",size=10,face="bold"),
                          axis.title.y=element_text(color="black",size=10),
                          axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          strip.text.x = element_text(size = 8, colour = "black", face="bold"), 
                          strip.background  = element_blank(),
                          legend.position = "none",
                          legend.text = element_text(size = 6, colour = "black"),
                          legend.title = element_text(size = 6),
                          legend.key.width = unit(0.3, 'cm'),
                          legend.key.size = unit(0, 'lines'),
                          panel.grid = element_blank(),
                          panel.background = element_blank()),
                  legend = "none",title = "",xlab = '', ylab = '# Family',width = 0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))+
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)




pdf("../data//family-number.antismash.pdf",width=1.2,height=2,onefile = FALSE)
bar4
dev.off()
png("../data//family-number.antismash.png",width = 1.2,height = 2,units='in',res=600)
bar4
dev.off()

