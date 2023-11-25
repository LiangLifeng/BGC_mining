library("ggpubr")
library("ggplot2")
library("reshape")
library("ggprism")
library("rstatix")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("aplot")
library("ggbreak")


# antismash -----------------------------------------
group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")

bgc_class = c("RiPPs","NRPS","PKS","Terpene","Others")
bgc_colors =c ("#e78ac3","#8da0cb","#fc8d62","#66c2a5","#cccccc") 

profile <- as.data.frame(t(read.table("../data/GEMs-OMD-FWL-PS.all.antisMash.3Group.plan3.Freeliving.xls.S",header=T, sep = "\t",check.name = F,row.names=1)))
classInfo <- read.table("../data/antisMash.Merge.BGC.class.txt",header=T,  sep = "\t", check.name = F,quote="")
totalBGC <- as.data.frame(t(read.table("../data/GEMs-OMD-FWL-PS.all.antisMash.3Group.plan3.Freeliving.total.xls",header=T, sep = "\t",check.name = F,row.names=1)))


totalBGC$ID <- rownames(totalBGC)
Ldata2 <-reshape2::melt(totalBGC,id.vars = c("ID"),variable.name = "Group", value.name = "total") 
Ldata2$newID <- paste(Ldata2$ID,Ldata2$Group,sep = "-")

needID <- intersect(classInfo$ID, rownames(profile))
classInfo <- subset(classInfo, (ID %in% needID))
profile$ID <- rownames(profile)
data1 <- merge(classInfo,profile, sort=F, by=c("ID"))
Ldata1 <-  reshape2::melt(data1,id.vars = c("ID", "Class"),variable.name = "Group", value.name = "abundance")
Ldata1$newID <- paste(Ldata1$ID,Ldata1$Group,sep = "-")


data <-  merge(Ldata1,Ldata2, sort=F, by=c("newID"))



plot <- data %>%
  group_by(ID.x)%>%
  summarize(Group = Group.x,
            Class = Class,
            abundance = abundance,
            total = total,
            Prevalence = (total /  sum(total))*100)


subGroup <- subset(plot, (Group %in% group_list[3]))

subGroup$Class <- factor(subGroup$Class ,levels=bgc_class)

orderPlot <- subGroup %>% 
  group_by(Class) %>% 
  arrange(Class,desc(Prevalence),ID.x)



plot$ID <- factor(plot$ID.x ,levels=orderPlot$ID.x)

theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   
    legend.text = element_text(color="black",size=10),
    legend.spacing.x=unit(0.1,'cm'), 
    legend.key.width=unit(0.5,'cm'), 
    legend.key.height=unit(0.5,'cm'), 
    legend.background=element_blank())
}

bar <- ggplot(plot,aes(ID,Prevalence,fill=Group),show.legend = FALSE)+
  geom_col(position="stack",width=0.8) +
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)+
  labs(y="% BGCs")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color="black",size=10),
        axis.title.y=element_text(color="black",size=10),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank(),
        legend.position="none")+
  theme_niwot()

bar

bar2 <- ggplot(plot,aes(ID,log10(total),fill=Group))+
  geom_bar(stat="identity", position=position_dodge(0.7),width=0.8)+
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)+
  labs(y="log10(# total BGC)")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y.left=element_text(color="black",size=10,),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right=element_blank(),
        axis.title.y=element_text(color="black",size=10,angle=90),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank(),
        legend.position="none")+
  theme_niwot()

bar2

classInfo$ID <- factor(classInfo$ID ,levels=orderPlot$ID.x)

group1 <- classInfo %>% select(1,2) %>% mutate(group="BGC Class") %>% 
  ggplot(aes(ID,group,fill=Class))+
  geom_tile()+
  scale_y_discrete(expand = c(0,0),position="right")+
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=bgc_colors,labels=bgc_class,limits = bgc_class, breaks = bgc_class)+
  theme_void()+
  theme(axis.text.x = element_text(color="black",size=8,angle=-90,hjust= 0 ,vjust = 0.5 ),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color="black",size=10),
        axis.title.y=element_text(color="black",size=10),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank(),
        legend.position="none")+
  theme_niwot()

group1


group2 <- classInfo %>% select(1,2) %>% mutate(group="BGC Class") %>% 
  
  ggplot(plot,aes(ID,Group,fill=Class))+
  geom_tile()+
  scale_y_discrete(expand = c(0,0),position="right")+
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=bgc_colors,labels=bgc_class,limits = bgc_class, breaks = bgc_class)+
  theme_void()+
  theme(axis.text.x = element_text(color="black",size=8,angle=-90,hjust= 0 ,vjust = 0.1 ),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color="black",size=10),
        axis.title.y=element_text(color="black",size=10),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank(),
        legend.position="none")+
  theme_niwot()

group2


mulitPlot <- bar2 %>% insert_bottom(bar,height = 0.5) %>% 
  insert_bottom(group1,height = 0.05)

mulitPlot

pdf("../data/BGC-3Group.Precentation.log.pdf",width=6.5,height=3.5,onefile = FALSE)
mulitPlot
dev.off()
png("../data/BGC-3Group.Precentation.log.png",width = 6.5,height =3.5,units='in',res=900)
mulitPlot
dev.off()

# antismash Ripps ----

subGroup <- subset(plot, (Class %in% c("RiPPs")))

bar <- ggplot(subGroup,aes(ID,Prevalence,fill=Group),show.legend = FALSE)+
  geom_col(position="stack",width=0.8) +
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)+
  labs(y="% BGCs")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_text(color="black",size=8,angle=-90,hjust= 0 ,vjust = 0.5 ),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color="black",size=10),
        axis.title.y=element_text(color="black",size=10),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank(),
        legend.position="none")+
  theme_niwot()

bar

bar2 <- ggplot(subGroup,aes(ID,log10(total),fill=Group))+
  geom_bar(stat="identity", position=position_dodge(0.7),width=0.8)+
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)+
  labs(y="log10(# total BGC)")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y.left=element_text(color="black",size=10,),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right=element_blank(),
        axis.title.y=element_text(color="black",size=10,angle=90),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank(),
        legend.position="none")+
  theme_niwot()

bar2

mulitPlot <- bar2 %>% insert_bottom(bar,height = 0.5) 

mulitPlot

pdf("../data/BGC-3Group.Ripps.Precentation.pdf",width=2.5,height=3,onefile = FALSE)
mulitPlot
dev.off()
png("../data/BGC-3Group.Ripps.Precentation.png",width = 2.5,height =3,units='in',res=900)
mulitPlot
dev.off()

# BAGEL4 ----

group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")

bgc_class = c("RiPPs","NRPS","PKS","Terpene","Others")
bgc_colors =c ("#e78ac3","#8da0cb","#fc8d62","#66c2a5","#cccccc") 

profile <- as.data.frame(t(read.table("../data/BAGEL4/GEMs-OMD-FWL-PS.BAGEL4.3Group.plan3.Freeliving.xls",header=T, sep = "\t",check.name = F,row.names=1)))
totalBGC <- as.data.frame(t(read.table("../data/BAGEL4/GEMs-OMD-FWL-PS.BAGEL4.3Group.plan3.Freeliving.total.xls",header=T, sep = "\t",check.name = F,row.names=1)))


totalBGC$ID <- rownames(totalBGC)
Ldata2 <-reshape2::melt(totalBGC,id.vars = c("ID"),variable.name = "Group", value.name = "total") 
Ldata2$newID <- paste(Ldata2$ID,Ldata2$Group,sep = "-")



profile$ID <- rownames(profile)
Ldata1 <-  reshape2::melt(profile,id.vars = c("ID"),variable.name = "Group", value.name = "abundance")
Ldata1$newID <- paste(Ldata1$ID,Ldata1$Group,sep = "-")


data <-  merge(Ldata1,Ldata2, sort=F, by=c("newID"))



plot <- data %>%
  group_by(ID.x)%>%
  summarize(Group = Group.x,
            abundance = abundance,
            total = total,
            Prevalence = (total /  sum(total))*100)


subGroup <- subset(plot, (Group %in% group_list[3]))

orderPlot <- subGroup %>% 
  arrange(desc(Prevalence),ID.x)



plot$ID <- factor(plot$ID.x ,levels=orderPlot$ID.x)

theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   
    legend.text = element_text(color="black",size=10), 
    legend.spacing.x=unit(0.1,'cm'), 
    legend.key.width=unit(0.5,'cm'), 
    legend.key.height=unit(0.5,'cm'), 
    legend.background=element_blank())
}

bar <- ggplot(plot,aes(ID,Prevalence,fill=Group),show.legend = FALSE)+
  geom_col(position="stack",width=0.8) +
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)+
  labs(y="% BGCs")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  theme(
        axis.text.x = element_text(color="black",size=8,angle=-90,hjust= 0 ,vjust = 0.5 ),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color="black",size=10),
        axis.title.y=element_text(color="black",size=10),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank(),
        legend.position="none")+
  theme_niwot()

bar

bar2 <- ggplot(plot,aes(ID,total,fill=Group))+
  geom_bar(stat="identity", position=position_dodge(0.7),width=0.8)+
  scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)+
  labs(y="# total BGC")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y.left=element_text(color="black",size=10,),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right=element_blank(),
        axis.title.y=element_text(color="black",size=10,angle=90),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank(),
        legend.position="none")+
  theme_niwot()

bar2

mulitPlot <- bar2 %>% insert_bottom(bar,height = 0.5)

mulitPlot

pdf("../data/BAGEL4/BAGEL4.BGC-3Group.Precentation.total.pdf",width = 3,height=3,onefile = FALSE)
mulitPlot
dev.off()
png("../data/BAGEL4/BAGEL4.BGC-3Group.Precentation.total.png",width = 3,height =3,units='in',res=900)
mulitPlot
dev.off()
