library("dplyr")
library("reshape")
library("ggplot2")
library("ggpubr")

group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")

#phylum Class   Order   Family  Genus   Species
phylum <- read.table("../data/3Group.represent.genome.info.all.split.txt.group.phylum.stat.xls",header=T,sep = "	",check.names=F)
Class <- read.table("../data/3Group.represent.genome.info.all.split.txt.group.Class.stat.xls",header=T,sep = "	",check.names=F)
Order <- read.table("../data/3Group.represent.genome.info.all.split.txt.group.Order.stat.xls",header=T,sep = "	",check.names=F)
Family <- read.table("../data/3Group.represent.genome.info.all.split.txt.group.Family.stat.xls",header=T,sep = "	",check.names=F)
Genus <- read.table("../data/3Group.represent.genome.info.all.split.txt.group.Genus.stat.xls",header=T,sep = "	",check.names=F)
Species <- read.table("../data/3Group.represent.genome.info.all.split.txt.group.Species.stat.xls",header=T,sep = "	",check.names=F)



phylumStat <- phylum %>%
  group_by(Group)%>%
  summarize(level = level,num = num,Prevalence = (num /  sum(num))*100 )%>%
  mutate(taxonomy=if_else(level == "p__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  summarize(Phylum = sum(Prevalence) )

ClassStat <- Class %>%
  group_by(Group)%>%
  summarize(level = level,num = num,Prevalence = (num /  sum(num))*100 )%>%
  mutate(taxonomy=if_else(level == "c__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  summarize(Class = sum(Prevalence) )

OrderStat <- Order %>%
  group_by(Group)%>%
  summarize(level = level,num = num,Prevalence = (num /  sum(num))*100 )%>%
  mutate(taxonomy=if_else(level == "o__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  summarize(Order = sum(Prevalence) )

FamilyStat <- Family %>%
  group_by(Group)%>%
  summarize(level = level,num = num,Prevalence = (num /  sum(num))*100 )%>%
  mutate(taxonomy=if_else(level == "f__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  summarize(Family = sum(Prevalence) )

GenusStat <- Genus %>%
  group_by(Group)%>%
  summarize(level = level,num = num,Prevalence = (num /  sum(num))*100 )%>%
  mutate(taxonomy=if_else(level == "g__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  summarize(Genus = sum(Prevalence) )

SpeciesStat <- Species %>%
  group_by(Group)%>%
  summarize(level = level,num = num,Prevalence = (num /  sum(num))*100 )%>%
  mutate(taxonomy=if_else(level == "s__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  summarize(Species = sum(Prevalence) )

P1 <- melt(phylumStat)
C1 <- melt(ClassStat)
O1 <- melt(OrderStat)
F1 <- melt(FamilyStat)
G1 <- melt(GenusStat)
S1 <- melt(SpeciesStat)

plot <- rbind(P1,C1,O1,F1,G1,S1)
plot$newID <- paste0(plot$variable, ":", plot$taxonomy)
plot$newGroup <- paste0(plot$Group, ":", plot$variable)
plot$newColor <- paste0(plot$Group, ":", plot$taxonomy)

group_list =c("Marine","Freshwater","Terrestrial")
taxoLevel  <- c("Phylum", "Class"  , "Order" ,  "Family" , "Genus" ,  "Species")
color_var = c(rep(c("#bdbdbd","#0f7ba2","#bdbdbd","#43b284","#bdbdbd","#fab255"),6))

newGroupOrder <- c()

  for (t in taxoLevel ){
    for (g in group_list) {
      newGroupOrder <- c(newGroupOrder,paste0(g,":",t))
      }
  }
variableOrder <-c()
  for (gg in group_list ){
    for (cc in c("Unclassified","Classified")) {
      variableOrder <- c(variableOrder,paste0(gg,":",cc))
    }
  }


newGroupOrder
variableOrder
plot$newGroup <- factor(plot$newGroup ,levels=newGroupOrder)
plot$newColor <- factor(plot$newColor ,levels=variableOrder)
plot$Group <- factor(plot$Group ,levels=group_list)

#先得到柱形图
p <- ggplot(plot, aes(Group, value, fill = newColor)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_fill_manual(values = c(color_var)) +
  labs(x = '', y = '', fill = 'taxonomy classified %') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank())+
  facet_grid(.~ variable)
p
#再经坐标系变换得到圆环图
p + coord_polar(theta = 'y')


p1 <- ggplot(plot, aes(newGroup, value, fill = newColor)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_fill_manual(values = c(color_var)) +
  labs(x = '', y = '', fill = 'taxonomy classified %') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank())

p1
#再经坐标系变换得到圆环图
p1 + coord_polar(theta = 'y')




p <- ggbarplot(plot, "Group", "value",fill = "newColor", color = "newColor", 
          ggtheme = theme_bw() +
            theme(axis.text.x = element_text(color="black",size=0,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                  axis.text.y = element_text(color="black",size=10,face="bold"),
                  axis.title.y=element_text(color="black",size=10,face="bold"),
                  axis.line = element_line(color="black"),
                  axis.ticks = element_line(color="black"),
                  strip.text.x = element_text(size = 10, colour = "black", face="bold"), 
                  strip.background  = element_blank(),
                  legend.position = "none",
                  legend.text = element_text(size = 10, colour = "black"),
                  legend.title = element_text(size = 10),
                  legend.key.width = unit(0.3, 'cm'),
                  legend.key.size = unit(0, 'lines'),
                  panel.grid = element_blank(),
                  panel.background = element_blank()),
          legend = "right",title = "",xlab = '', ylab = 'taxonomy classified %',width = 0.7)+ 
  scale_fill_manual(values=color_var,labels=variableOrder,limits = variableOrder, breaks = variableOrder) +
  scale_color_manual(values=color_var,labels=variableOrder,limits = variableOrder, breaks = variableOrder) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank())
pp <- facet(p, facet.by = "variable",nrow =1,short.panel.labs = TRUE,strip.position = "top", 
            panel.labs.font = list(size = 10, angle = 0),panel.grid = element_blank(), panel.labs.background = list(fill = NA, color = NA))

pdf("../data/3Group.taxonomy.classified.rate.pdf",width= 6,height=3,onefile = FALSE)
pp
dev.off()
png("../data/3Group.taxonomy.classified.rate.png",width = 6,height =3,units='in',res=900)
pp
dev.off()

# level number -----------------------------------------
PhylumNum <- phylum %>%
  group_by(Group) %>%
  dplyr::summarize(Phylum = n() )

ClassNum <- Class %>%
  group_by(Group) %>%
  dplyr::summarize(Class = n() )

OrderNum <- Order %>%
  group_by(Group) %>%
  dplyr::summarize(Order = n() )

FamilyNum <- Family %>%
  group_by(Group) %>%
  dplyr::summarize(Family = n() )

GenusNum <- Genus %>%
  group_by(Group) %>%
  dplyr::summarize(Genus = n() )

SpeciesNum <- Species %>%
  group_by(Group) %>%
  dplyr::summarize(Species = n() )

P1 <- melt(PhylumNum)
C1 <- melt(ClassNum)
O1 <- melt(OrderNum)
F1 <- melt(FamilyNum)
G1 <- melt(GenusNum)
S1 <- melt(SpeciesNum)

plot <- rbind(P1,C1,O1,F1,G1,S1)

plot$Group <- factor(plot$Group ,levels=group_list)

group_list =c("Marine","Freshwater","Terrestrial")
color_var1 = c("#0f7ba2","#43b284","#fab255")

p2 <- ggbarplot(plot, "Group", "value",fill = "Group", color = "Group", 
                ggtheme = theme_bw() +
                  theme(axis.text.x = element_text(color="black",size=0,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                        axis.text.y = element_text(color="black",size=10,face="bold"),
                        axis.title.y=element_text(color="black",size=10,face="bold"),
                        axis.line = element_line(color="black"),
                        axis.ticks = element_line(color="black"),
                        strip.text.x = element_text(size = 10, colour = "black", face="bold"), 
                        strip.background  = element_blank(),
                        legend.position = "none",
                        legend.text = element_text(size = 10, colour = "black"),
                        legend.title = element_text(size = 10),
                        legend.key.width = unit(0.3, 'cm'),
                        legend.key.size = unit(0, 'lines'),
                        panel.grid = element_blank(),
                        panel.background = element_blank()),
                legend = "none",title = "",xlab = '', ylab = 'numbers',width = 0.7)+ 
  scale_fill_manual(values=color_var1,labels=group_list,limits = group_list, breaks = group_list) +
  scale_color_manual(values=color_var1,labels=group_list,limits = group_list, breaks = group_list) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank())
pp2 <- facet(p2, facet.by = "variable",nrow =1,short.panel.labs = TRUE,strip.position = "top", scales = "free_y",
             panel.labs.font = list(size = 10, angle = 0),panel.grid = element_blank(), panel.labs.background = list(fill = NA, color = NA))

pp2
pdf("../data/3Group.taxonomy.level.numbers.pdf",width= 6,height=3,onefile = FALSE)
pp2
dev.off()
png("../data/3Group.taxonomy.level.numbers.png",width = 6,height =3,units='in',res=900)
pp3
dev.off()





# genome number -----------------------------------------
PhylumNum <- phylum %>%
  mutate(taxonomy=if_else(level == "p__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  dplyr::summarize(Phylum = sum(num) )

ClassNum <- Class %>%
  mutate(taxonomy=if_else(level == "c__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  dplyr::summarize( Class = sum(num) )

OrderNum <- Order %>%
  mutate(taxonomy=if_else(level == "o__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  dplyr::summarize( Order = sum(num) )

FamilyNum <- Family %>%
  mutate(taxonomy=if_else(level == "f__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  dplyr::summarize( Family = sum(num) )

GenusNum <- Genus %>%
  mutate(taxonomy=if_else(level == "g__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  dplyr::summarize( Genus = sum(num) )

SpeciesNum <- Species %>%
  mutate(taxonomy=if_else(level == "s__","Unclassified","Classified")) %>%
  group_by(Group,taxonomy) %>%
  dplyr::summarize( Species = sum(num) )

P1 <- melt(PhylumNum)
C1 <- melt(ClassNum)
O1 <- melt(OrderNum)
F1 <- melt(FamilyNum)
G1 <- melt(GenusNum)
S1 <- melt(SpeciesNum)

plot <- rbind(P1,C1,O1,F1,G1,S1)
plot$newID <- paste0(plot$variable, ":", plot$taxonomy)
plot$newGroup <- paste0(plot$Group, ":", plot$variable)
plot$newColor <- paste0(plot$Group, ":", plot$taxonomy)

plot$newGroup <- factor(plot$newGroup ,levels=newGroupOrder)
plot$newColor <- factor(plot$newColor ,levels=variableOrder)
plot$Group <- factor(plot$Group ,levels=group_list)

p3 <- ggbarplot(plot, "Group", "value",fill = "newColor", color = "newColor", 
               ggtheme = theme_bw() +
                 theme(axis.text.x = element_text(color="black",size=0,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                       axis.text.y = element_text(color="black",size=10,face="bold"),
                       axis.title.y=element_text(color="black",size=10,face="bold"),
                       axis.line = element_line(color="black"),
                       axis.ticks = element_line(color="black"),
                       strip.text.x = element_text(size = 10, colour = "black", face="bold"), 
                       strip.background  = element_blank(),
                       legend.position = "none",
                       legend.text = element_text(size = 10, colour = "black"),
                       legend.title = element_text(size = 10),
                       legend.key.width = unit(0.3, 'cm'),
                       legend.key.size = unit(0, 'lines'),
                       panel.grid = element_blank(),
                       panel.background = element_blank()),
               legend = "none",title = "",xlab = '', ylab = 'numbers',width = 0.7)+ 
  scale_fill_manual(values=color_var,labels=variableOrder,limits = variableOrder, breaks = variableOrder) +
  scale_color_manual(values=color_var,labels=variableOrder,limits = variableOrder, breaks = variableOrder) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank())
pp3 <- facet(p3, facet.by = "variable",nrow =1,short.panel.labs = TRUE,strip.position = "top", scales = "free_y",
            panel.labs.font = list(size = 10, angle = 0),panel.grid = element_blank(), panel.labs.background = list(fill = NA, color = NA))

pp3
pdf("../data/3Group.taxonomy.genome.numbers.pdf",width= 6,height=3,onefile = FALSE)
pp3
dev.off()
png("../data/3Group.taxonomy.genome.numbers.png",width = 6,height =3,units='in',res=900)
pp3
dev.off()
