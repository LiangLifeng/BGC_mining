library(UpSetR)


group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")

bgc_class = c("RiPPs","NRPS","PKS","Terpene","Others")
bgc_colors =c ("#e78ac3","#8da0cb","#fc8d62","#66c2a5","#cccccc") 

bgc <- read.table("../data/GEMs-OMD-FWL-PS.all.antisMash.3Group.Merge_BGC.txt",header=T,sep = "	",check.names=F,row.names=1)
genomeInfo <- read.table("../data/3Group.represent.genome.info.all.split.txt.group.Freeliving.txt",header=T,sep = "	",check.names=F)

bgc[bgc > 1]  <- 1


genomeInfo <- subset(genomeInfo, (ID %in% rownames(bgc)))

bgc$ID <- rownames(bgc)
plot <-  merge(bgc,genomeInfo, by=c("ID"))


upset(plot)

bgc_class = c("RiPPs","NRPS","PKS","Terpene","Others")
bgc_colors =c ("#e78ac3","#8da0cb","#fc8d62","#66c2a5","#cccccc") 
group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")






#Fig 1a 
library("dplyr")
library("ggpubr")
df <- as.data.frame(matrix(nrow=3, ncol=0))

group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")
df$Group <- group_list
df$genomeNum <-c(8471,2379,2479)
plot <- df %>%
  summarize(Group = Group,
            genomeNum = genomeNum,
            Prevalence = (genomeNum /  sum(genomeNum))*100,
            labs = paste0(genomeNum, " (", round(Prevalence,2), "%)")
            )

plot$Group <- factor(plot$Group ,levels=group_list)

gn <- ggpubr::ggpie(plot, "Prevalence", label = labs,
      lab.pos = "in",
      fill = "Group", color = "white",
      lab.font = c(2, "bold", "black"),
      palette = color_var)+theme(legend.position="none")
gn
pdf("../data/3Group.genomeNum.pdf",width=1.5,height=1.5,onefile = FALSE)
gn
dev.off()
png("../data/3Group.genomeNum.png",width = 1.5,height =4.5,units='in',res=900)
gn
dev.off()


ASphylum <- read.table("../data/antisMash.phylum.exit.txt",header=T,sep = "	",check.names=F,row.names=1)
ASfamily <- read.table("../data/antisMash.family.exit.txt",header=T,sep = "	",check.names=F,row.names=1)

phylum <- read.table("../data/all.phylum.exit.txt",header=T,sep = "	",check.names=F,row.names=1)
family <- read.table("../data/all.family.exit.txt",header=T,sep = "	",check.names=F,row.names=1)



group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")
upset(phylum,order.by = "degree")

color_var = c("#dd5129","#4daf4a","#984ea3")


pu <- upset(phylum, sets = rev(group_list), sets.bar.color = rev(color_var) ,order.by = "degree",decreasing = TRUE,
      empty.intersections = "on",keep.order = TRUE,number.angles = 0, mainbar.y.label = "# Phylum", sets.x.label="# Phylum",
      queries = list(list(query = intersects, params = 'Marine', color = '#984ea3', active = F),
                     list(query = intersects, params = 'Freshwater', color = '#984ea3', active = F),
                     list(query = intersects, params = 'Terrestrial', color = '#984ea3', active = F),
                     list(query = intersects, params = list("Marine", "Freshwater"), color = "#4daf4a", active = F),
                     list(query = intersects, params = list("Marine", "Terrestrial"), color = "#4daf4a", active = F),
                     list(query = intersects, params = list("Terrestrial", "Freshwater"), color = "#4daf4a", active = F),
                     list(query = intersects, params = list("Terrestrial", "Freshwater","Marine"), color = "#dd5129", active = T)
                     
      )
      )
pu
pdf("../data/antismash.3Group.phylum.upset.pdf",width= 3,height=2.5,onefile = FALSE)
pu
dev.off()
png("../data/antismash.3Group.phylum.upset.png",width = 3,height =2.5,units='in',res=900)
pu
dev.off()

fu <- upset(family, sets = rev(group_list), sets.bar.color = rev(color_var) ,order.by = "degree",decreasing = TRUE,
            empty.intersections = "on",keep.order = TRUE,number.angles = 0, mainbar.y.label = "# Family", sets.x.label="# Family",
            queries = list(list(query = intersects, params = 'Marine', color = '#984ea3', active = F),
                           list(query = intersects, params = 'Freshwater', color = '#984ea3', active = F),
                           list(query = intersects, params = 'Terrestrial', color = '#984ea3', active = F),
                           list(query = intersects, params = list("Marine", "Freshwater"), color = "#4daf4a", active = F),
                           list(query = intersects, params = list("Marine", "Terrestrial"), color = "#4daf4a", active = F),
                           list(query = intersects, params = list("Terrestrial", "Freshwater"), color = "#4daf4a", active = F),
                           list(query = intersects, params = list("Terrestrial", "Freshwater","Marine"), color = "#dd5129", active = T)
                           
            )
            
            )
fu
pdf("../data/antismash.3Group.family.upset.pdf",width= 3.5,height=2.5,onefile = FALSE)
fu
dev.off()
png("../data/antismash.3Group.family.upset.png",width = 3.5,height =2.5,units='in',res=900)
fu
dev.off()



phylum <- read.table("../data/phylum.exit.otu.stat.txt",header=T,sep = "	",check.names=F)
library("dplyr")
library("ggpubr")
library("ggbreak")
library("ggplot2")


plot <- phylum %>%
  group_by(IDGroup,Type)%>%
  summarize(    total = sum(num))

color_var = c("#dd5129","#4daf4a","#984ea3")
group_list =c("Marine","Freshwater","Terrestrial")

nicelist <- c("3niches","2niches","1niches")
plot$Type <- factor(plot$Type ,levels=nicelist)
plot$IDGroup <- factor(plot$IDGroup ,levels=group_list)


write.table(plot,"../data/3Group.antisMash.3niches.txt" , sep = "\t"    , row.names = T, col.names = NA, quote = F)


  bar2 <- ggplot(plot,aes(IDGroup,total,fill=Type))+
    geom_bar(stat="identity", position=position_dodge(0.7),width=0.8)+
    
    #geom_text(label=plot$total, vjust=1.6, color="black", size=3.5,angle=90, position=position_dodge(0.7))+
    scale_fill_manual(values=color_var,labels=nicelist,limits = nicelist, breaks = nicelist)+
    labs(y="# genome",x="")+
    scale_y_continuous(expand = c(0,0))+
    
    #scale_y_break(breaks=c(1500,2500),scale=0.2,ticklabels = seq(2500,3000,200),)+  
    scale_y_break(breaks=c(2500,6000),scale=0.5,ticklabels = seq(6000,7000,1000))+
    theme_bw() +
    theme(axis.text.x = element_text(color="black",size=10,angle=-45,hjust= 0 ,vjust = 0 ,face="bold"),
          axis.text.y = element_text(color="black",size=6,face="bold"),
          axis.title.y=element_text(color="black",size=10),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(size = 8, colour = "black", face="bold"), 
          strip.background  = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 6, colour = "black"),
          legend.title = element_text(size = 6),
          legend.key.width = unit(0.3, 'cm'),
          legend.key.size = unit(0, 'lines'),
          panel.grid = element_blank(),
          panel.background = element_blank())
    
  bar2

  
  pdf("../data/3Group.antisMash.3niches.pdf",width=2.8,height=3,onefile = FALSE)
  bar2
  dev.off()
  png("../data/3Group.antisMash.3niches.png",width = 2.8,height =3,units='in',res=900)
  bar2
  dev.off()
  
  
  
  write.table(plot,"../data/3Group.all.3niches.txt" , sep = "\t"    , row.names = T, col.names = NA, quote = F)
