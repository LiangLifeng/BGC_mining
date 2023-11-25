library("ggtree")
library("ggtreeExtra")
library("ggplot2")
library("MicrobiotaProcess")
library("tidytree")
library("ggstar")
library("forcats")
library("RColorBrewer")
library("MetBrewer")
library("ape")
library("phylobase")
library('cowplot')
library("treeio")
library("ggnewscale")
library("tidyverse")
library("dplyr")
library("ggpubr")
library("reshape2")


sumRow_color <- function(onedata,level){
  dd <- as.data.frame(rowSums(onedata))
  onedataSum <- as.data.frame(t(matrix(unlist(strsplit(rownames(dd),":")),nrow=2)))
  colnames(onedataSum) <- c("Group","taxonmy")
  onedataSum$bgcTotal <- dd$`rowSums(onedata)`

  onedataSum <- onedataSum %>%group_by(taxonmy)%>%
                summarize(alltotal = sum(bgcTotal) , 
                          numGroup = ifelse(alltotal >1000, ">1000",
                                            ifelse(alltotal>500,"500-1000",
                                                   ifelse(alltotal>100,"100-500",
                                                          ifelse(alltotal>10,"10-100","<10")))))
  rownames(onedataSum) <- onedataSum$taxonmy
  return(onedataSum)
}


# input ----
# parsing the output of MetaPhlAn
mpse3 <- mp_import_metaphlan(profile=".//bac120.gtdb_representative.family.most.txt", mapfilename=".//sample_info.txt")
# The result is stored to the taxatree or otutree slot, you can use mp_extract_tree to extract the specific slot.
taxa.tree <- mpse3 %>% 
  mp_extract_tree(type="taxatree")

bgPhylum <- read.table(".//data/phylum.marker.txt.2",header=T, sep = "\t", check.name = F,quote="")
Phylumlabel <- read.table(".//data/phylum.marker.txt.3.txt",header=T, sep = "\t", check.name = F,quote="")

## Family ----
antisMashFamilyGN <- read.table(".//GEMs-OMD-FWL-PS.all.antisMash.3Group.plan3.Family.anno.genoeme.number.xls",header=T,  sep = "\t", check.name = F,quote="",row.names=1)
antisMashFamilyGNSum <- sumRow_color(antisMashFamilyGN,"Family")
antisMashFamilyGNSum$numGroup
#"<10","10-100","100-500","500-1000",">1000"
showFamilyID <-  subset(antisMashFamilyGNSum, (numGroup %in% c("10-100","100-500","500-1000",">1000")))$taxonmy
showFamilyID
MergeBGC <- read.table(".//data/GEMs-OMD-FWL-PS.all.antisMash.3Group.Merge_BGC.plan3.Family.plot.txt",header=T,  sep = "\t", check.name = F,quote="")

totalBGCAbundance <- read.table(".//data/GEMs-OMD-FWL-PS.all.antisMash.3Group.Merge_BGC.plan3.Family.groupMean.Level.xls",header=T,  sep = "\t", check.name = F,quote="",row.names=1)

totalgenomeNum <- read.table(".//data/GEMs-OMD-FWL-PS.all.antisMash.3Group.Merge_BGC.plan3.Family.genome.number.xls",header=T,  sep = "\t", check.name = F,quote="")



RiPPs <- read.table(".//data/RiPPs.plot.matrix.txt",header=T,  sep = "\t", check.name = F,quote="")
NRPS <- read.table(".//data/NRPS.plot.matrix.txt",header=T,  sep = "\t", check.name = F,quote="")
PKS <- read.table(".//data/PKS.plot.matrix.txt",header=T,  sep = "\t", check.name = F,quote="")
Terpene <- read.table(".//data/Terpene.plot.matrix.txt",header=T,  sep = "\t", check.name = F,quote="")
Others <- read.table(".//data/Others.plot.matrix.txt",header=T,  sep = "\t", check.name = F,quote="")




## Order ----
antisMashOrderGN <- read.table(".//GEMs-OMD-FWL-PS.all.antisMash.3Group.plan3.Order.anno.genoeme.number.xls",header=T,  sep = "\t", check.name = F,quote="",row.names=1)
antisMashOrderGNSum <- sumRow_color(antisMashOrderGN,"Order")

## Class ----
antisMashClassGN <- read.table(".//GEMs-OMD-FWL-PS.all.antisMash.3Group.plan3.Class.anno.genoeme.number.xls",header=T,  sep = "\t", check.name = F,quote="",row.names=1)
antisMashClassGNSum <- sumRow_color(antisMashClassGN,"Class")

## Phylum ----
antisMashPhylumGN <- read.table(".//GEMs-OMD-FWL-PS.all.antisMash.3Group.plan3.Phylum.anno.genoeme.number.xls",header=T,  sep = "\t", check.name = F,quote="",row.names=1)
antisMashPhylumGNSum <- sumRow_color(antisMashPhylumGN,"Phylum")

taxonmyGNSum <- rbind(antisMashFamilyGNSum,antisMashOrderGNSum,antisMashClassGNSum,antisMashPhylumGNSum)

#  color info----
colorlist1 <-c("#d9d9d9","#95c36e","#74c8c3",   "#5a97c1", "#295384", "#0a2e57")
totalBGCGroup <- c("NoUse","<10","10-100","100-500","500-1000",">1000")
bgc_class = c("RiPPs","NRPS","PKS","Terpene","Others")
bgc_colors =c ("#e78ac3","#8da0cb","#fc8d62","#66c2a5","#cccccc") 
group_list =c("Marine","Freshwater","Terrestrial")
#color_var = c("#081d58","#004529","#662506")
color_var = c("#0f7ba2","#43b284","#fab255")

phylumcolorVar<- c(
  "Proteobacteria"="#e41a1c",
  "Bacteroidota"="#0a2d46",
  "Actinobacteriota"="#4daf4a",
  "Verrucomicrobiota"="#984ea3",
  "Chloroflexota"="#ff7f00",
  "Planctomycetota"="#1c9d7c",
  "Patescibacteria"="#a65628",
  "Acidobacteriota"="#f781bf",
  "Marinisomatota"="#8dd3c7",
  "Cyanobacteria"="#80b1d3",
  "Desulfobacterota"="#bebada",
  "Firmicutes"="#fb8072",
  "Firmicutes_A"="#735852",
  "Myxococcota"="#fdb462",
  "Omnitrophota"="#b3de69",
  "Bdellovibrionota"="#fccde5",
  "Gemmatimonadota"="#bc80bd",
  "Nitrospirota"="#ccebc5",
  "Spirochaetota"="#ffed6f",
  "Desulfobacterota_F"="#bebada",
  "Desulfobacterota_D"="#bebada",
  "Desulfobacterota_B"="#bebada",
  "Firmicutes_B"="#fb8072",
  "Firmicutes_D"="#fb8072",
  "Bdellovibrionota_C"="#fccde5",
  "Desulfobacterota_E"="#bebada",
  "Firmicutes_C"="#fb8072",
  "Desulfobacterota_G"="#bebada",
  "Firmicutes_E"="#fb8072",
  "Firmicutes_H"="#fb8072",
  "Firmicutes_G"="#fb8072",
  "Firmicutes_F"="#fb8072"

)


tree0 <- as_tibble(taxa.tree)
tree1 <- as_tibble(taxa.tree) %>%  
        mutate(numGroup = ifelse(label %in% taxonmyGNSum$taxonmy , taxonmyGNSum[label,]$numGroup,"NoUse")) 

write.table (tree1, ".//tree.info.txt", sep = "\t"    , row.names = T, col.names = NA, quote = F)



# Layer 0: tree ----
#fan radial as.phylo  as.treedata
p0<- ggtree(as.treedata(tree1),aes(color = numGroup),layout="fan",size = 0.3,alpha = 1) + 
  scale_color_manual(values=colorlist1,labels=totalBGCGroup,limits = totalBGCGroup, breaks =totalBGCGroup, name = "#annoGenomes",
                     guide=guide_legend(keywidth=0.5, keyheight=0.5, order=1))+
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0,0,0,0))
p0

#p2 <- collapse(p0, node=2192)  + geom_point2(aes(subset=(node==c(2192))), shape=21, size=1.5, color="#737373", fill = "#737373")
#p3 <-   collapse(p2, node=2197)+ geom_point2(aes(subset=(node==2197)), shape=21, size=1.5, color="#737373", fill = "#737373") 
#p4 <-   collapse(p3, node=2418)+ geom_point2(aes(subset=(node==2418)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p5 <-   collapse(p4, node=2267)+ geom_point2(aes(subset=(node==2267)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p6 <-   collapse(p5, node=2244)+ geom_point2(aes(subset=(node==2244)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p7 <-   collapse(p6, node=2206)+ geom_point2(aes(subset=(node==2206)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p8 <-   collapse(p7, node=2209)+ geom_point2(aes(subset=(node==2209)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p9 <-   collapse(p8, node=3107)+ geom_point2(aes(subset=(node==3107)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p10 <-  collapse(p9, node=3098)+ geom_point2(aes(subset=(node==3098)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p11 <-  collapse(p10,node=3111) + geom_point2(aes(subset=(node==3111)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p12 <-    collapse(p11,node=2214) + geom_point2(aes(subset=(node==2214)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p13 <-    collapse(p12,node=2716) + geom_point2(aes(subset=(node==2716)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p14 <-    collapse(p13,node=2696) + geom_point2(aes(subset=(node==2696)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p15 <-    collapse(p14,node=2695) + geom_point2(aes(subset=(node==2695)), shape=21, size=1.5, color="#737373", fill = "#737373")
#p00 <-    collapse(p15,node=2725) + geom_point2(aes(subset=(node==2725)), shape=21, size=1.5, color="#737373", fill = "#737373")



#p160 <- p00+ new_scale_fill() + new_scale_color()


p161 <- open_tree(p0,5)+
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0,0,0,0))
p16 <- rotate_tree(p161, 90)

p16
# Layer 1: Product Classes -------------------------Others--------------------------------------

  OthersPD.1 <- subset(MergeBGC, select = c( "Level","Others")) 
  OthersPD <- subset(OthersPD.1, (Level %in% showFamilyID))
  OthersPD.S = filter(OthersPD, Others > 0)
  OthersPD.S$variableNum <- rep(0, times=length(OthersPD.S$Level))
  colnames(OthersPD.S) <- c("Level"  ,     "bgcClass"  ,    "variableNum")

  p17  <- p16 + new_scale_fill() + new_scale_color()
  p18 <- p17 +new_scale_fill() +  #ggnewscale::new_scale("fill") +
    geom_fruit(data=OthersPD.S,geom=geom_point,offset=0.05,mapping=aes(x= variableNum,y=Level, size = bgcClass) , shape = 21, fill ="#cccccc",color = "#cccccc" ,#
               pwidth=0.5,alpha=1, 
               axis.params=list(axis="y", text.angle=-45, hjust=0 ),
               grid.params=list()) + 
    scale_size_continuous(range = c(0, 5),name="Total BGC couunts of each Class",guide=guide_legend(keywidth=0.5, keyheight=0.5,override.aes=list(starshape=15),order=2))
  
  p18

  Others.1 <- melt(Others)
  Othersmap <- subset(Others.1, (id %in% showFamilyID))
  Othersmap$variable <- factor(Othersmap$variable ,levels=group_list)
  OthersmapF = filter(Othersmap, value > 0)
  
  p19  <- p18 + ggnewscale::new_scale("size") 
  p20 <- p19 + ggnewscale::new_scale("fill") +
    geom_fruit(data=OthersmapF,geom=geom_tile,offset=0.01,mapping=aes(y=id, x=variable, alpha=value, fill=variable),
               pwidth=0.1,
               axis.params=list(axis="none", text.angle=-45, hjust=0)) +
    scale_fill_manual(values=color_var,name="niches",guide=guide_legend(keywidth=0.5, keyheight=0.5, order=4)) +
    scale_alpha_continuous(range=c(0.5, 1),name="abundance", guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3))
    
  p20

# Layer 3: Product Classes -------------------------Terpene--------------------------------------
  TerpenePD.1 <- subset(MergeBGC, select = c( "Level","Terpene")) 
  TerpenePD <- subset(TerpenePD.1, (Level %in% showFamilyID))
  TerpenePD.S = filter(TerpenePD, Terpene > 0)
  TerpenePD.S$variableNum <- rep(0, times=length(TerpenePD.S$Level))
  colnames(TerpenePD.S) <- c("Level"  ,     "bgcClass"  ,    "variableNum")
  
  p21  <- p20 + new_scale_fill() + new_scale_color()  
  p22 <- p21 + new_scale_fill() +
    geom_fruit(data=TerpenePD.S,geom=geom_point,offset=0.07,mapping=aes(x= variableNum,y=Level , size = bgcClass) , shape = 21,fill ="#66c2a5",color = "#66c2a5" ,
               pwidth=0.5,alpha=1,
               axis.params=list(axis="none", text.angle=-45, hjust=0 ),
               grid.params=list() )  +
    scale_size_continuous(range = c(0, 5) ,guide="none")
  
  p22
  
  Terpene.1 <- melt(Terpene)
  Terpenemap <- subset(Terpene.1, (id %in% showFamilyID))
  Terpenemap$variable <- factor(Terpenemap$variable ,levels=group_list)
  TerpenemapF = filter(Terpenemap, value > 0)
  
  p23  <- p22 + new_scale_fill() + new_scale_color()
  p24 <- p23 +ggnewscale::new_scale("alpha") +
    geom_fruit(data=TerpenemapF,geom=geom_tile,offset=0.02,mapping=aes(y=id, x=variable, alpha=value, fill=variable),
               pwidth=0.1,
               axis.params=list(axis="none", text.angle=-45, hjust=0 )) +
    scale_fill_manual(values=color_var,guide="none") +
    scale_alpha_continuous(range=c(0.5, 1),name="Terpene abundance", guide=guide_legend(keywidth=0.5, keyheight=0.5, order=5))
  p24

# Layer 5: Product Classes -------------------------PKS--------------------------------------
  PKSPD.1 <- subset(MergeBGC, select = c( "Level","PKS")) 
  PKSPD <- subset(PKSPD.1, (Level %in% showFamilyID))
  PKSPD.S = filter(PKSPD, PKS > 0)
  PKSPD.S$variableNum <- rep(0, times=length(PKSPD.S$Level))
  colnames(PKSPD.S) <- c("Level"  ,     "bgcClass"  ,    "variableNum")
  
  p25  <- p24 + new_scale_fill() + new_scale_color()
  p26 <- p25 + new_scale_fill() +
    geom_fruit(data=PKSPD.S,geom=geom_point,offset=0.07,mapping=aes(x= variableNum,y=Level , size = bgcClass) , shape = 21,fill ="#fc8d62",color = "#fc8d62" ,
               pwidth=0.5,alpha=1,stat="identity",orientation="y", 
               axis.params=list(axis="none", text.angle=-45, hjust=0 ),
               grid.params=list() )  
  p26
  
  
  PKS.1 <- melt(PKS)
  PKSmap <- subset(PKS.1, (id %in% showFamilyID))
  PKSmap$variable <- factor(PKSmap$variable ,levels=group_list)
  PKSmapF = filter(PKSmap, value > 0)
  
  p27  <- p26 + new_scale_fill() + new_scale_color()
  p28 <- p27 + new_scale_fill() +ggnewscale::new_scale("alpha") +
    geom_fruit(data=PKSmapF,geom=geom_tile,offset=0.02,mapping=aes(y=id, x=variable, alpha=value, fill=variable),
               pwidth=0.1,
               axis.params=list(axis="none", text.angle=-45, hjust=0 )) +
    scale_fill_manual(values=color_var,guide="none") +
    scale_alpha_continuous(range=c(0.5, 1),name="PKS abundance", guide=guide_legend(keywidth=0.5, keyheight=0.5, order=7))
  p28

# Layer 7: Product Classes -------------------------NRPS--------------------------------------
  NRPSPD.1 <- subset(MergeBGC, select = c( "Level","NRPS")) 
  NRPSPD <- subset(NRPSPD.1, (Level %in% showFamilyID))
  NRPSPD.S = filter(NRPSPD, NRPS > 0)
  NRPSPD.S$variableNum <- rep(0, times=length(NRPSPD.S$Level))
  colnames(NRPSPD.S) <- c("Level"  ,     "bgcClass"  ,    "variableNum")
  
  p29  <- p28 + new_scale_fill() + new_scale_color()
  p30 <- p29 + new_scale_fill() +
    geom_fruit(data=NRPSPD.S,geom=geom_point,offset=0.07,mapping=aes(x= variableNum,y=Level , size = bgcClass) , shape = 21,fill ="#8da0cb",color = "#8da0cb" ,
               pwidth=0.5,alpha=1,stat="identity",orientation="y", 
               axis.params=list(axis="none", text.angle=-45, hjust=0 ),
               grid.params=list() )  
  p30

  NRPS.1 <- melt(NRPS)
  NRPSmap <- subset(NRPS.1, (id %in% showFamilyID))
  NRPSmap$variable <- factor(NRPSmap$variable ,levels=group_list)
  NRPSmapF = filter(NRPSmap, value > 0)
  
  p31  <- p30 + new_scale_fill() + new_scale_color()
  p32 <- p31 + new_scale_fill() +ggnewscale::new_scale("alpha") +
    geom_fruit(data=NRPSmapF,geom=geom_tile,offset=0.02,mapping=aes(y=id, x=variable, alpha=value, fill=variable),
               pwidth=0.1,
               axis.params=list(axis="none", text.angle=-45, hjust=0 )) +
    scale_fill_manual(values=color_var,guide="none") +
    scale_alpha_continuous(range=c(0.5, 1),name="NRPS abundance", guide=guide_legend(keywidth=0.5, keyheight=0.5, order=7))
  p32


# Layer 9: Product Classes -------------------------RiPPs--------------------------------------
  RiPPsPD.1 <- subset(MergeBGC, select = c( "Level","RiPPs")) 
  RiPPsPD <- subset(RiPPsPD.1, (Level %in% showFamilyID))
  RiPPsPD.S = filter(RiPPsPD, RiPPs > 0)
  RiPPsPD.S$variableNum <- rep(0, times=length(RiPPsPD.S$Level))
  colnames(RiPPsPD.S) <- c("Level"  ,     "bgcClass"  ,    "variableNum")
  
  p33  <- p32 + new_scale_fill() + new_scale_color()
  p34 <- p33 + new_scale_fill() +
    geom_fruit(data=RiPPsPD.S,geom=geom_point,offset=0.07,mapping=aes(x= variableNum,y=Level , size = bgcClass) , shape = 21,fill ="#e78ac3",color = "#e78ac3" ,
               pwidth=0.5,alpha=1,stat="identity",orientation="y", 
               axis.params=list(axis="none", text.angle=-45, hjust=0 ),
               grid.params=list() )  
  
  p34
  
  RiPPs.1 <- melt(RiPPs)
  RiPPsmap <- subset(RiPPs.1, (id %in% showFamilyID))
  RiPPsmap$variable <- factor(RiPPsmap$variable ,levels=group_list)
  RiPPsmapF = filter(RiPPsmap, value > 0)
  
  p35  <- p34 + new_scale_fill() + new_scale_color()
  p36 <- p35 + new_scale_fill() +ggnewscale::new_scale("alpha") +
    geom_fruit(data=RiPPsmapF,geom=geom_tile,offset=0.02,mapping=aes(y=id, x=variable, alpha=value, fill=variable),
               pwidth=0.1,
               axis.params=list(axis="none", text.angle=-45, hjust=0 )) +
    scale_fill_manual(values=color_var,guide="none") +
    scale_alpha_continuous(range=c(0.5, 1),name="RiPPs abundance", guide=guide_legend(keywidth=0.5, keyheight=0.5, order=7))
  p36

  # Layer 11: BGC counts ------------------------------
  totalBGCAbundanceSum <- as.data.frame(rowSums(totalBGCAbundance))
  colnames(totalBGCAbundanceSum) <- c("Total")
  totalBGCAbundanceSum$Fid <- rownames(totalBGCAbundanceSum)
  totalBGCAbundanceSum.1 <- subset(totalBGCAbundanceSum,Fid %in% showFamilyID)
  totalBGCAbundanceSum.S <- subset(totalBGCAbundanceSum.1, select = c( "Total"))
  
  p37  <- p36 + new_scale_fill() + new_scale_color() +ggnewscale::new_scale("fill")
  p38 <-  gheatmap(p37, totalBGCAbundanceSum.S,offset=5.1,width=0.2,colnames_angle=0, colnames_offset_y = 0 ,font.size = 0,hjust=1)+
    #scale_fill_viridis_c(option="A", name="log(BGC counts)")
    scale_fill_gradient(low = "#c6dbef", high = "#08306b", na.value = NA,name="Total BGC counts",
                        guide=guide_colourbar(barwidth = 4, barheight = 0.5,
                                              title.theme= element_text(size = 10),
                                              label.theme = element_text(size = 6),
                                              direction = "horizontal",order=9, title.position = "top") )
  
  p38
  
  # Layer 12: Genome counts ----------------------------------------
  # "Group"      "Level"      "genommeNum"
  totalgenomeNumS <- subset(totalgenomeNum,Level %in% showFamilyID)
  totalgenomeNumS$Group <- factor(totalgenomeNumS$Group ,levels=group_list)
  totalgenomeNumS$GenomeCounts <- log(totalgenomeNumS$genommeNum+1)

  p39  <- p38 + new_scale_fill() + new_scale_color()
  p40 <- p39 + new_scale_fill() +
    geom_fruit(data=totalgenomeNumS,geom=geom_bar,offset=0.4,mapping=aes(y=Level, x=GenomeCounts,fill = Group), 
               pwidth=0.4,alpha=1,stat="identity",orientation="y", 
               axis.params=list(axis="x", text.angle=0, hjust=0 , text.size=2),
               grid.params=list() 
    ) + 
    
    scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list,
                      guide="none") +
    theme(
      legend.background=element_rect(fill=NA), 
      legend.title=element_text(size=7), 
      legend.text=element_text(size=6), 
      legend.spacing.y = unit(0.02, "cm") ) 
  
  p40 

  # Layer 13: label ----------------------------------------
 
  p41  <- p40 + new_scale_fill() + new_scale_color()
  p42 <- p41 + new_scale_fill() +
    #geom_hilight(data=bgPhylum, mapping=aes(node=node,fill = phylum),extendto=12, alpha=0.05, color=NA,size=0.05,show.legend=FALSE)+
    geom_cladelab(data=Phylumlabel, mapping=aes(node=node, label=phylum,color = phylum , offset.text=pos),
                  offset.text=0.5,hjust='center',angle="auto",offset=9,barsize=1,horizontal=FALSE, fontsize=2.5,show.legend=FALSE)+
    scale_color_manual(values=phylumcolorVar,guide="none" ) +
    scale_fill_manual(values=phylumcolorVar,guide="none" ) 

  p42  
    
  
pdf(".//more10.collapse.plot.pdf",width=15,height=10,onefile = FALSE)
p42
dev.off()






