library("ggtree")
library("ape")
library("MicrobiotaProcess")
library("phylobase")
library("tidytree")
library("RColorBrewer")
library('cowplot')
library("ggtreeExtra")
library("treeio")
library("ggnewscale")
library("tidyverse")

colorVar <- c( "#66c2a5","#8da0cb","#fc8d62","grey",
               "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#8dd3c7","#80b1d3","#bebada","#fb8072","#fb8072","#fdb462","#b3de69","#fccde5","#bc80bd","#ccebc5","#ffed6f","#bebada","#bebada","#bebada","#fb8072","#fb8072","#fccde5","#bebada","#fb8072","#bebada","#fb8072","#fb8072","#fb8072","#fb8072")

orderList <- c("Freshwater","Marine","Terrestrial","Other",
               "Proteobacteria","Bacteroidota","Actinobacteriota","Verrucomicrobiota","Chloroflexota","Planctomycetota","Patescibacteria","Acidobacteriota","Marinisomatota","Cyanobacteria","Desulfobacterota","Firmicutes","Firmicutes_A","Myxococcota","Omnitrophota","Bdellovibrionota","Gemmatimonadota","Nitrospirota","Spirochaetota","Desulfobacterota_F","Desulfobacterota_D","Desulfobacterota_B","Firmicutes_B","Firmicutes_D","Bdellovibrionota_C","Desulfobacterota_E","Firmicutes_C","Desulfobacterota_G","Firmicutes_E","Firmicutes_H","Firmicutes_G","Firmicutes_F")

RiPPs <- c("bacteriocin","sactipeptide","microviridin","lanthipeptide","lassopeptide","thiopeptide","LAP","RaS_RiPP","linaridin")
PKS <- c("T1PKS","T3PKS","transAT_PKS","transAT_PKS_like")
Others <- c("siderophore","betalactone","NAGGN","PUFA","ectoine","hserlactone","indole","resorcinol","TfuA_related","furan","butyrolactone","oligosaccharide","arylpolyene","acyl_amino_acids","ladderane","hglE_KS")
NRPS <- c("NRPS","NRPS_like")

tree <-read.tree("../phylophlan_result/fna.tre.treefile")
nichGroup <- read.table("../anno//cluster.nich.txt",header=T,  sep = "\t", check.name = F,quote="")
phylumAnontxt1 <- read.table("../anno/cluster.phylum.anno.txt",header=T, sep = "\t", check.name = F,quote="")
phylumAnontxt <- read.table("../anno/cluster.phylum.anno.txt",header=T, sep = "\t", check.name = F,quote="",row.names = 1)
genomeStat <- read.table("../anno/cluster.genoeme.number.txt",header=T, sep = "\t", check.name = F,quote="")
BGCStat <- read.table("../anno/cluster.BGC.txt.2",header=T, sep = "\t", check.name = F,quote="",row.names = 1)
mergeBGC <- read.table("../anno/cluster.Meger.BGC.txt",header=T, sep = "\t", check.name = F,quote="",row.names = 1)
bgPhylum <- read.table("../anno/phylum.marker.txt.2",header=T, sep = "\t", check.name = F,quote="")
Phylumlabel <- read.table("../anno/phylum.marker.txt.3",header=T, sep = "\t", check.name = F,quote="")
mergeBGC <- read.table("../anno/cluster.Meger.BGC.txt.3",header=T, sep = "\t", check.name = F,quote="",row.names = 1)


tt <- merge(phylumAnontxt1,nichGroup,by="id",all=TRUE)
rownames(tt) <- tt$id
phylumInfo <-subset(tt,select=-c(id))


RiPPsM <-subset(mergeBGC,select=c(RiPPs)) 
NRPSM <-subset(mergeBGC,select=c(NRPS)) 
PKSM <-subset(mergeBGC,select=c(PKS)) 
TerpeneM <-subset(mergeBGC,select=c(Terpene)) 
OthersM <-subset(mergeBGC,select=c(Others)) 


p0 <- ggtree(tree,size = 0.1 ,layout="circular")   + 
  geom_hilight(data=bgPhylum, mapping=aes(node=node,fill = phylum),extendto=1.35, alpha=0.2, color=NA,size=0.05,show.legend=FALSE)+
  geom_cladelab(data=Phylumlabel, mapping=aes(node=node, label=phylum,offset.text=pos),hjust=0.5,angle="auto",barsize=NA,horizontal=FALSE, fontsize=2,fontface="italic",show.legend=FALSE)+
  theme(legend.position = c(.9, .6),
        legend.background = element_rect(),
        legend.key = element_blank(), # removes the border
        legend.key.size = unit(0.2, 'cm'), # sets overall area/size of the legend
        legend.text = element_text(size = 8), # text size
        title = element_text(size = 10))

p <- open_tree(p0,5)
#添加门水平分类标记

p1 <- p %<+% phylumAnontxt1  + geom_tippoint(aes(color=Phylum),size = 0.15,show.legend=FALSE)

p2 <- gheatmap(p1, phylumInfo,offset=-0.05, width=0.1,colnames_angle=45, colnames_offset_y = 0 ,font.size = 2,color =NA,hjust=1)+
  scale_color_manual(values=colorVar,limits = orderList ) +
  scale_fill_manual(values=colorVar,limits = orderList ) 

#RiPPs #00441b "#1b9e77","#66c2a5","#b3e2cd","#ffffff"
p3  <- p2 + new_scale_fill() + new_scale_color()
p5 <- gheatmap(p3, RiPPsM,offset=0.05, width=0.06,colnames_angle=45, colnames_offset_y = 0 ,font.size = 2,color =NA,hjust=1)+
  scale_color_manual(values=c("#00441b","#238b45","#41ab5d","#74c476","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0"),name="RiPPs",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=2) ) +
  scale_fill_manual(values=c("#00441b","#238b45","#41ab5d","#74c476","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0"),name="RiPPs",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=2) ) 

#NRPS #543005 #d95f03	#fc8d62	#fdcdac
p6  <- p5 + new_scale_fill() + new_scale_color()
p7 <- gheatmap(p6, NRPSM,offset=0.12, width=0.06,colnames_angle=45, colnames_offset_y = 0 ,font.size = 2,color =NA,hjust=1)+
  scale_color_manual(values=c("#7f3b08","#b35806","#e08214","#fdb863","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0"),name="NRPS",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=3) ) +
  scale_fill_manual(values=c("#7f3b08","#b35806","#e08214","#fdb863","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0"),name="NRPS",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=3) ) 



#PKS #3f007d #7570b4	#8da0cb	#cbd5e8
p8  <- p7 + new_scale_fill() + new_scale_color()
p9 <- gheatmap(p8, PKSM,offset=0.19, width=0.06,colnames_angle=45, colnames_offset_y = 0 ,font.size = 2,color =NA,hjust=1)+
  scale_color_manual(values=c("#2d004b","#542788","#8073ac","#b2abd2","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0"),name="PKS",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=4) ) +
  scale_fill_manual(values=c("#2d004b","#542788","#8073ac","#b2abd2","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0"),name="PKS",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=4) ) 


#Terpene #8e0152 #e7299a	#e78ac3	#f4cae4
p10  <- p9 + new_scale_fill() + new_scale_color()
p11 <- gheatmap(p10, TerpeneM,offset=0.26, width=0.06,colnames_angle=45, colnames_offset_y = 0 ,font.size = 2,color =NA,hjust=1)+
  scale_color_manual(values=c("#8e0152","#c51b7d","#de77ae","#f1b6da","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0") ,name="Terpene",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=5) ) +
  scale_fill_manual(values=c("#8e0152","#c51b7d","#de77ae","#f1b6da","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0"),name="Terpene",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=5)  ) 


#Others #000000 #252525	#636363	#cccccc
p12  <- p11 + new_scale_fill() + new_scale_color()
p13 <- gheatmap(p12, OthersM,offset=0.33, width=0.06,colnames_angle=45, colnames_offset_y = 0 ,font.size = 2,color =NA,hjust=1)+
  scale_color_manual(values=c("#1a1a1a","#4d4d4d","#878787","#bababa","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0") ,name="Others",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=6) ) +
  scale_fill_manual(values=c("#1a1a1a","#4d4d4d","#878787","#bababa","#ffffff"	),limits = c(">10","2-10","1-2","0-1","0") ,name="Others",guide=guide_legend(legend.key.size = unit(0.2, 'cm'),order=6) ) 


#genome number
p14  <- p13 + new_scale_fill() + new_scale_color()
p15 <- p14 + geom_fruit(data=genomeStat, geom=geom_bar,
                        mapping=aes(y=ID, x=GenomeNumber),
                        offset=0.46,
                        orientation="y", 
                        color="grey30",
                        stat="identity")



pdf("tree.Merga.BGC.pdf",width=12,height=10,onefile = FALSE)
p15
dev.off()




