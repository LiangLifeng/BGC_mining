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
library("dplyr")
library("ggpubr")
library("reshape2")
library("ggbreak")
library("patchwork")
library("magick")


# input  -----------------------------------------
group_list =c("Marine","Freshwater","Terrestrial")
color_var = c("#0f7ba2","#43b284","#fab255")

allinfo <- read.table("../data/GCC/GEMs-OMD-FWL-PS.all.antisMash.3Group.region.GCC-GCF.add.info.txt",header=T,  sep = "\t", check.name = F,quote="")
#GCF
gcf <- read.table("../data/GCC/GEMs-OMD-FWL-PS.all.antisMash.3Group.region.GCC-GCF.groupby.mean.GCF.groupby-mean.xls",header=T,check.names=F,row.names=1)
gcfBGCMean <- read.table("../data/GCC/GEMs-OMD-FWL-PS.all.antisMash.3Group.region.GCC-GCF.groupby.mean.GCF.groupby-mean.xls.mergeBGC.txt",header=T,check.names=F,row.names=1)
gcf_groupPre <- read.table("../data/GCC/GEMs-OMD-FWL-PS.all.antisMash.3Group.region.GCC-GCF.stat_GCF_Group.txt",header=T,  sep = "\t", check.name = F,quote="")
gcf_phylumPre <- read.table("../data/GCC/GEMs-OMD-FWL-PS.all.antisMash.3Group.region.GCC-GCF.stat_GCF_Phylum.M.txt",header=T,  sep = "\t", check.name = F,quote="")
gcf_OTUPre <- read.table("../data/GCC/GEMs-OMD-FWL-PS.all.antisMash.3Group.region.GCC-GCF.stat_GCF_OTU.txt",header=T,  sep = "\t", check.name = F,quote="")

# cluster----------------------------------------------------------------------
  hc <- hclust(dist(gcf))

# Tree ----------------------------------------------------------------------
  gcf1 <- ggtree(hc, layout="fan", open.angle = 5, size = 0.15, alpha = 1) + #branch.length = 'none'
    geom_tiplab(align=TRUE, linetype='dotted', linesize=.2, mapping = aes(label = NA)) +
    #scale_color_viridis_c() +
  
    #coord_polar(theta = 'y', start = 0, direction = -1) +
    theme(panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.margin = margin(0,0,0,0))
  
  #gcf1 <- ggtree(hc, layout="fan", open.angle = 5, size = 0.15, alpha = 1) 
  gcf1
  gcf1 <- rotate_tree(gcf1, 90)
  gcf1
  


# Layer 1: lable ----------------------------------------------------------------------
  gcfLable <- as.data.frame(matrix(nrow=length(rownames(gcf)), ncol=0))
  gcfLable$gcfID <- rownames(gcf)
  gcfLable$gcflabel <- gsub("gcf_", "", rownames(gcf))
  rownames(gcfLable) <- rownames(gcf)

  gcf1
  gcf2 <- gcf1 %<+% gcfLable + geom_tiplab(aes(label=gcflabel),size=1,offset=0,fontface="bold")
  gcf2


    
# Layer 2: Group  Prevalence ----------------------------------------------------------------------
  gcf_groupPre$Prevalence <- gcf_groupPre$Prevalence*100
  gcf_groupPre$Group <- factor(gcf_groupPre$Group,levels=group_list)
  
  gcf3  <- gcf2 + new_scale_fill() + new_scale_color()

  gcf4 <- gcf3 + new_scale_fill() +
    geom_fruit(data=gcf_groupPre,geom=geom_col,mapping=aes(y=GID, x=Prevalence, fill=Group), 
      pwidth=0.1,alpha=0.6,stat="identity",orientation="y", offset = 0.09,
      axis.params=list(axis="x", text.angle=-45, hjust=0 ),
      grid.params=list() ) + 
    scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list,
                      guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
    theme(
      legend.background=element_rect(fill=NA), 
      legend.title=element_text(size=7), 
      legend.text=element_text(size=6), 
      legend.spacing.y = unit(0.02, "cm") ) 
  
  gcf4

# Layer 3: Product Classes -------------------------gcfBGCMean--------------------------------------
  
  bgc_class = c("RiPPs","NRPS","PKS","Terpene","Others")
  bgc_colors =c ("#e78ac3","#8da0cb","#fc8d62","#66c2a5","#cccccc")    
  gcfBGCMean$GID <- rownames(gcfBGCMean)
  gcfBGCMeanL <- melt(gcfBGCMean,id.vars = c("GID"))
  gcfBGCMeanL$variable <- factor(gcfBGCMeanL$variable,levels=bgc_class)
  
  gcfBGCMeanL$variableNum <- as.numeric(gcfBGCMeanL$variable)
  

  
  gcfBGCMean.S = filter(gcfBGCMeanL, value > 0)
  
  
  gcf5  <- gcf4 + new_scale_fill() + new_scale_color()
  
  gcf6 <- gcf5 + new_scale_fill() +
    geom_fruit(data=gcfBGCMean.S,geom=geom_point,mapping=aes(x= variableNum,y=GID,  fill = variable, color = variable , size = value) , shape = 21,
               pwidth=0.2,alpha=1,stat="identity",orientation="y", 
               axis.params=list(axis="x", text.angle=-45, hjust=0 ),
               grid.params=list() ) + 
    
    scale_fill_manual(values=bgc_colors,labels=bgc_class,limits = bgc_class, breaks = bgc_class,
                      guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
    scale_color_manual(values=bgc_colors) +
    scale_size_continuous(range = c(0, 3))+
    theme(
      legend.background=element_rect(fill=NA), 
      legend.title=element_text(size=7), 
      legend.text=element_text(size=6), 
      legend.spacing.y = unit(0.02, "cm") ) 
  
  
  gcf6
  

  
  
  # Layer 4: phylum  Prevalence -----------------------------------------gcf_phylumPreN-----------------------------
  gcf_phylumPre$Group <- factor(gcf_phylumPre$Group,levels=group_list)
  gcf_phylumPreN <- gcf_phylumPre %>% mutate(PT=if_else(Prevalence >10 , "B" ,"S"))
  gcf_phylumPreN <- subset(gcf_phylumPreN, (PT %in% c("B")))
  gcf_phylumPreLess10 <- gcf_phylumPreN %>%group_by(GID)%>%summarize(Prevalence = 100- sum(Prevalence))%>%mutate(level = "Other")
  gcf_phylumPreB <- subset(gcf_phylumPreN, select = c("GID","Prevalence","level"))
  gcf_phylumMergeLess <-rbind(gcf_phylumPreB,gcf_phylumPreLess10)
  
  
  phyla_order = c(
    "Marine:p__Proteobacteria","Terrestrial:p__Proteobacteria","Freshwater:p__Proteobacteria","Marine:p__Bacteroidota",
    "Freshwater:p__Bacteroidota","Terrestrial:p__Bacteroidota","Marine:p__Marinisomatota","Marine:p__Actinobacteriota","Terrestrial:p__Actinobacteriota",
    "Freshwater:p__Actinobacteriota","Marine:p__Chloroflexota","Marine:p__Cyanobacteria","Marine:p__Verrucomicrobiota","Freshwater:p__Verrucomicrobiota",
    "Marine:p__Planctomycetota","Freshwater:p__Planctomycetota","Terrestrial:p__Firmicutes","Terrestrial:p__Firmicutes_A",
    "Freshwater:p__Patescibacteria","Marine:p__SAR324","Freshwater:p__Chloroflexota","Marine:Other","Freshwater:Other","Terrestrial:Other","Other" 
    )

  phyla_colors = c(
    "Marine:p__Proteobacteria" = "#fdcdac",
    "Marine:p__Bacteroidota" = "#a6d854",
    "Terrestrial:p__Proteobacteria" = "#fc8d62",
    "Marine:p__Marinisomatota" = "#b2df8a",
    "Marine:p__Actinobacteriota" = "#8da0cb",
    "Marine:p__Chloroflexota" = "#e6ab02",
    "Freshwater:p__Proteobacteria" = "#b64f32",
    "Freshwater:p__Bacteroidota" = "#66a61e",
    "Terrestrial:p__Bacteroidota" = "#1b9e77",
    "Marine:p__Cyanobacteria" = "#e41a1c",  #811e18
    "Marine:p__Verrucomicrobiota" = "#a6cee3",
    "Marine:p__Planctomycetota" = "#ffd92f",
    "Freshwater:p__Verrucomicrobiota" = "#1f78b4",
    "Terrestrial:p__Firmicutes" = "#e78ac3",
    "Terrestrial:p__Actinobacteriota" = "#7570b3",
    "Freshwater:p__Actinobacteriota" = "#984ea3",
    "Freshwater:p__Planctomycetota" = "#e5c494",
    "Freshwater:p__Patescibacteria" = "#4292c6",
    "Marine:p__SAR324" = "#ffff33",
    "Freshwater:p__Chloroflexota" = "#a6761d",
    "Terrestrial:p__Firmicutes_A" = "#f4cae4",
    "Marine:Other" = "#bababa",
    "Freshwater:Other" = "#878787",
    "Terrestrial:Other" = "#4d4d4d",
    "Other" = "#e0e0e0"
  )
  
  gcf_phylumMergeLess$level <- factor(gcf_phylumMergeLess$level,levels=phyla_order)
  
  gcf7  <- gcf6 + new_scale_fill() + new_scale_color()
  
  gcf8 <- gcf7 + new_scale_fill() +
    geom_fruit(data=gcf_phylumMergeLess,geom=geom_col,mapping=aes(y=GID, x=Prevalence,fill = level), 
               pwidth=0.4,alpha=1,stat="identity",orientation="y", 
               axis.params=list(axis="x", text.angle=-45, hjust=0 ),
               grid.params=list() ) + 

    scale_fill_manual(values=phyla_colors,labels=phyla_order,limits = phyla_order, breaks = phyla_order,
                     guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
    theme(
      legend.background=element_rect(fill=NA), 
      legend.title=element_text(size=7), 
      legend.text=element_text(size=6), 
      legend.spacing.y = unit(0.02, "cm") ) 
  
  
  gcf8
  
  # Layer 5: BGC -----------------------------------------gcf-----------------------------
  
  sbgc_order = c("glycocin", "sactipeptide", "cyanobactin", "proteusin", "microviridin", "LAP", "lanthipeptide", "linaridin", 
                 "bacteriocin", "thiopeptide", "bottromycin", "head_to_tail", "RaS-RiPP", "lassopeptide", 
                 "lipolanthine", "NRPS", "NRPS-like", "PKS-like", "T1PKS", "T2PKS", "T3PKS", "transAT-PKS", "transAT-PKS-like", 
                 "terpene", "acyl_amino_acids", "amglyccycl", "arylpolyene", "betalactone", "blactam", "butyrolactone", "CDPS", 
                 "ectoine", "furan", "hglE-KS", "hserlactone", "indole", "ladderane", "melanin", "NAGGN", "nucleoside", "oligosaccharide",
                 "other", "PBDE", "phenazine", "phosphonate", "PpyS-KS", "PUFA", "resorcinol", "siderophore", "TfuA-related", "tropodithietic-acid", 
                 "phosphoglycolipid", "thioamide-NRP", "fused"
    
  )
    sbgc_colors = c(
    "glycocin" = "#e78ac3",
    "sactipeptide" = "#e78ac3",
    "cyanobactin" = "red",
    "proteusin" = "blue",
    "microviridin" = "#e78ac3",
    "LAP" = "red",
    "lanthipeptide" = "green",
    "linaridin" = "blue",
    "bacteriocin" = "#e78ac3",
    "thiopeptide" = "red",
    "bottromycin" = "#e78ac3",
    "head_to_tail" = "green",
    "RaS-RiPP" = "#e78ac3",
    "lassopeptide" = "red",
    "lipolanthine" = "#e78ac3",
    "NRPS" = "#8da0cb",
    "NRPS-like" = "#8da0cb",
    "PKS-like" = "#fc8d62",
    "T1PKS" = "#fc8d62",
    "T2PKS" = "#fc8d62",
    "T3PKS" = "#fc8d62",
    "transAT-PKS" = "#fc8d62",
    "transAT-PKS-like" = "#fc8d62",
    "terpene" = "#66c2a5",
    "acyl_amino_acids" = "#b3b3b3",
    "amglyccycl" = "#b3b3b3",
    "arylpolyene" = "#b3b3b3",
    "betalactone" = "#b3b3b3",
    "blactam" = "#b3b3b3",
    "butyrolactone" = "#b3b3b3",
    "CDPS" = "#b3b3b3",
    "ectoine" = "#b3b3b3",
    "furan" = "#b3b3b3",
    "hglE-KS" = "#b3b3b3",
    "hserlactone" = "#b3b3b3",
    "indole" = "#b3b3b3",
    "ladderane" = "#b3b3b3",
    "melanin" = "#b3b3b3",
    "NAGGN" = "#b3b3b3",
    "nucleoside" = "#b3b3b3",
    "oligosaccharide" = "#b3b3b3",
    "other" = "#b3b3b3",
    "PBDE" = "#b3b3b3",
    "phenazine" = "#b3b3b3",
    "phosphonate" = "#b3b3b3",
    "PpyS-KS" = "#b3b3b3",
    "PUFA" = "#b3b3b3",
    "resorcinol" = "#b3b3b3",
    "siderophore" = "#b3b3b3",
    "TfuA-related" = "#b3b3b3",
    "tropodithietic-acid" = "#b3b3b3",
    "phosphoglycolipid" = "#b3b3b3",
    "thioamide-NRP" = "#b3b3b3",
    "fused" = "#b3b3b3"
  )
  
  
  gcf_filter <- gcf
  gcf_filter$GID <- rownames(gcf_filter)
  gcf_filterL <- melt(gcf_filter,id.vars = c("GID"))
  
  gcf_filterL$variable <- factor(gcf_filterL$variable,levels=sbgc_order)
  gcf_filterL$variableNum <- as.numeric(gcf_filterL$variable)
  gcf_filterL.S = filter(gcf_filterL, value > 0)
  
  gcf9  <- gcf8 + new_scale_fill() + new_scale_color()
  
  gcf10 <- gcf9 + new_scale_fill() +
    geom_fruit(data=gcf_filterL.S,geom=geom_point,mapping=aes(x= variableNum,y=GID,  fill = variable,color = variable ,  size = value) , shape = 21, #
               pwidth=1.2,alpha=1,stat="identity",orientation="y", 
               axis.params=list(axis="x", text.angle=-45, hjust=0 ,nbreak = 38),
               grid.params=list() ) + 
    
    scale_fill_manual(values=sbgc_colors,labels=sbgc_order,limits = sbgc_order, breaks = sbgc_order,
                      guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
    scale_color_manual(values=sbgc_colors,labels=sbgc_order,limits = sbgc_order, breaks = sbgc_order,
                       guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
    #scale_color_manual(values=sbgc_order) +
    
    scale_size_continuous(range = c(0, 1))+
    theme(
      legend.background=element_rect(fill=NA), 
      legend.title=element_text(size=7), 
      legend.text=element_text(size=6), 
      legend.spacing.y = unit(0.02, "cm") ) 
  
  
  gcf10
  
  
  
  
  
  # Layer 6: BGC counts -----------------------------------------gcf_OTUPre-----------------------------
  gcf_OTUPre$Group <- factor(gcf_OTUPre$Group,levels=group_list)
  gcf_OTUPreCounts <- subset(gcf_OTUPre, select = c( "GID","Total"))
  gcf_OTUPreCountsU <- unique(gcf_OTUPreCounts)
  rownames(gcf_OTUPreCountsU) <- gcf_OTUPreCountsU$GID
  gcf_OTUPreCountsU$Total <- log(gcf_OTUPreCountsU$Total)
  gcf_OTUPreCountsUPH <- subset(gcf_OTUPreCountsU, select = c("Total"))
  
  #gcf_OTUPreCountsUPH <- gcf_OTUPreCountsU %>%
  #  mutate(counts=
  #           if_else(Total < 13,"< 13",
  #                   if_else(Total < 133,"< 133",
  #                           if_else(Total < 1300,"< 1300","> 1300")))
  #         )
  
  #gcf_OTUPreCountsUPH <- subset(gcf_OTUPreCountsUPH, select = c("counts"))
  
  gcf11  <- gcf10 + new_scale_fill() + new_scale_color()

  gcf12 <- gheatmap(gcf11, gcf_OTUPreCountsUPH,offset=3.05,width=0.04,colnames_angle=0, colnames_offset_y = 0 ,font.size = 2,hjust=1)+
    #scale_fill_viridis_c(option="A", name="log(BGC counts)")
    scale_fill_gradient(low = "grey90", high = "black", na.value = NA,name="log(BGC counts)",
                        guide=guide_colourbar(barwidth = 4, barheight = 0.5,
                                              title.theme= element_text(size = 10),
                                              label.theme = element_text(size = 6),
                                              direction = "horizontal",order=9) )
  
  gcf12
  
  # Layer 7: Genome counts -----------------------------------------allinfo-----------------------------
  
  gcf_genome <- subset(allinfo, select = c( "GCF.fit_predict.clusterID","Genome_id","Freeliving"))
  gcf_genomeS <- unique(gcf_genome)
  gcf_genomeCounts <- gcf_genomeS %>%
    group_by(GCF.fit_predict.clusterID,Freeliving) %>%
    summarise(GenomeCounts = n())
  
  gcf_genomeCounts$Group <- factor(gcf_genomeCounts$Freeliving,levels=group_list)
  gcf_genomeCounts$GenomeCounts <- log(gcf_genomeCounts$GenomeCounts+1)
  

  gcf13  <- gcf12 + new_scale_fill() + new_scale_color()
  
  gcf14 <- gcf13 + new_scale_fill() +
    geom_fruit(data=gcf_genomeCounts,geom=geom_bar,offset=0.09,mapping=aes(y=GCF.fit_predict.clusterID, x=GenomeCounts,fill = Group), 
               pwidth=0.2,alpha=1,stat="identity",orientation="y", 
               #axis.params=list(axis="x", text.angle=-45, hjust=0 ),
               #grid.params=list() 
               ) + 

    scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list,
                      guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
    theme(
      legend.background=element_rect(fill=NA), 
      legend.title=element_text(size=7), 
      legend.text=element_text(size=6), 
      legend.spacing.y = unit(0.02, "cm") ) 
  
  
  gcf14
  
  
  
  
# output  -----------------------------------------
  
#ggsave("../data/GCC/GCF-overview.pdf", width = 20, height = 10, unit = "cm")

  pdf("../data/GCC/GCF-overview.pdf",width=15,height=10,onefile = FALSE)
  gcf14
  dev.off()
  png("../data/GCC/GCF-overview.png",width = 15,height = 10,units='in',res=900)
  gcf14
  dev.off()

