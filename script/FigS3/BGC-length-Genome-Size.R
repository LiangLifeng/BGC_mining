library("ggpubr")
library("ggplot2")
library("reshape")
library("ggprism")
library("getopt")
library("ggpubr")
library("rstatix")
library("gtools")
library("tidyverse")
library("multcomp")
library("nparcomp")
library("multcompView")
#library("do")
library("dplyr")
library("PMCMRplus")
library("PMCMR")

spec <- matrix(c(
	'help','h',0,'logical',"help info author: 多维组学 Momics\tv1:2021-04-27 \nV2:2021-12-02 add leter sign\nv3:2022-03-11 plot each Variable ",
	'profile','i','1','character','plot data, long data type',
	'cGroup','g','1','character','the mapping column ID of you want to draw group by,default ="Group"',
	'color','c','2','character','color list ,link by - , use the HEX coding which remove the #. default = 1b9e77:d95f02:7570b3:e7298a:66a61e:e6ab02:a6761d:666666:999999',
	'group','k','2','character','drow group name list ,link by : , eg, group1:group2 ,default="All", if you only want to draw a part, need to give -s',
	'testMethod','t','2','character',' statistical tests method, we only support wilcox_test, but "rstatix" support more kind of method, you can see in https://rpkgs.datanovia.com/rstatix/ and replace the "wilcox_test" which method you need,default="wilcox_test" ',
	'padjust','q','2','character',' adjust pvalue , [Y,N],default="Y" ',
	'plabel','l','2','character',' pvalue label form, it can be show as "value" or "signif" or "value:signif" or "signif:value",default="signif" ',
	'eachplotYN','e','2','character','plot each variable or not [Y/N],default="N"',
	'outdir','o','2','character','output dir,default=./',
	'prefix','p','2','character','prefix,default=All'

),ncol=5,byrow=TRUE)


opt <- getopt(spec)
# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help)) {
        cat(getopt(spec, usage=TRUE));
        q(status=1);
}

if( is.null(opt$outdir) ) { opt$outdir <- getwd()}
if( is.null(opt$prefix) ) { opt$prefix <- 'All'}
if( is.null(opt$color) ) { opt$color <- '1b9e77:d95f02:7570b3:e7298a:66a61e:e6ab02:a6761d:666666:999999'}
if( is.null(opt$part) ) { opt$part <- 'N'}
if( is.null(opt$group) ) { opt$group <- 'All'}
if( is.null(opt$specified) ) { opt$specified <- 'All'}
if( is.null(opt$padjust) ) { opt$padjust <- 'Y'}
if( is.null(opt$eachplotYN) ) { opt$eachplotYN <- 'N'}
if( is.null(opt$testMethod) ) { opt$testMethod <- 'wilcox_test'}
if( is.null(opt$plabel) ) { opt$plabel <- 'signif'}
if( is.null(opt$cGroup) ) { opt$cGroup <- 'Group'}


group_list = unlist(strsplit(opt$group, ":"))
plabel  <-'p.adj.signif'
color_var = unlist(strsplit(opt$color, ":"))
color_var = c(paste("#",color_var,sep=""))

#Function
nparcomp_compare1 <- function(data,group,compare,value){
  a <- data.frame(stringsAsFactors = F)#做一个空的数据框
  pv <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])#统计需要运行多重比较的次数
  for (i in type)#进行type次多重比较
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]#根据指定的i去取相应的数据集出来
    names(sub_dat)[names(sub_dat)==compare] <- 'g1' ## 重命名方便后面使用
    names(sub_dat)[names(sub_dat)==value] <- 'value' ## 重命名方便后面使用
    sub_dat$g1 <- factor(sub_dat$g1)#将列转化成因子以进行多重比较
    k <- nparcomp::nparcomp(value ~ g1,data=sub_dat, alternative = "two.sided")
    b <- k$Analysis
    dif <- b$p.Value
    names(dif) <- b$Comparison
    difname <- names(dif) %>% 
      Replace(from = ' , ',to='-') %>% 
      Replace(from=c('p\\(','\\)'),to='')#正则表达式替换
    temp_name <- data.frame(tn=difname) %>% 
      separate(col = 'tn',sep = '-',into = c('n1','n2'))
    difname <- paste(temp_name$n2,temp_name$n1,sep = '-') %>%
      Replace(from = ' - ',to='-')#正则表达式替换
    names(dif) <- difname
    difL <- multcompLetters(dif)
    dif[is.na(dif)]<-1
    labels <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    labels$compare = rownames(labels)
    labels$factor <- i
    mean_sd <- as.matrix(aggregate(sub_dat$value, by=list(sub_dat$g1), FUN=function(x)c(std=sd(x),mean=mean(x),Max = max(x) )))
    mean_sd <- as.data.frame(mean_sd)
    names(mean_sd) <- c('compare','std','mean','Max')
	    
    a <- rbind(a,merge(mean_sd,labels,by='compare'))
    b$factor <- i
    pv <- rbind(pv,b)
    
  }
  
  names(a) <- c(compare,'std','mean','Max','Letters',group)
  names(pv) <- c("Comparison","Estimator" ,"Lower" ,"Upper" ,"Statistic","p.Value",group)
  return(list(leter = a, pvalue = pv))
}


##2
PMCMR_compare1 <- function(data,group,compare,value){
  a <- data.frame(stringsAsFactors = F)
  pv <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)
    
    k <- kwAllPairsNemenyiTest(value ~ g1,data=sub_dat)
    n <- as.data.frame(k$p.value)
    h <- n %>%
      mutate(compare=rownames(n)) %>%
      gather(group,p,-compare,na.rm = TRUE) %>%
      unite(compare,group,col="G",sep="-")
    dif <- h$p
    names(dif) <- h$G
    b <- h
    difL <- multcompLetters(dif)
    K.labels <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    K.labels$compare = rownames(K.labels)
    K.labels$type <- i
    mean_sd <- as.matrix(aggregate(sub_dat$value, by=list(sub_dat$g1), FUN=function(x)c(std=sd(x),mean=mean(x),Max = max(x) )))
    mean_sd <- as.data.frame(mean_sd)
    names(mean_sd) <- c('compare','std','mean','Max')
    a <- rbind(a,merge(mean_sd,K.labels,by='compare'))
    b$factor <- i
    pv <- rbind(pv,b)

  }
  names(a) <- c(compare,'std','mean','Max','Letters',group)
  names(pv) <- c("Comparison","p.Value",group)
  return(list(leter = a, pvalue = pv))
  
}


##1 profile
profile <- read.table(opt$profile,header=T, sep = "\t", row.names=1,check.name = F)
profile$ID <- rownames(profile)
profile$Group <- profile[,opt$cGroup]
data <- subset(profile,Group %in% group_list)
data$Group <- factor(data$Group,levels=group_list)
data$Genome_size <- as.numeric(data$Genome_size)
data$Genome_sizeM <- data$Genome_size / 1000000
data$regionLength <- as.numeric(data$regionLength)
data$regionLengthK <- data$regionLength / 1000

# remove not exit in 3 Group BGC 
  unpervasiveBGC = c("amglyccycl","bottromycin","furan","head_to_tail","oligosaccharide","PBDE","phenazine","phosphoglycolipid","PpyS-KS","RaS-RiPP","tropodithietic-acid")
  removeBGCdata <- subset(profile, !(BGC %in% unpervasiveBGC))
  removeBGCdata$Group <- factor(removeBGCdata$Group,levels=group_list)
  removeBGCdata$Genome_size <- as.numeric(removeBGCdata$Genome_size)
  removeBGCdata$Genome_sizeM <- removeBGCdata$Genome_size / 1000000
  removeBGCdata$regionLength <- as.numeric(removeBGCdata$regionLength)
  removeBGCdata$regionLengthK <- removeBGCdata$regionLength / 1000

dim(removeBGCdata)
dim(data)

#################### @ 1 ####################
##3 group Genome_size
    box1 <- ggboxplot(data, x = "Group", y = "Genome_sizeM", fill = "Group",size = 0.2, outlier.size = 0.02, alpha=0.7,outlier.shape = 20,
                    ggtheme = theme_bw() +
                        theme(
                            axis.text.x = element_text(color="black",size=0,angle=-45,hjust= 0.1 ,vjust = 0 ),
                            axis.text.y = element_text(color="black",size=10),
                            axis.title.y=element_text(color="black",size=10),
                            axis.line = element_line(color="black"),
                            axis.ticks = element_line(color="black"),
                            legend.position = "none",
                            legend.text = element_text(size = 6, colour = "black"),
                            legend.title = element_text(size = 6),
                            legend.key.width = unit(0.3, 'cm'),
                            legend.key.size = unit(0, 'lines'),
                            panel.grid = element_blank(),
                            panel.background = element_blank()),
                    legend = "none",title = "",xlab = '', ylab = 'Mbp',width = 0.5)+ 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)
    
    stat.test <- data %>%
            wilcox_test(Genome_sizeM  ~ Group) %>%
            adjust_pvalue(method = "BH") %>%
            add_xy_position(x = "Group", dodge = 0.8)

    Sbox1 <- box1 +stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01,y.position = c(12, 12.5, 12), bracket.shorten = 0.05, size = 1.5)

    write.table(as.matrix(stat.test),file=paste(opt$outdir,"/",opt$prefix,"Genome_size.wilcox_test.qvalue.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-GenomeSize.pdf",sep=""),width = 1.2,height = 2)
    box1
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-GenomeSize.png",sep=""),width = 1.2,height = 2,units='in',res=600)
    box1
    dev.off()

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-GenomeSize.sign.pdf",sep=""),width = 1.2,height = 2)
    Sbox1
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-GenomeSize.sign.png",sep=""),width = 1.2,height = 2,units='in',res=600)
    Sbox1
    dev.off()

#################### @ 1 ####################

#################### @ 2 ####################
#BGC class 3 group comparing
    box2 <- ggboxplot(data, x = "BGC_Class", y = "CDSNum", fill = "Group",size = 0.2, outlier.size = 0.01, alpha=0.8, outlier.shape = 20,
                    ggtheme = theme_bw() +
                      theme(axis.text.x = element_text(color="black",size=10,angle=0, hjust= 0.5 ,vjust = 1.5 ),
                            axis.text.y = element_text(color="black",size=10),
                            axis.title.y=element_text(color="black",size=10),
                            axis.line = element_line(color="black"),
                            axis.ticks = element_line(color="black"),
                            strip.text.x = element_text(size = 8, colour = "black"), 
                            strip.background  = element_blank(),
                            legend.position = "none",
                            legend.text = element_text(size = 6, colour = "black"),
                            legend.title = element_text(size = 6),
                            legend.key.width = unit(0.3, 'cm'),
                            legend.key.size = unit(0, 'lines'),
                            panel.grid = element_blank(),
                            panel.background = element_blank()),
                    legend = "none",title = "",xlab = 'BGC Class', ylab = '# of BGC genes',width = 0.5)+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)

    #order
    order_test <- data %>%
      group_by(BGC_Class) %>%
      kruskal_test(CDSNum ~ Group) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj")

    BGC_Class_order = as.character(order_test[order(order_test$p.adj),]$BGC_Class)
    #leter sign
    df3 <- PMCMR_compare1(data,'BGC_Class','Group','CDSNum')
    write.table(df3$pvalue,file=paste(opt$outdir,"/",opt$prefix,".BGC_Class.CDSNum.nonparametric.pvalue.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)
    write.table(df3$leter,file=paste(opt$outdir,"/",opt$prefix,".BGC_Class.CDSNum.nonparametric.leter.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)
    df3$leter$std <- as.numeric(df3$leter$std)
    df3$leter$mean <- as.numeric(df3$leter$mean)
    df3$leter$Max <- as.numeric(df3$leter$Max)
    df3$leter$Group <- factor(df3$leter$Group,levels=group_list)
    df3$leter$BGC_Class <- factor(df3$leter$BGC_Class,levels=BGC_Class_order)
    df3$leter <- arrange(df3$leter, order(df3$leter$BGC_Class),order(df3$leter$Group))
    splot <- df3$leter[order(df3$leter[,'BGC_Class'],df3$leter[,'Group']),]

    Sbox2 <- box2 + geom_text(data=splot,aes(x=BGC_Class,y=max(Max*1.01),label=Letters,colour=Group),angle=90,size = 3,position = position_dodge2(0.9),show.legend = FALSE)+scale_color_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-BGC_Class.CDSNum.pdf",sep=""),width = 3.5,height = 2)
    box2
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-BGC_Class.CDSNum.png",sep=""),width = 3.5,height = 2,units='in',res=600)
    box2
    dev.off()

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-BGC_Class.CDSNum.sign.pdf",sep=""),width = 3.5,height = 2)
    Sbox2
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-BGC_Class.CDSNum.sign.png",sep=""),width = 3.5,height = 2,units='in',res=600)
    Sbox2
    dev.off()

#################### @ 2 ####################



#################### @ 3 S ####################
#BGC class 3 group  regionLengthK comparing
    box3 <- ggboxplot(data, x = "BGC_Class", y = "regionLengthK", fill = "Group",size = 0.2, outlier.size = 0.01, alpha=0.8, outlier.shape = 20,
                    ggtheme = theme_bw() +
                      theme(axis.text.x = element_text(color="black",size=10,angle=0, hjust= 0.5 ,vjust = 1.5 ),
                            axis.text.y = element_text(color="black",size=10),
                            axis.title.y=element_text(color="black",size=10),
                            axis.line = element_line(color="black"),
                            axis.ticks = element_line(color="black"),
                            strip.text.x = element_text(size = 8, colour = "black"), 
                            strip.background  = element_blank(),
                            legend.position = "none",
                            legend.text = element_text(size = 6, colour = "black"),
                            legend.title = element_text(size = 6),
                            legend.key.width = unit(0.3, 'cm'),
                            legend.key.size = unit(0, 'lines'),
                            panel.grid = element_blank(),
                            panel.background = element_blank()),
                    legend = "none",title = "",xlab = 'BGC Class', ylab = 'BGC region length (kb)',width = 0.5)+
    #scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))+
    scale_y_continuous(limits = c(0, 180), breaks = seq(0, 200, 50)) +
    scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)

    #order
    order_test3 <- data %>%
      group_by(BGC_Class) %>%
      kruskal_test(regionLengthK ~ Group) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj")

    BGC_Class_order3 = as.character(order_test3[order(order_test3$p.adj),]$BGC_Class)
    #leter sign
    df4 <- PMCMR_compare1(data,'BGC_Class','Group','regionLengthK')
    write.table(df4$pvalue,file=paste(opt$outdir,"/",opt$prefix,".BGC_Class.regionLength.nonparametric.pvalue.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)
    write.table(df4$leter,file=paste(opt$outdir,"/",opt$prefix,".BGC_Class.regionLength.nonparametric.leter.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)
    df4$leter$std <- as.numeric(df4$leter$std)
    df4$leter$mean <- as.numeric(df4$leter$mean)
    df4$leter$Max <- as.numeric(df4$leter$Max)
    df4$leter$Group <- factor(df4$leter$Group,levels=group_list)
    df4$leter$BGC_Class <- factor(df4$leter$BGC_Class,levels=BGC_Class_order3)
    df4$leter <- arrange(df4$leter, order(df4$leter$BGC_Class),order(df4$leter$Group))
    splot4 <- df4$leter[order(df4$leter[,'BGC_Class'],df4$leter[,'Group']),]

    Sbox3 <- box3 + geom_text(data=splot4,aes(x=BGC_Class,y=170,label=Letters,colour=Group),angle=90,size = 3,position = position_dodge2(0.9),show.legend = FALSE)+scale_color_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-BGC_Class.regionLength.pdf",sep=""),width = 3.5,height = 2)
    box3
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-BGC_Class.regionLength.png",sep=""),width = 3.5,height = 2,units='in',res=600)
    box3
    dev.off()

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-BGC_Class.regionLength.sign.pdf",sep=""),width = 3.5,height = 2)
    Sbox3
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-BGC_Class.regionLength.sign.png",sep=""),width = 3.5,height = 2,units='in',res=600)
    Sbox3
    dev.off()


#################### @ 3 E ####################


#################### @ 4 S ####################

# small BGC type 3 group  regionLengthK comparing
    Fbox4 <- facet(ggboxplot(removeBGCdata, x = "Group", y = "regionLengthK", fill = "Group",facet.by = "BGC",alpha=0.8,outlier.size = 0.01,size = 0.2,outlier.shape = 20,
            ggtheme = theme_bw() +
              theme(axis.text.x = element_text(color="black",size=10,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                                   axis.text.y = element_text(color="black",size=6,face="bold"),
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
                           legend = "none",title = "",xlab = 'BGC', ylab = 'BGC region length (kb)',width = 0.7)+ 

                   scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
                   scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list), 
        facet.by = "BGC", ncol =7,short.panel.labs = TRUE,strip.position = "top",scales = "free_y", panel.labs.font = list(size = 6, angle = 0),panel.grid = element_blank(), panel.labs.background = list(fill = NA, color = NA) )

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-BGC.regionLength.pdf",sep=""),width = 6,height = 8)
    Fbox4
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-BGC.regionLength.png",sep=""),width = 6,height = 8,units='in',res=600)
    Fbox4
    dev.off()



    #order
    order_test4 <- removeBGCdata %>%
      group_by(BGC) %>%
      kruskal_test(regionLengthK ~ Group) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj")

    BGC_order4 = as.character(order_test4[order(order_test4$p.adj),]$BGC)
    #leter sign
    df5 <- PMCMR_compare1(removeBGCdata,'BGC','Group','regionLengthK')
    write.table(df5$pvalue,file=paste(opt$outdir,"/",opt$prefix,".BGC.regionLength.nonparametric.pvalue.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)
    write.table(df5$leter,file=paste(opt$outdir,"/",opt$prefix,".BGC.regionLength.nonparametric.leter.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)
    df5$leter$std <- as.numeric(df5$leter$std)
    df5$leter$mean <- as.numeric(df5$leter$mean)
    df5$leter$Max <- as.numeric(df5$leter$Max)
    df5$leter$Group <- factor(df5$leter$Group,levels=group_list)
    df5$leter$BGC <- factor(df5$leter$BGC,levels=BGC_order4)
    df5$leter <- arrange(df5$leter, order(df5$leter$BGC),order(df5$leter$Group))
    splot5 <- df5$leter[order(df5$leter[,'BGC'],df5$leter[,'Group']),]

    SFbox4 <- Fbox4 + geom_text(data=splot5,aes(x=Group,y=Max*1.2 ,label=Letters,colour=Group),size = 3,position = position_dodge2(0.9),show.legend = FALSE)+scale_color_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-BGC.regionLength.sign.pdf",sep=""),width = 6,height = 8)
    SFbox4
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-BGC.regionLength.sign.png",sep=""),width = 6,height = 8,units='in',res=600)
    SFbox4
    dev.off()

print('$$$$$$$$$4444444444 DONE')

#################### @ 4 E ####################

#################### @ 5 S ####################

#BGC 3 group comparing
  print("##############")
    Fbox5 <- facet(ggboxplot(removeBGCdata, x = "Group", y = "CDSNum", fill = "Group",facet.by = "BGC",alpha=0.8,outlier.size = 0.01,size = 0.2,outlier.shape = 20,
            ggtheme = theme_bw() +
              theme(axis.text.x = element_text(color="black",size=10,angle=-45,hjust= 0.1 ,vjust = 0 ,face="bold"),
                                   axis.text.y = element_text(color="black",size=6,face="bold"),
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
                           legend = "none",title = "",xlab = 'BGC', ylab = '# of BGC genes',width = 0.7)+ 

                   scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
                   scale_fill_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list), 
        facet.by = "BGC", ncol =7,short.panel.labs = TRUE,strip.position = "top",scales = "free_y", panel.labs.font = list(size = 6, angle = 0),panel.grid = element_blank(), panel.labs.background = list(fill = NA, color = NA) )






    #order
    print('&&&&&&&&&&&&&&&&&&&')

    BGCorder_test <- removeBGCdata %>%
      group_by(BGC) %>%
      kruskal_test(CDSNum ~ Group) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj")
    print('%%%%%%%%%%%%%%%%%%')
    BGC_order = as.character(BGCorder_test[order(BGCorder_test$p.adj),]$BGC)
    #leter sign
    df6 <- PMCMR_compare1(removeBGCdata,'BGC','Group','CDSNum')
    write.table(df6$pvalue,file=paste(opt$outdir,"/",opt$prefix,".BGC.CDSNum.nonparametric.pvalue.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)
    write.table(df6$leter,file=paste(opt$outdir,"/",opt$prefix,".BGC.CDSNum.nonparametric.leter.txt",sep=""),quote = F, sep = "\t", row.names= T, col.names = NA)
    df6$leter$std <- as.numeric(df6$leter$std)
    df6$leter$mean <- as.numeric(df6$leter$mean)
    df6$leter$Max <- as.numeric(df6$leter$Max)
    df6$leter$Group <- factor(df6$leter$Group,levels=group_list)
    df6$leter$BGC <- factor(df6$leter$BGC,levels=BGC_order)
    df6$leter <- arrange(df6$leter, order(df6$leter$BGC),order(df6$leter$Group))
    splot6 <- df6$leter[order(df6$leter[,'BGC'],df6$leter[,'Group']),]
    
    SFbox5 <- Fbox5 + geom_text(data=splot6,aes(x=Group,y=Max*1.2 ,label=Letters,colour=Group),size = 3,position = position_dodge2(0.9),show.legend = FALSE)+scale_color_manual(values=color_var,labels=group_list,limits = group_list, breaks = group_list)

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-BGC.CDSNum.pdf",sep=""),width = 6,height = 8)
    Fbox5
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-BGC.CDSNum.png",sep=""),width = 6,height = 8,units='in',res=600)
    Fbox5
    dev.off()

    pdf(paste(opt$outdir,"/",opt$prefix,".3Group-BGC.CDSNum.sign.pdf",sep=""),width = 6,height = 8)
    SFbox5
    dev.off()
    png(paste(opt$outdir,"/",opt$prefix,".3Group-BGC.CDSNum.sign.png",sep=""),width = 6,height = 8,units='in',res=600)
    SFbox5
    dev.off()

#################### @ 5 E ####################
