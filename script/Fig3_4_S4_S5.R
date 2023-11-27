library ("ComplexHeatmap")
library ("RColorBrewer")
library("circlize")
library('cowplot')
library('grid')
library('dendsort')
Args <- commandArgs(TRUE)


corr <- read.table(Args[1],header=T,check.names=F,row.names=1)
pvalue <- read.table(Args[2],header=T,check.names=F,row.names=1)
enrich <- read.table(Args[3],header=T,check.names=F,row.names=1)
levelInfo <- read.table(Args[4],header=T,check.names=F,row.names=1)
gramlInfo <- read.table(Args[5],header=T,check.names=F,row.names=1)


#sample_order <- Heatmap(as.matrix(corr),show_row_names = T , cluster_rows = TRUE,  cluster_columns = TRUE)
	
bgcOrder <- c("lanthipeptide",  "lassopeptide",  "linaridin",  "LAP",  "microviridin",  "proteusin",  "thiopeptide",  "bacteriocin",  "NRPS",  "NRPS-like",  "T1PKS",  "T3PKS",  "terpene",  "indole",  "acyl_amino_acids",  "arylpolyene",  "betalactone",  "butyrolactone",  "ectoine",  "furan",  "hglE-KS",  "hserlactone",  "ladderane",  "NAGGN",  "phenazine",  "resorcinol",  "siderophore",  "transAT-PKS",  "transAT-PKS-like",  "TfuA-related" )

#bgcOrder <- c( "sactipeptide",  "cyanobactin",  "proteusin",  "microviridin",  "LAP",  "lanthipeptide",  "linaridin",  "bacteriocin",  "thiopeptide",  "bottromycin",  "head_to_tail",  "RaS-RiPP",  "lassopeptide",  "NRPS",  "NRPS-like",  "PKS-like",  "T1PKS",  "T2PKS",  "T3PKS",  "transAT-PKS",  "transAT-PKS-like",  "terpene",  "acyl_amino_acids",  "amglyccycl",  "arylpolyene",  "betalactone",  "butyrolactone",  "CDPS",  "ectoine",  "furan",  "hglE-KS",  "hserlactone",  "indole",  "ladderane",  "NAGGN",  "nucleoside",  "oligosaccharide",  "other",  "PBDE",  "phenazine",  "phosphonate",  "PpyS-KS",  "PUFA",  "resorcinol",  "siderophore",  "TfuA-related",  "tropodithietic-acid",  "phosphoglycolipid")


#levelOrder <- c("Marine:p__Bacteroidota",  "Freshwater:p__Bacteroidota",  "Terrestrial:p__Bacteroidota", 
# "Marine:p__Planctomycetota",  "Freshwater:p__Planctomycetota",  "Terrestrial:p__Planctomycetota",
#   "Marine:p__Proteobacteria",  "Freshwater:p__Proteobacteria",  "Terrestrial:p__Proteobacteria", 
#    "Marine:p__Acidobacteriota",  "Freshwater:p__Acidobacteriota",  "Terrestrial:p__Acidobacteriota", 
#	 "Marine:p__Actinobacteriota",  "Freshwater:p__Actinobacteriota",  "Terrestrial:p__Actinobacteriota", 
#	  "Marine:p__Chloroflexota",  "Freshwater:p__Chloroflexota",  "Terrestrial:p__Chloroflexota",  
#	  "Marine:p__Cyanobacteria",  "Freshwater:p__Cyanobacteria",  "Terrestrial:p__Cyanobacteria",  
#	  "Marine:p__Firmicutes_A",  "Freshwater:p__Firmicutes_A",  "Terrestrial:p__Firmicutes_A",  
#	  "Marine:p__Gemmatimonadota",  "Freshwater:p__Gemmatimonadota",  "Terrestrial:p__Gemmatimonadota", 
#	   "Marine:p__Verrucomicrobiota",  "Freshwater:p__Verrucomicrobiota",  "Terrestrial:p__Verrucomicrobiota"
#	    )


levelOrder <- c("Marine:p__Actinobacteriota",  "Freshwater:p__Actinobacteriota",  "Terrestrial:p__Actinobacteriota",  "Marine:p__Chloroflexota",  "Freshwater:p__Chloroflexota",  "Terrestrial:p__Chloroflexota",   "Marine:p__Firmicutes_A",  "Freshwater:p__Firmicutes_A",  "Terrestrial:p__Firmicutes_A",  "Marine:p__Acidobacteriota",  "Freshwater:p__Acidobacteriota",  "Terrestrial:p__Acidobacteriota",  "Marine:p__Bacteroidota",  "Freshwater:p__Bacteroidota",  "Terrestrial:p__Bacteroidota",  "Marine:p__Cyanobacteria",  "Freshwater:p__Cyanobacteria",  "Terrestrial:p__Cyanobacteria",  "Marine:p__Gemmatimonadota",  "Freshwater:p__Gemmatimonadota",  "Terrestrial:p__Gemmatimonadota",  "Marine:p__Planctomycetota",  "Freshwater:p__Planctomycetota",  "Terrestrial:p__Planctomycetota",  "Marine:p__Proteobacteria",  "Freshwater:p__Proteobacteria",  "Terrestrial:p__Proteobacteria",  "Marine:p__Verrucomicrobiota",  "Freshwater:p__Verrucomicrobiota",  "Terrestrial:p__Verrucomicrobiota")
#colOrder <- colnames(corr)[column_order(sample_order)]
#rowOrder <- rownames(corr)[row_order(sample_order)]



pvalue_order <- as.matrix(pvalue[levelOrder,bgcOrder])
plotdata <- as.matrix(corr[levelOrder,bgcOrder])
enrich_order <- as.matrix(enrich[bgcOrder,])
levelInfo_order <- as.matrix(levelInfo[levelOrder,])
gramInfo_order <- as.matrix(gramlInfo[levelOrder,])

annotation_top <-  HeatmapAnnotation(
	"BGC Class" = enrich_order[,1],
	col=list("BGC Class" = c("RiPPs" = "#e78ac3", "NRPS" = "#8da0cb", "PKS" = "#fc8d62", "Terpene" = "#66c2a5", "Others" = "#cccccc") ), 
	show_legend=TRUE,
	annotation_name_gp = gpar(fontsize = 10),gap=unit(0.1,"cm"),
	annotation_legend_param=list(title_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 8)),
	simple_anno_size = unit(0.2, "cm")
	)

print(gramInfo_order[,1])
annotation_row <- rowAnnotation(
	Group=levelInfo_order[,2],
	Gram = gramInfo_order[,1],
	"log(# OTU)" = anno_barplot(log(as.numeric(levelInfo_order[,3])), width = unit(2.5, "cm"), gp = gpar(fill = "grey50", col = "grey50")),
	col=list(Phylum= c("Bacteroidota"="#b3de69", "Proteobacteria"="#80b1d3", "Planctomycetota"="#7570b3",  "Acidobacteriota" = "#efc86e" , "Actinobacteriota"= "#97c684", "Chloroflexota"= "#6f9969",  "Cyanobacteria"= "#aab5d5", "Firmicutes_A"= "#808fe1",  "Gemmatimonadota"= "#5c66a8", "Verrucomicrobiota"= "#454a74"),
	Group= c("Marine" = "#0f7ba2", "Freshwater" = "#43b284", "Terrestrial" = "#fab255"),
	Gram = c("G+"="#d95f02","G-"="#1b9e77","CPR"="#7570b3","Cand."="#66a61e")

	),
	show_legend=TRUE,
	annotation_name_gp = gpar(fontsize = 8),gap=unit(0.1,"cm"),
	annotation_legend_param=list(title_gp = gpar(fontsize = 10),labels_gp = gpar(fontsize = 10)),
	simple_anno_size = unit(0.2, "cm")
	)

Width <- ncol(corr)*0.5 + 6
Height <- nrow(corr)*0.5 + 6

#col_fun = colorRamp2(c(min(plotdata), 0, max(plotdata)), c("navy", "white", "firebrick3"))

#heat_lgd = Legend(col_fun = col_fun, title = "Abundance",labels_gp=gpar(fontsize=8),title_gp=gpar(fontsize=8),direction = "horizontal",border = "grey50")
#anno_lgd = lgd = Legend(labels = c( "m_DM_0","m_DM1_0","m_DM1_7","m_DM2_0","m_HC","m_HLP_0","m-HLP_18","HLP_0","HLP_1","HLP_4","HLP_7","HLP_16","HLP_18","HC","DM_0","DM_1","DM_4","DM_7","HC_0","HC_1","HC_4","HC_7","preFMT","FMT" ), legend_gp = gpar(fill = color,fontsize = 8), title = "Enrichment", title_gp = gpar(fontsize = 8))
top90Value = max(plotdata)*0.8


#color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100) YlGnBu
print("%%%%%")
f3 <-c("#252525", colorRampPalette(rev(brewer.pal(6, "YlOrRd")))(100))
heatmap <- Heatmap(plotdata,
	top_annotation = annotation_top,
	right_annotation = annotation_row,
	row_order = levelOrder, 
	column_order = bgcOrder,
	#row_dend_reorder = TRUE,
	#column_dend_reorder=TRUE,
	cluster_rows = FALSE,
	cluster_columns = FALSE,
	show_row_names = FALSE,
	col = f3,
	column_split = factor(enrich_order[,1], levels = c("RiPPs","NRPS","PKS","Terpene","Others")) ,
	row_split = factor(levelInfo_order[,1], levels = c("Actinobacteriota",  "Chloroflexota",  "Firmicutes",  "Firmicutes_A",  "Acidobacteriota",  "Bacteroidota",  "Bdellovibrionota",  "Chlamydiota",  "Cyanobacteria",  "Desulfobacterota",  "Desulfobacterota_F",  "Gemmatimonadota",  "Myxococcota",  "Nitrospirota",  "Planctomycetota",  "Proteobacteria",  "Spirochaetota",  "Verrucomicrobiota",  "Omnitrophota",  "Patescibacteria",  "Zixibacteria")),
	row_title_rot = 0,
	column_title_rot = 90,
	column_title_gp = gpar(fontsize=10),
	row_title_gp = gpar(fontsize=10),
	heatmap_legend_param = list(border = "grey50",title = "Relative abundance",labels_gp=gpar(fontsize=8),legend_height = unit(2,"cm"),grid_width= unit(0.3, "cm"),direction = "horizontal",title_gp = gpar(fontsize = 8),legend_width = unit(2.5, "cm"),title_gp=gpar(fontsize=8)),
	width = unit(0.3*ncol(plotdata), "cm"), height = unit(0.3*nrow(plotdata), "cm"),
	row_names_gp = gpar(fontsize = 8),
	column_names_gp = gpar(fontsize = 8),
	cell_fun = function(j, i, x, y, width, height, fill) 
	{
		if (pvalue_order[i,j] == 1 ){
			if (  plotdata[i,j]  >=  top90Value ){
				grid.points(x,y, pch="*",gp=gpar(fontsize = 12,col="black"))
			}else{
				grid.points(x,y, pch="*",gp=gpar(fontsize = 12,col="white"))
			}
			}
	}
)
pdf(Args[6],width=10,height=10)

draw(heatmap, merge_legend = TRUE)

