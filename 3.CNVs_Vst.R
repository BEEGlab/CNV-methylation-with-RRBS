##R script
### 3. Calculates Vst
#input files: Samples metadata, sex-specific CNVs, CNV genotypes (rounded)
#output files: Vst estimates. Figures: Vst genome scans and heatmap

library(RColorBrewer)
library(dplyr)
library(stringr)
library(plyr)
library(ggplot2)
library(reshape)
library(pheatmap)

###############################################################################################################
samples<-read.table("Samples96.out",col.names=c("sample","sex","pop"))
#order
poporder<-c("LET","BAR","NYN","FAL","KIE","SYL") ###order by increasing PSU when plotting
sampleorder<-samples$sample[order(match(samples$pop,poporder),samples$sex)]
samples$sample<-factor(samples$sample,levels=sampleorder)

###############################################################################################################
## exclude
SexSpecific<-read.csv("SexSpecificCNVs.csv")
SexSpecific<-SexSpecific[,2:4]

###############################################################################################################
###############################################################################################################
getVst <- function(dat, groups, comparison) {
  groupLevels <- levels(groups)
  dat1 <- na.omit(dat[groups==groupLevels[groupLevels==comparison[1]]])
  dat2 <- na.omit(dat[groups==groupLevels[groupLevels==comparison[2]]])
  Vtotal <- var(c(dat1, dat2))
  Vgroup <- ((var(dat1)*length(dat1)) + (var(dat2)*length(dat2))) /
    (length(dat1)+length(dat2))
  Vst <- c((Vtotal-Vgroup) / Vtotal)
  if (Vst == "NaN" | Vst < 0){
    Vst <- 0
  }
  return(Vst)
}
###############################################################################################################
###############################################################################################################

## CNV Regions
vstregions<-read.csv("CNVs.genotypes.combined.CN.csv")
vstregions<-anti_join(vstregions,SexSpecific)

dd <- vstregions[,4:99] ##just read depth
group <- factor(samples$pop)

##Vst pairwise
vstregions$Vst_SvsL <- apply(dd, 1, function(x) getVst(x, group, c("SYL","LET")))
vstregions$Vst_SvsB <- apply(dd, 1, function(x) getVst(x, group, c("SYL","BAR")))
vstregions$Vst_SvsN <- apply(dd, 1, function(x) getVst(x, group, c("SYL","NYN")))
vstregions$Vst_SvsF <- apply(dd, 1, function(x) getVst(x, group, c("SYL","FAL")))
vstregions$Vst_SvsK <- apply(dd, 1, function(x) getVst(x, group, c("SYL","KIE")))

vstregions$Vst_KvsL <- apply(dd, 1, function(x) getVst(x, group, c("KIE","LET")))
vstregions$Vst_KvsB <- apply(dd, 1, function(x) getVst(x, group, c("KIE","BAR")))
vstregions$Vst_KvsN <- apply(dd, 1, function(x) getVst(x, group, c("KIE","NYN")))
vstregions$Vst_KvsF <- apply(dd, 1, function(x) getVst(x, group, c("KIE","FAL")))

vstregions$Vst_FvsL <- apply(dd, 1, function(x) getVst(x, group, c("FAL","LET")))
vstregions$Vst_FvsB <- apply(dd, 1, function(x) getVst(x, group, c("FAL","BAR")))
vstregions$Vst_FvsN <- apply(dd, 1, function(x) getVst(x, group, c("FAL","NYN")))

vstregions$Vst_NvsL <- apply(dd, 1, function(x) getVst(x, group, c("NYN","LET")))
vstregions$Vst_NvsB <- apply(dd, 1, function(x) getVst(x, group, c("NYN","BAR")))

vstregions$Vst_BvsL <- apply(dd, 1, function(x) getVst(x, group, c("BAR","LET")))

###############################################################################################################
## CNV Genes
vstgenes<-read.csv("CNVs.genes.combined.CN.csv") ##already excludes sex-specific

dd <- vstgenes[,5:100] ##just read depth
group <- factor(samples$pop)

##Vst pairwise
vstgenes$Vst_SvsL <- apply(dd, 1, function(x) getVst(x, group, c("SYL","LET")))
vstgenes$Vst_SvsB <- apply(dd, 1, function(x) getVst(x, group, c("SYL","BAR")))
vstgenes$Vst_SvsN <- apply(dd, 1, function(x) getVst(x, group, c("SYL","NYN")))
vstgenes$Vst_SvsF <- apply(dd, 1, function(x) getVst(x, group, c("SYL","FAL")))
vstgenes$Vst_SvsK <- apply(dd, 1, function(x) getVst(x, group, c("SYL","KIE")))

vstgenes$Vst_KvsL <- apply(dd, 1, function(x) getVst(x, group, c("KIE","LET")))
vstgenes$Vst_KvsB <- apply(dd, 1, function(x) getVst(x, group, c("KIE","BAR")))
vstgenes$Vst_KvsN <- apply(dd, 1, function(x) getVst(x, group, c("KIE","NYN")))
vstgenes$Vst_KvsF <- apply(dd, 1, function(x) getVst(x, group, c("KIE","FAL")))

vstgenes$Vst_FvsL <- apply(dd, 1, function(x) getVst(x, group, c("FAL","LET")))
vstgenes$Vst_FvsB <- apply(dd, 1, function(x) getVst(x, group, c("FAL","BAR")))
vstgenes$Vst_FvsN <- apply(dd, 1, function(x) getVst(x, group, c("FAL","NYN")))

vstgenes$Vst_NvsL <- apply(dd, 1, function(x) getVst(x, group, c("NYN","LET")))
vstgenes$Vst_NvsB <- apply(dd, 1, function(x) getVst(x, group, c("NYN","BAR")))

vstgenes$Vst_BvsL <- apply(dd, 1, function(x) getVst(x, group, c("BAR","LET")))

###plot
chrom_order<-c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX",
               "chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX",
               "chrXXI","chrUn")
chrom_order2<-gsub("chr","",chrom_order)

pvst<-vstregions[,c(1:3,118:(ncol(vstregions)-6))]
names(pvst)[1]<-"Chrom"
pvst[,1]<-gsub("chr","",pvst[,1])
names(pvst)[4:ncol(pvst)]<-gsub("Vst_","",names(pvst[4:ncol(pvst)]))
pvst<-melt(pvst,id.vars=c("Chrom","Start","End"))
pvst$Chrom<-factor(pvst$Chrom,levels=chrom_order2)

ggplot(pvst) + 
  geom_smooth(aes(x=Start,y=value,color=variable),size=0.5,se=F) +
  geom_point(aes(x=Start,y=value,color=variable),size=0.75) +
  facet_grid(variable~ Chrom,scales="free_x", space="free_x") + # separate plots for each chromosome, stretch smaller chromosomes
  xlab("Chrom") + 
  ylab("Vst") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.margin = unit(0.05, "lines"))
ggsave("Vst.regions.pairwise.scan.pdf",width=8,height=6)

pvstg<-vstgenes[,c(1:3,127:(ncol(vstgenes)-6))]
pvstg[,1]<-gsub("chr","",pvstg[,1])
names(pvstg)[1:3]<-c("Chrom","Start","End")
names(pvstg)[4:ncol(pvstg)]<-gsub("Vst_","",names(pvstg[4:ncol(pvstg)]))
pvstg<-melt(pvstg,id.vars=c("Chrom","Start","End"))
pvstg$Chrom<-factor(pvstg$Chrom,levels=chrom_order2)

ggplot(pvstg) + 
  geom_smooth(aes(x=Start,y=value,color=variable),size=0.5,se=F) +
  geom_point(aes(x=Start,y=value,color=variable),size=0.75) +
  #facet_wrap(~ Chrom,scales="free_x",nrow=1) + # separate plots for each chromosome, stretch smaller chromosomes
  facet_grid(variable~ Chrom,scales="free_x", space="free_x") + # separate plots for each chromosome, stretch smaller chromosomes
  xlab("Chrom") +
  #ylim(0,0.45) + 
  ylab("Vst") + 
  #facet_grid(~Chromosome, scales = 'free_x', space = 'free_x', switch = 'x') + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.margin = unit(0.05, "lines"))
ggsave("Vst.genes.pairwise.scan.pdf",width=8,height=6)

write.table(pvst, 'CNVs.genotypes.combined.CN.Vst.pairwise.csv', row.names = F, col.names=F,quote=F,sep="\t")
write.table(pvstg, 'CNVs.genes.combined.CN.Vst.pairwise.csv', row.names = F,col.names=F,quote=F,sep="\t")

############# ############# ############# ############# ############# ############# 
############# ############# ############# ############# ############# ############# 
## Plot average Vst vs 5 other pops
vstregions$LET<-apply(vstregions[,grepl("vsL",names(vstregions))], 1, mean)
vstregions$BAR<-apply(vstregions[,grep(paste(c("vsB","Bvs"),collapse="|"),names(vstregions))], 1, mean)
vstregions$NYN<-apply(vstregions[,grep(paste(c("vsN","Nvs"),collapse="|"),names(vstregions))], 1, mean)
vstregions$FAL<-apply(vstregions[,grep(paste(c("vsF","Fvs"),collapse="|"),names(vstregions))], 1, mean)
vstregions$KIE<-apply(vstregions[,grep(paste(c("vsK","Kvs"),collapse="|"),names(vstregions))], 1, mean)
vstregions$SYL<-apply(vstregions[,grepl("Svs",names(vstregions))], 1, mean)

pvst<-vstregions[,c(1:3,(ncol(vstregions)-5):ncol(vstregions))]
names(pvst)[1]<-"Chrom"
pvst[,1]<-gsub("chr","",pvst[,1])
names(pvst)[4:ncol(pvst)]<-gsub("Vst_","",names(pvst[4:ncol(pvst)]))
pvst<-melt(pvst,id.vars=c("Chrom","Start","End"))
pvst$Chrom<-factor(pvst$Chrom,levels=chrom_order2)

ggplot(pvst) + 
  geom_smooth(aes(x=Start,y=value,color=variable),size=0.5,se=F) +
  geom_point(aes(x=Start,y=value,color=variable),size=0.75) +
  #facet_wrap(~ Chrom,scales="free_x",nrow=1) + # separate plots for each chromosome, stretch smaller chromosomes
  facet_grid(variable~ Chrom,scales="free_x", space="free_x") + # separate plots for each chromosome, stretch smaller chromosomes
  xlab("Chrom") +
  ylab("Vst") + 
  #facet_grid(~Chromosome, scales = 'free_x', space = 'free_x', switch = 'x') + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.margin = unit(0.05, "lines"))
ggsave("Vst.regions.scan.pdf",width=8,height=6)

vstgenes$LET<-apply(vstgenes[,grepl("vsL",names(vstgenes))], 1, mean)
vstgenes$BAR<-apply(vstgenes[,grep(paste(c("vsB","Bvs"),collapse="|"),names(vstgenes))], 1, mean)
vstgenes$NYN<-apply(vstgenes[,grep(paste(c("vsN","Nvs"),collapse="|"),names(vstgenes))], 1, mean)
vstgenes$FAL<-apply(vstgenes[,grep(paste(c("vsF","Fvs"),collapse="|"),names(vstgenes))], 1, mean)
vstgenes$KIE<-apply(vstgenes[,grep(paste(c("vsK","Kvs"),collapse="|"),names(vstgenes))], 1, mean)
vstgenes$SYL<-apply(vstgenes[,grepl("Svs",names(vstgenes))], 1, mean)

pvstg<-vstgenes[,c(1:3,(ncol(vstgenes)-5):ncol(vstgenes))]
pvstg[,1]<-gsub("chr","",pvstg[,1])
names(pvstg)[1:3]<-c("Chrom","Start","End")
names(pvstg)[4:ncol(pvstg)]<-gsub("Vst_","",names(pvstg[4:ncol(pvstg)]))
pvstg<-melt(pvstg,id.vars=c("Chrom","Start","End"))
pvstg$Chrom<-factor(pvstg$Chrom,levels=chrom_order2)

ggplot(pvstg) + 
  geom_smooth(aes(x=Start,y=value,color=variable),size=0.5,se=F) +
  geom_point(aes(x=Start,y=value,color=variable),size=0.75) +
  #facet_wrap(~ Chrom,scales="free_x",nrow=1) + # separate plots for each chromosome, stretch smaller chromosomes
  facet_grid(variable~ Chrom,scales="free_x", space="free_x") + # separate plots for each chromosome, stretch smaller chromosomes
  xlab("Chrom") +
  ylim(0,0.45) + 
  ylab("Vst") + 
  #facet_grid(~Chromosome, scales = 'free_x', space = 'free_x', switch = 'x') + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.margin = unit(0.05, "lines"))
ggsave("Vst.genes.scan.pdf",width=8,height=6)

write.csv(vstregions, 'CNVs.genotypes.combined.CN.Vst.csv', row.names = F)
write.csv(vstgenes, 'CNVs.genes.combined.CN.Vst.csv', row.names = F)

################################################################################
################################################################################
################################################################################
## Heatmap of Vst outliers, top 5% of each pairwise
apply(vstregions[,118:ncol(vstregions)],2,function(x) quantile(x,0.99))
apply(vstgenes[,127:ncol(vstgenes)],2,function(x) quantile(x,0.99))

### Keeps columns where any is > 95th quantile. Useful for now doing heatmaps
vstr_o<-vstregions[,118:ncol(vstregions)] %>%
  dplyr::filter(if_any(everything(), ~.x >= quantile(.x,.99)))
vstg_o<-vstgenes[,127:ncol(vstgenes)] %>%
  dplyr::filter(if_any(everything(), ~.x >= quantile(.x,.98) & .x > 0))

### Format for heatmap
vstr_o<-merge(vstr_o,vstregions)
vstg_o<-merge(vstg_o,vstgenes)

###group for plotting
popgroup<-dplyr::select(samples, pop)
row.names(popgroup)<-samples$sample
popgroup$pop<-factor(popgroup$pop,levels=poporder)

###Heatmap
myColorsO <- brewer.pal(11,"RdBu")

############################################################################
###remove repeated
## regions
vstr_o$id<-paste(vstr_o$X.Chrom,vstr_o$Start,vstr_o$End,sep="_")
vstr_m<-vstr_o[!duplicated(vstr_o$id),25:120]
row.names(vstr_m)<-vstr_o$id[!duplicated(vstr_o$id)]
vstr_m<-vstr_m[levels(sampleorder)]

## genes
vstg_m<-vstg_o[!duplicated(vstg_o$Gene),26:121]
row.names(vstg_m)<-vstg_o$Gene[!duplicated(vstg_o$Gene)]
vstg_m<-vstg_m[levels(sampleorder)]

myColors<-myColorsO[c(1,5,6,7,8,10,11)] ##7 depths, with scale="none"
pheatmap(vstr_m,
         filename = "Heatmap.CN.Vstoutliers.1pct.regions.pdf",
         color = myColors,
         cluster_rows = T,
         show_rownames = T,
         annotation_col = popgroup,
         scale = "none",
         gaps_col=seq(16, 96, by=16),
         cluster_cols=F
)
############################################################################
## GENES
myColors<-myColorsO[c(1,5,6,7,8,10)] ##6 depths, with scale="none"

anno_colors = list(pop = c(LET="#F8766D",BAR="#B79F00",NYN="#00BA38",FAL="#00BFC4",KIE="#619CFF",SYL="#F564E3"))
pheatmap(vstg_m,
         filename = "Heatmap.CN.Vstoutliers.2pct.genes.pdf",
         color = myColors,
         cluster_rows = T,
         show_rownames = T,
         annotation_col = popgroup,
         annotation_colors = anno_colors,
         scale = "none",
         gaps_col=seq(16, 96, by=16),
        cellheight = 6,
         cluster_cols=F,
fontsize=6,
width=10,height=8,
border_color = "white"
)
