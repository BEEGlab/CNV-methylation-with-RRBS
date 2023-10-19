##R script
### 5. Differentially Methylated Sites
#input files: Samples metadata, global and pairwise DMS results, Methylation estimates at filtered CpG sites and the order of the samples in the file, gene CNVs, DMS overlapping genes and CNVs
#output files: DMS information reformatted. Figures: DMS PCA and genome scan

library(RColorBrewer)
library(dplyr)
library(stringr)
library(plyr)
library(readr)
library(ggplot2)
library(reshape)
library(pheatmap)

###############################################################################################################
###############################################################################################################
samples<-read.table("Samples96.out",col.names=c("sample","sex","pop")) 
#order
poporder<-c("LET","BAR","NYN","FAL","KIE","SYL") ###order by increasing PSU when plotting
sampleorder<-samples$sample[order(match(samples$pop,poporder),samples$sex)]
samples$sample<-factor(samples$sample,levels=sampleorder)
###############################################################################################################
###############################################################################################################

##import DMS
dir_dms<-"."
files_dms <- dir(dir_dms,pattern="*bed$") ##files that end with ...
files_dms <- files_dms[grepl("myDiff15p",files_dms)] ##files that contain
DMS<-data.frame()
print(files_dms)
metadata<-data.frame(comp=0:15,pop=c("global","NYN_SYL","KIE_SYL","FAL_SYL","BAR_SYL","LET_SYL","NYN_KIE","FAL_KIE","BAR_KIE","LET_KIE","FAL_NYN","BAR_NYN","LET_NYN","BAR_FAL","LET_FAL","LET_BAR"))

for(i in files_dms){
  if(grepl("myDiff",i)){
    samplename<-strsplit( strsplit(i,"_")[[1]][2] , "\\.")[[1]][1] ##extracts name from filename
    if(samplename=="global"){
      samplename=0
    }
    data<-read.delim(paste0(dir_dms,i))
    data$comp<-samplename
    if(empty(DMS)){
      DMS<-data
    }else{
      DMS<-rbind(DMS,data)
    }
  }
}
DMS$comp<-as.numeric(DMS$comp)
DMS$comp<-factor(DMS$comp,levels=sort(unique(DMS$comp)))##make factor and order
DMS25<-(DMS[(DMS$scores.meth.diff>25|DMS$scores.meth.diff<(-25)),])

################################################################################
## PCA of these DMS sites

##already filtered for coverage of 10x
meth0<-read.csv("MethylationProportions.txt") ##all data!!
names(meth0)[9]<-"chrom"
methsample<-read.table("MethylationSampleOrder.txt",sep=" ")
##remove decimals from transcript IDs
meth0$feature.name<-gsub("(.*)\\..*","\\1",meth0$feature.name)

temp<-meth0[,1:11]
methcov<-meth0[,1:11]
##calculate C/T ratios per sample.
offset2=3
for(n in 1:length(methsample)){
  ##offset by 13 in columns
  off=11
  ##increment by 3 each time in the loop
  inc<-off+(n*offset2)
  temp[,(n+11)]<-meth0[,inc]/meth0[,(inc-1)]
  methcov[,(n+11)]<-meth0[,(inc-1)]
}
names(temp)[12:ncol(temp)]<-as.character(methsample[1,])
names(methcov)[12:ncol(methcov)]<-as.character(methsample[1,])

mat<-temp[,9:ncol(temp)]
names(mat)[1:3]<-c("seqnames","starts","ends")
mat[,2]<-mat[,2]-1
# remove rows containing NA values, they might be introduced at unite step
mat=mat[ rowSums(is.na(mat))==0, ]

##get only DMS regions
mat_dms<-merge(mat,DMS)
##get rid of duplicates
mat_dms<-unique(mat_dms[,1:90])

##only keep methylation values
mat_dms<-mat_dms[,4:90]
mat<-mat[,4:90]

################################################################################
samples2<-join(samples,data.frame(sample=names(mat)),type="inner")
popsorder<-as.data.frame(samples2$pop)
row.names(popsorder)<-colnames(mat)
colnames(popsorder)<-"pop"
popsorder$pop<-factor(popsorder$pop,levels=poporder)

################
library(ggfortify)
pca_res <- prcomp(t(mat_dms), scale. = TRUE)
#pca_res <- prcomp(t(mat_dms))
mat_dms_info<-cbind(popsorder,t(mat_dms))
theme_set(theme_light(16))
autoplot(pca_res, data=mat_dms_info,colour='pop') +
  geom_text(aes(label=rownames(mat_dms_info),color=pop,), size=3,
            nudge_y = 0.01,
            check_overlap = T) +
  theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank())
ggsave(paste0("DMS.Methylation.PCA.pdf"))

### threshold 25
DMS<-DMS25
###

##merge data
##same as: cat DMS.bed | cut -f 1-3,9 | sortBed -i -  | mergeBed -i - -c 4 -o count,distinct_sort_num > DMS.merged.bed
DMS_num<-aggregate(DMS[,c(1:3,9)],comp~seqnames+starts+ends,length)
names(DMS_num)[4]<-"num"
DMS2<-aggregate(DMS[,c(1:3,9)],comp~seqnames+starts+ends,function(x)paste(sort(x), collapse = ','))
DMS2<-merge(DMS_num,DMS2)
DMS2<-DMS2[order(DMS2$seqnames,DMS2$starts),]
write.table(DMS2,"DMS.merged.bed",sep="\t",row.names=F,col.names=T,quote=F)
###############################################################################################################
cnvs<-read.csv('CNVs.genotypes.combined.CN.csv')
write.table(cnvs[,c(1:3,100:101,104:115)],"CNVs.genotypes.combined.CN.bed",sep="\t",row.names=F,col.names=T,quote=F)
genes<-read.csv('CNVs.genes.combined.CN.csv')
write.table(genes[,c(103:105,4,101,102,108,111:126)],"CNVs.genes.combined.CN.bed",sep="\t",row.names=F,col.names=T,quote=F)

##use bedtools to look at overlaps
#sed -i.bak 's/^X\./#/' CNVs.genotypes.combined.CN.bed
#sed -i.bak 's/^chrom_gene/#Chrom/' CNVs.genes.combined.CN.bed
#sed -i.bak 's/^seqnames/#Chrom/' DMS.merged.bed

#intersectBed -a DMS.merged.bed -b CNVs.genotypes.combined.CN.bed -wao > DMS.CNVs.bed
#intersectBed -a DMS.merged.bed -b CNVs.genes.combined.CN.bed -wao > DMS.CNVs.genes.bed
#intersectBed -a DMS.merged.bed -b CNVs.genotypes.combined.CN.bed | wc -l

#cat CNVs.genes.combined.CN.bed | sed 1,1d | awk '{$2=$2-1500;OFS="\t";print $1,$2,$3,$4}' > CNVs.genes.combined.CN.bed.promoters
#intersectBed -a DMS.merged.bed -b CNVs.genes.combined.CN.bed.promoters -wao > DMS.CNVs.genes.bed.promoters
#intersectBed -a DMS.merged.bed -b CNVs.genes.combined.CN.bed.promoters | wc -l

dmsg<-read.delim("DMS.CNVs.genes.bed",header=F)
names(dmsg)<-c(names(DMS2),names(genes[,c(103:105,4,101,102,108,111:126)]),"cnv")

dms<-read.delim("DMS.CNVs.bed",header=F)
names(dms)<-c(names(DMS2),names(cnvs[,c(1:3,100:101,104:115)]),"cnv")

##randomize in bedtools to count overlaps, determine whether observed is within range of random permutations
#shuffleBed -i CNVs.genotypes.combined.CN.bed -g gasAcu1.chrom.sizes |intersectBed -a DMS.merged.bed -b - | wc -l

#for i in `seq 1 1000`;do 
#c=`shuffleBed -i CNVs.genotypes.combined.CN.bed -g gasAcu1.chrom.sizes |intersectBed -a DMS.merged.bed -b - | wc -l`;
#echo $c >> DMS.CNVs.bed.random;
#done
##also for genes
#for i in `seq 1 1000`;do c=`shuffleBed -i CNVs.genes.combined.CN.bed.promoters -g gasAcu1.chrom.sizes |intersectBed -a DMS.merged.bed -b - | wc -l`; echo $c >> DMS.CNVs.genes.bed.random;done

##################################################################################################################
##################################################################################################################
### DMS scans
chrom_order<-c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX",
               "chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX",
               "chrXXI","chrUn")
chrom_order2<-gsub("chr","",chrom_order)
poporder<-c("LET","BAR","NYN","FAL","KIE","SYL") ###order by increasing PSU when plotting

scandms<-DMS[,c(1:3,7,9)]
names(scandms)[1:3]<-c("Chrom","Start","End")
scandms[,1]<-gsub("chr","",scandms[,1])
scandms$Chrom<-factor(scandms$Chrom,levels=chrom_order2)
scandms$comp<-factor(scandms$comp)

## gets all pairwise from each pop into one, make sure to swap signal depending on order of analysis
test<-scandms[scandms$comp %in% c(1,2,3,4,5),]
test$scores.meth.diff<-test$scores.meth.diff*(-1)
test$comp<-"SYL"
avgs<-test
test<-scandms[scandms$comp %in% c(2,6,7,8,9),]
test$scores.meth.diff[test$comp %in% c(6,7,8,9)]<-test$scores.meth.diff[test$comp %in% c(6,7,8,9)]*(-1)
test$comp<-"KIE"
avgs<-rbind(avgs,test)
test<-scandms[scandms$comp %in% c(3,7,10,13,14),]
test$scores.meth.diff[test$comp %in% c(13,14)]<-test$scores.meth.diff[test$comp %in% c(13,14)]*(-1)
test$comp<-"FAL"
avgs<-rbind(avgs,test)
test<-scandms[scandms$comp %in% c(1,6,10,11,12),]
test$scores.meth.diff[test$comp %in% c(10,11,12)]<-test$scores.meth.diff[test$comp %in% c(10,11,12)]*(-1)
test$comp<-"NYN"
avgs<-rbind(avgs,test)
test<-scandms[scandms$comp %in% c(4,8,11,13,15),]
test$scores.meth.diff[test$comp %in% c(15)]<-test$scores.meth.diff[test$comp %in% c(15)]*(-1)
test$comp<-"BAR"
avgs<-rbind(avgs,test)
test<-scandms[scandms$comp %in% c(5,9,12,14,15),]
test$comp<-"LET"
avgs<-rbind(avgs,test)
avgs$comp<-factor(avgs$comp,levels=poporder)

################################################################################
#Default ggplot colors
gg_color_hue <- function(n) {  hues = seq(15, 375, length = n + 1); hcl(h = hues, l = 65, c = 100)[1:n] }
n = 6
cols = gg_color_hue(n)
cols2 = c("#FFC5C1","#E1D68E","#A8EABC","#A0F2F5","#BFD7FF","#FFC3F8") #faded
#LET BAR NYN FAL KIE SYL

avgs$num <- as.numeric(avgs$comp)

avgs2<-avgs[abs(avgs$scores.meth.diff)>=25,] ##add color, depending on sign
avgs2$color1 <- cols[avgs2$num]
avgs2$color2 <- cols2[avgs2$num]
avgs2$color <- ifelse(avgs2$scores.meth.diff>0,avgs2$color1,avgs2$color2)

avgs<-merge(avgs,avgs2[,c(1:6,9)],all.x=T)
avgs$color <- ifelse(is.na(avgs$color),"#F0EDEC",avgs$color)
avgs$color <- factor(avgs$color)

ggplot(avgs,aes(x=Start,y=scores.meth.diff)) +
  geom_point(aes(color=color),size=0.15) +
  scale_color_manual(values = levels(avgs$color), 
                     labels = levels(avgs$color)) +
  facet_grid(comp~ Chrom,scales="free_x", space="free_x") + # separate plots for each chromosome, stretch smaller chromosomes
  xlab("Chrom") +
  ylab("DMS") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.margin = unit(0.05, "lines"))
ggsave("DMS.scan.pdf",width=8,height=6)