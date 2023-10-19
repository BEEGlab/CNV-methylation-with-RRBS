##R script
###Gene ontology enrichment of lists of genes
#input files: list of gene IDs of interest in the first column, and full list of genes with associated gene ontologies
#output files: GO enrichment analysis


##########################################################################################
genesIDs<-read.csv("genes.list.csv") ##list of gene IDs (in the first column) to look for enrichment
ensemblgo<-read.delim("Ensembl110_Stickleback_GO.txt",header=T,sep="\t")
##########################################################################################
# Gene Ontology Analysis
library(topGO)
db<-"Genes.IDs.db" ##Universe db file from biomart

# read custom annotation, all GO IDs separated by commas
GOmap <- readMappings(db)
refset <- names(GOmap) # save the geneID in the universe
  
  #create vector describing the interesting genes with 1, the rest with 0
  genes_of_interest = factor(as.integer(refset %in% geneIDs[,1]))
  names(genes_of_interest) <- refset
  
  #Set the gene ontology categories to loop through
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)){
    ##create the topGO data object
    tgData = new("topGOdata", ontology = ontology[i], allGenes = genes_of_interest, annot = annFUN.gene2GO, gene2GO = GOmap)
    ##sets statistics 
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method="fdr")
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method="fdr")
    allRes    = GenTable(tgData, classic=fisherRes, weight=weightRes, orderBy="weight", ranksOf="classic", topNodes=100, numChar=200)
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(db,"topGO.over",ontology[i],"csv",sep="."))
  }
