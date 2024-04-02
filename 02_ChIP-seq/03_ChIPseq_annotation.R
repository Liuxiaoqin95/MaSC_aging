###annotate peaks

library(clusterProfiler)
library(ChIPseeker)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggplot2)
k=keys(org.Mm.eg.db,keytype = "ENSEMBL")
list=select(org.Mm.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- promoters(txdb, upstream = 3000, downstream = 3000)

filename=list.files(pattern = '.narrowPeak')
peaks=list()
for(i in filename){
  temp=readPeakFile(i)
  peaks[[i]]=temp
  
}
names(peaks)=sapply(filename,function(x){unlist(strsplit(x,'[.]'))[1]})

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=10000,annoDb="org.Mm.eg.db")
plotAnnoPie(peakAnnoList[[1]])
plotAnnoPie(peakAnnoList[[2]])

peakAnnoList_df=lapply(peakAnnoList,as.data.frame)
gene_promoter = lapply(peakAnnoList_df, function(x){inter=abs(x$distanceToTSS)<3000;x$geneId[inter]})
gene_enhancer = lapply(peakAnnoList_df, function(x){inter1=which(abs(x$distanceToTSS)<3000);inter2=setdiff(which(abs(x$distanceToTSS)<10000),inter1);x$geneId[inter2]})

gene_promoter_df= lapply(peakAnnoList_df, function(x){inter=abs(x$distanceToTSS)<3000;x[inter,]})
gene_enhancer_df = lapply(peakAnnoList_df, function(x){inter1=which(abs(x$distanceToTSS)<3000);inter2=setdiff(which(abs(x$distanceToTSS)<10000),inter1);x[inter2,]})

for(i in names(gene_promoter)){
  print(i)
  xlsx::write.xlsx(gene_enhancer_df[[i]],file = "bcl11bbindingsite20201117.xlsx",sheetName = paste0(i,"_enhancer"),append = T,col.names = T)
  xlsx::write.xlsx(gene_promoter_df[[i]],file = "bcl11bbindingsite20201117.xlsx",sheetName = paste0(i,"_promoter"),append = T,col.names = T)
}


