###cell cycle analysis

library(dplyr)
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(monocle)
library(RColorBrewer)
library(Seurat)
library(colorRamps)
library(pheatmap)
library(stringr)
library(GSVA)
library(ComplexHeatmap)
library(dendsort)
mm10_gene_length=read.table("E:/databases/mm10_geneloc/mm10_geneloc_rmdp_exonlength",as.is = T,header = T)

library(scater)
seuratData@meta.data$stage=pData(HSMM)$stage[match(colnames(seuratData),colnames(HSMM))]
table(seuratData@meta.data$stage)/ncol(seuratData)
gene_threshold=1000
seuratData=seuratData
seuratData=seuratData[,seuratData@meta.data$nFeature_RNA>gene_threshold]
seuratData=seuratData[rownames(seuratData)%in%mm10_gene_length$symbol,]#####
seuratData@meta.data$stage=pData(HSMM)$stage[match(colnames(seuratData),colnames(HSMM))]
table(seuratData@meta.data$stage)/ncol(seuratData)

tpm<- log2(scater::calculateTPM(seuratData@assays$RNA@counts,
                                mm10_gene_length$exonlength[match(rownames(seuratData),mm10_gene_length$sym)]) + 1)

cellcycle_Hela=read.table("E:/course/paper/from20200601/technique/cell_cycle/2015.genomeresearch.HSCagingcellcycle/2002. Whitfield et al. 2002/CellCycleGeneList_1134.txt",header = T,sep = "\t",as.is = T)
cellcycle_Hela=cellcycle_Hela[,c("PHASE","NAME","Annotation","Name_Annot")]
cellcycle_Hela$symbol=sapply(cellcycle_Hela$NAME,function(x){unlist(strsplit(x,"[ ]"))[1]})
cellcycle_Hela$mouse=human2mice$m[match(cellcycle_Hela$symbol,human2mice$h)]
cellcycle_Hela=cellcycle_Hela[!is.na(cellcycle_Hela$mouse),]

cellcycle_Hela_list=lapply(c("G1/S","G2","G2/M","M/G1","S phase"),function(x){cellcycle_Hela$mouse[cellcycle_Hela$PHASE==x]})
names(cellcycle_Hela_list)=c("G1/S","G2","G2/M","M/G1","S phase")
cellcycle_Hela_list=lapply(cellcycle_Hela_list,function(x){x[x%in%rownames(tpm)]})
cellcycle_Hela_list=lapply(cellcycle_Hela_list,unique)
cellcycle_Hela_list=lapply(cellcycle_Hela_list,function(x){x[rowSums(tpm[x,])>0]})
cellcycle_Hela_list=lapply(cellcycle_Hela_list,function(x){a=sapply(x,function(y){sum(tpm[y,]!=0)});return(x[a>10])})
cellcycle_Hela_list0=cellcycle_Hela_list


cor_threshold=0.25
print(cor_threshold)

cellcycle_Hela_list=cellcycle_Hela_list0
for(stage in names(cellcycle_Hela_list)){
  print(stage)
  a=apply(tpm[cellcycle_Hela_list[[stage]],],2,mean)
  m=psych::corr.test(a,t(as.matrix(tpm[cellcycle_Hela_list[[stage]],])),adjust = 'fdr',ci=F)
  # hist(m$r)
  cellcycle_Hela_list[[stage]]=colnames(m$r)[m$r[1,]>cor_threshold]
}
t=sapply(cellcycle_Hela_list,length)

if(t["G1/S"]==0|t["G2/M"]==0){next()}
cellcycle_Hela_list_exp=lapply(cellcycle_Hela_list,function(x){t=as.matrix(tpm[x,]);t=t(scale(t(t),scale=F));apply(t,2,mean)})
cellcycle_Hela_list_exp_df=do.call("cbind",cellcycle_Hela_list_exp)
colnames(cellcycle_Hela_list_exp_df)=names(cellcycle_Hela_list_exp)
cellcycle_Hela_list_exp_df=cellcycle_Hela_list_exp_df[,c("G1/S","S phase","G2","G2/M","M/G1")]

df=data.frame(cellcycle_Hela_list_exp_df)
df$cellPhase=seuratData@meta.data$Phase[match(rownames(df),colnames(seuratData))]

p1=ggplot(df,aes(G1.S,G2.M,col=cellPhase))+geom_point()
p1

df$cellPhase2=NA
df$cellPhase2[(df$G1.S<0)&(df$G2.M<0)]="G0"
a=apply(df[is.na(df$cellPhase2),1:5],1,which.max)

df$cellPhase2[is.na(df$cellPhase2)]=colnames(df)[a]
p2.1=ggplot(df,aes(G1.S,G2.M,col=cellPhase2))+geom_point()
a=table(df$cellPhase2)/nrow(df)

# text(b,a/2,labels = round(a,2))
seuratData@meta.data$G0=df$cellPhase2[match(colnames(seuratData),rownames(df))]
data=as.data.frame(table(seuratData@meta.data[,c("G0","stage")]))
data=data%>%group_by(stage)%>%mutate(percent=Freq/sum(Freq))
p2=ggplot(data,aes(stage,Freq,fill=G0))+geom_bar(stat="identity")
p3=ggplot(data,aes(stage,percent,fill=G0))+geom_bar(stat="identity")
p4=ggplot(data)+geom_line(aes(x=stage,y=percent,group=G0,col=G0))+geom_text(aes(x=State2,y=percent,label=round(percent,3)))

t=table(seuratData@meta.data[c('stage',"G0")])/apply(table(seuratData@meta.data[c('stage',"G0")]),1,sum)
df=data.frame(t)
p5=ggplot(df,aes(stage,Freq))+geom_line(aes(group=G0,col=G0))+geom_text(aes(label=round(Freq,3)))
pdf(paste0("20210831_G0_corthreshold",cor_threshold,"gene_threshold",gene_threshold,".pdf"))
print(p1)
print(p2.1)
print(p2)
print(p3)
print(p4)
print(p5)
b=barplot(a)
text(b,a/2,labels = round(a,2))
dev.off()
