###create Seurat object and filter some contaminated cells

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

gzfile=list.files(path = 'E:/singlecell/20201104_mito_mutation/allcount/',pattern = '.tsv.gz')

countdata=list()
for(i in gzfile){
  print(i)
  countdata[[gsub('.count.*','',i)]]=read.table(paste0('E:/singlecell/20201104_mito_mutation/with_mito_allcount/',i),header = T)
  
}
names(countdata)=gsub("_pri","",names(countdata)) 


barcode=read.table('E:/singlecell/all_count/Barcode.txt',header = F)
allbarcode=barcode$V2

for(i in names(countdata)){
  countdata[[i]]=countdata[[i]][,c('gene',intersect(allbarcode,colnames(countdata[[i]])))]
}

for(i in names(countdata)){
  colnames(countdata[[i]])=c('gene',paste(i,colnames(countdata[[i]])[-1],sep = '_'))
}
lapply(countdata,dim)

data=countdata[[1]]
for(i in 2:length(gzfile)){
  
  data=merge(data,countdata[[i]],by='gene',all=T)
  
}

rownames(data)=data$gene
data=data[,-1]
data[is.na(data)]=0
data=as(as.matrix(data),'sparseMatrix')


dataname0="20211223_virgin"
seuratData <- CreateSeuratObject(data) 
seuratData@meta.data$month_version=sapply(rownames(seuratData@meta.data),function(x){paste0(unlist(strsplit(x,'_'))[1:2],collapse = '_')})
seuratData@meta.data$month_version=factor(seuratData@meta.data$month_version,
                                          levels = c(paste0(rep("m",30*21),paste0(rep(1:30,each=21),paste0("_v",c("",1:20)))),
                                          paste0(rep("m",30*21),paste0(rep(1:30,each=21),paste0("_p",c("",1:20)))),
                                          paste0("mut_v",1:4),
                                          c("DMBA_v1","DMBA_v2","TmutLum_v","TmutBasal_v","TwtLum_v","TwtBasal_v","17_s")
                                          ))

seuratData@meta.data$month_version=droplevels(seuratData@meta.data$month_version)
seuratData@meta.data$month=sapply(colnames(seuratData),function(x){unlist(strsplit(x,'_'))[1]})
seuratData@meta.data$month=factor(seuratData@meta.data$month,levels =c( paste0('m',c(1:30)),'mut',"DMBA","TmutLum","TmutBasal","TwtLum","TwtBasal"))
seuratData@meta.data$month=droplevels(seuratData@meta.data$month)
table(seuratData@meta.data$month)
seuratData@meta.data$virgin_parous="virgin"
seuratData@meta.data$virgin_parous[grep("p",seuratData@meta.data$month_version)]="parous"
seuratData@meta.data$virgin_parous=as.factor(seuratData@meta.data$virgin_parous)
seuratData@meta.data$virgin_parous=droplevels(seuratData@meta.data$virgin_parous)

seuratData@meta.data$batcheffect=sapply(seuratData@meta.data$month_version,function(x){if(x%in%c("m2_v1","m14_p","m19_p")){return('20191105&20191022')}else if(x%in%c("m4_v1","m4_v1_pri","mut_v1","mut_v1_pri")){return('20190728_merge')}else if(x%in%c('m4_v2',"mut_v2")){return('20190813_merge')}else if(x%in%c("m4_v3","m12_p","m17_p")){return('20190905&20190823')}else if(x%in%c("m2_v2","m4_v4","m9_v1","m9_v2","m11_v1","m11_v2","m15_v1","m15_v2","m30_p1","m30_p2","mut_v3","mut_v4")){return("20191221")}else if(x%in%c("m4_p1","m4_p2","m7_v1","m7_v1_pri","m7_v2_pri","m7_v2","DMBA_v1","DMBA_v2")){return("20200326")}else if(x%in%c("m13_v1","m13_v2","m15_v3","m15_v4")){return("20200421")}else if(x%in%c("m13_v3")){return("20200530_new_m13")}else if(x%in%"m13_v13"){return("20200616_m13_merged")}else if(x%in%c("m17_v1","m17_v2")){return("20200708")}else if(x%in%c("m19_v1","m19_v2")){return("20200830")}else if(x%in%c("TmutLum_v","TmutBasal_v","TwtLum_v","TwtBasal_v","m17_s")){return("20200801_merge")}else if(x%in%c("m22_v1","m22_v2","m24_v1")){return("20201228")}else if(x%in%c("m24_v2",paste0("m29_v",1:3))){return("20211222")}else{return(0)}})
seuratData@meta.data$batcheffect=as.factor(seuratData@meta.data$batcheffect)
seuratData@meta.data$batcheffect=droplevels(seuratData@meta.data$batcheffect)


seuratData[["percent.mt"]] = colSums(as.matrix(seuratData@assays$RNA@counts[grep('^mt-',rownames(seuratData@assays$RNA@counts)),]))/seuratData$nCount_RNA
seuratData[["percent.ercc"]] = colSums(as.matrix(seuratData@assays$RNA@counts[grep('^ERCC-',rownames(seuratData@assays$RNA@counts)),]))/seuratData$nCount_RNA
p1=VlnPlot(seuratData, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 2, group.by = "month")
p2=VlnPlot(seuratData, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 2, group.by = "month_version")
p3=VlnPlot(seuratData, c('percent.mt','percent.ercc'), pt.size = 0.1, ncol = 2, group.by = "month")
p4=VlnPlot(seuratData, c('percent.mt','percent.ercc'), pt.size = 0.1, ncol = 2, group.by = "month_version")

seuratData <- subset(seuratData, subset = nFeature_RNA > 200&nFeature_RNA < 6000 & percent.mt < 0.05 & percent.ercc<0.1)
p5=VlnPlot(seuratData, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 2, group.by = "month")
p6=VlnPlot(seuratData, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 2, group.by = "month_version")
p7=VlnPlot(seuratData, c('percent.mt','percent.ercc'), pt.size = 0.1, ncol = 2, group.by = "month")
p8=VlnPlot(seuratData, c('percent.mt','percent.ercc'), pt.size = 0.1, ncol = 2, group.by = "month_version")
pdf(paste0(dataname0,"gene_count.pdf"))
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
dev.off()
p1=ggplot(seuratData@meta.data,aes(x=month,y=nCount_RNA,col=virgin_parous))+geom_boxplot()+geom_jitter(size=0.1)+labs(title = "nCount_RNA")+
  theme_bw()+theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),plot.title = element_text(hjust=0.5),legend.position = "none")
p1.1=ggplot(seuratData@meta.data,aes(x=month,y=nFeature_RNA,col=virgin_parous))+geom_boxplot()+geom_jitter(size=0.1)+labs(title = "nFeature_RNA")+
  theme_bw()+theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),plot.title = element_text(hjust=0.5),legend.position = "none")

p2=ggplot(seuratData@meta.data,aes(x=month_version,y=nCount_RNA,col=virgin_parous))+geom_boxplot()+geom_jitter(size=0.1)+labs(title = "nCount_RNA")+
  theme_bw()+theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),plot.title = element_text(hjust=0.5),legend.position = "none")
p2.1=ggplot(seuratData@meta.data,aes(x=month_version,y=nFeature_RNA,col=virgin_parous))+geom_boxplot()+geom_jitter(size=0.1)+labs(title = "nFeature_RNA")+
  theme_bw()+theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),plot.title = element_text(hjust=0.5),legend.position = "none")

p3=ggplot(seuratData@meta.data,aes(x=month,y=percent.mt,col=virgin_parous))+geom_boxplot()+geom_jitter(size=0.1)+labs(title = "percent.mt")+
  theme_bw()+theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),plot.title = element_text(hjust=0.5),legend.position = "none")
p3.1=ggplot(seuratData@meta.data,aes(x=month,y=percent.ercc,col=virgin_parous))+geom_boxplot()+geom_jitter(size=0.1)+labs(title = "percent.ercc")+
  theme_bw()+theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),plot.title = element_text(hjust=0.5),legend.position = "none")

p4=ggplot(seuratData@meta.data,aes(x=month_version,y=percent.mt,col=virgin_parous))+geom_boxplot()+geom_jitter(size=0.1)+labs(title = "percent.mt")+
  theme_bw()+theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),plot.title = element_text(hjust=0.5),legend.position = "none")
p4.1=ggplot(seuratData@meta.data,aes(x=month_version,y=percent.ercc,col=virgin_parous))+geom_boxplot()+geom_jitter(size=0.1)+labs(title = "percent.ercc")+
  theme_bw()+theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),plot.title = element_text(hjust=0.5),legend.position = "none")
pdf(paste0(dataname0,"gene_count2.pdf"))
print(p1+p1.1)
print(p2+p2.1)
print(p3+p3.1)
print(p4+p4.1)
dev.off()
getwd()

seuratData=seuratData[grep("ERCC",rownames(seuratData),invert = T),]
seuratData=NormalizeData(object = seuratData, normalization.method = "LogNormalize", scale.factor = 10000)
seuratData <- CellCycleScoring(seuratData, s.features = s.genes$MGI.symbol, g2m.features =g2m.genes$MGI.symbol)#, set.ident = TRUE)

seuratData <- FindVariableFeatures(seuratData, selection.method = "vst", x.cutoff=c(0.01,Inf),y.cutoff=0.01,nfeatures = 5000)
seuratData <- ScaleData(object = seuratData, vars.to.regress = c("S.Score", "G2M.Score",'batcheffect','percent.mt','percent.ercc','nFeature_RNA'))

seuratData <- RunPCA(seuratData, npcs = 70, verbose = FALSE,features =VariableFeatures(seuratData) )

dims=10
seuratData <- RunTSNE(seuratData, dims = 1:dims, verbose = FALSE,check_duplicates = FALSE)
seuratData <- FindNeighbors(seuratData,reduction = 'pca', verbose = FALSE, dims = 1:dims)
seuratData <- FindClusters(seuratData, algorithm = 1, random.seed = 256, resolution = 1)
pdf(paste0(dataname0,"cluster.pdf"))
DimPlot(seuratData,  reduction = "tsne", group.by = "seurat_clusters", label = TRUE)+scale_fill_manual(values = getPalette(length(unique(seuratData@meta.data$seurat_clusters))))
DimPlot(seuratData,  reduction = "tsne", group.by = "month", label = TRUE)+scale_fill_manual(values = getPalette(length(unique(seuratData@meta.data$month))))
DimPlot(seuratData,  reduction = "tsne", group.by = "month_version", label = TRUE)+scale_fill_manual(values = getPalette(length(unique(seuratData@meta.data$month_version))))
DimPlot(seuratData,  reduction = "tsne", group.by = "batcheffect", label = TRUE)+scale_fill_manual(values = getPalette(length(unique(seuratData@meta.data$batcheffect))))
dev.off()


##rm endothelium, blood cell, stromal

dataname0=paste0(dataname0,"rm_stromal")

keep=apply(seuratData@assays$RNA@data[c("Pecam1","Ptprc","Lyve1","Col1a1"),],1,function(x){x<0.1})
keep=apply(keep,1,sum)
seuratData=seuratData[,keep==4]
seuratData <- RunPCA(seuratData, npcs = 70, verbose = FALSE,features =VariableFeatures(seuratData) )

dims=10
seuratData <- RunTSNE(seuratData, dims = 1:dims, verbose = FALSE,check_duplicates = FALSE)
seuratData <- FindNeighbors(seuratData,reduction = 'pca', verbose = FALSE, dims = 1:dims)

seuratData <- FindClusters(seuratData, algorithm = 1, random.seed = 256, resolution = 1)
pdf(paste0(dataname0,"cluster_rmstromal.pdf"))
DimPlot(seuratData,  reduction = "tsne", group.by = "seurat_clusters", label = TRUE)+scale_fill_manual(values = getPalette(length(unique(seuratData@meta.data$seurat_clusters))))
DimPlot(seuratData,  reduction = "tsne", group.by = "month", label = TRUE)+scale_fill_manual(values = getPalette(length(unique(seuratData@meta.data$month))))
DimPlot(seuratData,  reduction = "tsne", group.by = "month_version", label = TRUE)+scale_fill_manual(values = getPalette(length(unique(seuratData@meta.data$month_version))))
DimPlot(seuratData,  reduction = "tsne", group.by = "batcheffect", label = TRUE)+scale_fill_manual(values = getPalette(length(unique(seuratData@meta.data$batcheffect))))
dev.off()
save(seuratData,file = paste0(dataname0,"seuratDatarm1col1a1.2.RData"))
