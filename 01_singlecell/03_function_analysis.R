### perform function analysis
# enrichment analysis
# Bcl11b activity analysis
# pathway activity analysis
# transcription factor enrichment analysis along pseudotime


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
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
load("20220428_msigdblist_human&mice.RData")#msigdb
load("mouse.gene2id.RData")

#------------------clusterProfiler enrichment------------------

enrichfunciton=list()
gsea_list=list()
for(celltype in names(markerlist)){
  # print(celltype)
  
  # print(celltype2)
  
  print(paste(celltype,celltype2))
  m=markerlist[[celltype]][[celltype2]]
  genes=m[order(m$avg_log2FC,decreasing = T),]
  geneList=genes$avg_log2FC
  names(geneList)=rownames(genes)
  gsea= GSEA(geneList, TERM2GENE = pathwaylist)
  gsea_list[[celltype]][[celltype2]]=gsea
  
  m=m[m$p_val<0.05,]  ###no pval filtering
  for(j in c("high","low")){
    print(j)
    genes=rownames(m)[m$avg_log2FC>0]
    if(j=="low"){
      genes=rownames(m)[m$avg_log2FC<0]
    }
    
    geneID=gene2id$ENTREZID[match(genes,gene2id$SYMBOL,nomatch=0)]
    ego=enrichGO(gene = geneID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    ekegg <-enrichKEGG(gene = geneID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)
    
    if(!is.null(ego)){
      ego=ego@result[ego@result$pvalue<0.05,]
      enrichfunciton[[celltype]][[celltype2]][[j]][["ego"]]=ego
    }
    if(!is.null(ekegg)){
      ekegg=ekegg@result[ekegg@result$pvalue<0.05,]
      ekegg$genesymbol=sapply(ekegg$geneID,function(x){a=unlist(strsplit(x,"/"));b=gene2id$SYMBOL[match(a,gene2id$ENTREZID,nomatch = 0)];return(paste0(b,collapse = "/"))})
      enrichfunciton[[celltype]][[celltype2]][[j]][["ekegg"]]=ekegg

    }
  }
  
  
  
}
function_RNA=list(fun=enrichfunciton,gene=markerlist,gsea=gsea_list)
save(function_RNA,file = "function_RNA.RData")

#------------------pathway activity------------------
load("E:/databases/pathway_genesets/02_kegg_API_download/kegg_api_geneSets20221019.RData")

keggmmu=mmu_kegg_list_des
keggmmu=keggmmu[1:(length(keggmmu)-1)]
names(keggmmu)=gsub(" - Mus musculus \\(mouse\\)","_",names(keggmmu))
names(keggmmu)=gsub(" ","_",names(keggmmu))
names(keggmmu)=gsub("_$","",names(keggmmu))
names(keggmmu)=gsub("[()/,-]","",names(keggmmu))
names(keggmmu)=gsub("____|___|__","_",names(keggmmu))

s=seuratData
s=AddModuleScore(s,features = c(keggmmu),seed = 100)
colnames(s@meta.data)[match("Cluster1",colnames(s@meta.data)):ncol(s@meta.data)]=c(names(keggmmu))
pdf("compare_pathway_activity.pdf")
for(i in sort(names(keggmmu))){
  # p1=VlnPlot(s,i,group.by = "stage")
  
  data=s@meta.data
  colnames(data)[match(i,colnames(data))]="exp"
  p2=ggplot(data,aes(stage,exp))+geom_violin()+labs(title = i)+geom_jitter(aes(col=stage))+geom_boxplot(width=.1,outlier.shape = NA)+
    ggpubr::stat_compare_means(comparisons = list(c("s1","s2"),c("s2","s3"),c("s3","s4")))+
    # geom_hline(yintercept = 0)+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),panel.border = element_blank(),panel.grid = element_blank(),axis.line = element_line(size=1,colour = 'black'))+
    scale_colour_brewer(palette = "Set1")
  p3=ggplot(data,aes(Pseudotime,exp))+geom_point()+
    ggpubr::stat_cor(label.x.npc = "middle")+
    stat_smooth(se=T,size=1.2)+#method="lm"
    ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x.npc = "middle")+
    labs(title=i)+
    # scale_color_gradient(low = "grey",high = "red")+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),
                     panel.border = element_blank(),panel.grid = element_blank(),
                     axis.line = element_line(size=0.5,colour = 'black'))
  # print(p1)
  print(p2)
  print(p3)
}
dev.off()

###
#------------------bcl11b activity------------------
##reanalyze bcl11b activity along with pseudotime

seuratData=seuratData_virgin

load("20210113_WT_VS_mut_bcl11b_pos_neg_targetmarkergene.RData")
load("bcl11b_chipseq_targetgene.RData")
markergenes=markergenes[markergenes$p_val<0.05,]
Bcl11bPosGene=intersect(bcl11b_chipseq_targetgene$SYMBOL,rownames(markergenes)[markergenes$avg_logFC<0])
Bcl11bNegGene=intersect(bcl11b_chipseq_targetgene$SYMBOL,rownames(markergenes)[markergenes$avg_logFC>0])

bcl11beffect=list()
for(type in c("Bcl11bPosGenep0.05","Bcl11bNegGenep0.05")){
  
  
  print(type)
  bcl11btarget=Bcl11bPosGene
  if(type=="Bcl11bNegGenep0.05"){
    bcl11btarget=Bcl11bNegGene
  }
  cellGeneExp=c()
  
  for(i in rownames(pData(HSMM))){
    t1=match(bcl11btarget,rownames(seuratData@assays$RNA@data),nomatch = 0)
    bcl11btarget=bcl11btarget[t1!=0]
    targetAverageExp=sum(seuratData@assays$RNA@data[t1,match(i,colnames(seuratData))]*exp(markergenes[bcl11btarget,"avg_logFC"]))/sum(exp(markergenes[bcl11btarget,"avg_logFC"]))
    cellGeneExp=c(cellGeneExp,targetAverageExp)
  }
  bcl11beffect[[type]]=cellGeneExp

}


data=do.call("cbind",bcl11beffect)
colnames(data)=c("Pos","Neg")
pData(HSMM)$rank=rank(pData(HSMM)$Pseudotime)
data=cbind(data,pData(HSMM))
data=data[,c("Pos","Neg","Pseudotime","rank")]
data$name=rownames(data)
data2=data.frame(Pseudotime=c(data$Pseudotime,data$Pseudotime),name=c(data$name,data$name),bcl11beffect=c(data$Pos,data$Neg),type=rep(c("Pos","Neg"),each=nrow(data)),rank=c(data$rank,data$rank),cell=c(rownames(data),rownames(data)))

data2$bcl11beffect_1=-1*data2$bcl11beffect
data2$stage=seuratData$stage[match(data2$name,colnames(seuratData))]


p1=ggplot(data2,aes(x=Pseudotime,y=bcl11beffect,col=type))+geom_point(size=0.1,alpha=0.3)+
  stat_smooth(se=T,size=1.2)+
  facet_wrap(~type)+
  ggpubr::stat_cor(label.x.npc = "middle")+
    ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x.npc = "middle")+
  #+geom_line(stat = 'density',aes(col=month))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),panel.border = element_blank(),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = 'black'))+
  scale_colour_brewer(palette = "Set1")


p3=ggplot(data2,aes(x=Pseudotime,y=bcl11beffect_1,col=type))+geom_point(size=0.1,alpha=0.3)+
  stat_smooth(se=T,size=1.2)+
  facet_wrap(~type)+
  ggpubr::stat_cor(label.x.npc = "middle")+
    ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x.npc = "middle")+
  #+geom_line(stat = 'density',aes(col=month))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),panel.border = element_blank(),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = 'black'))+
  scale_colour_brewer(palette = "Set1")


p4.1=ggplot(data2[data2$type=="Neg",],aes(Pseudotime,bcl11beffect_1,col=type))+geom_point(size=0.1,alpha=0.3)+
  facet_wrap(~type,scales = "free")+
  ggpubr::stat_cor(label.x.npc = "middle",label.y.npc = "middle")+
  stat_smooth(se=T,size=1.2)+#coord_cartesian(ylim=c(0,3))+
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x=0,label.y=-0.8)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),axis.text.x = element_text(angle = 0),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = "black"))+
  scale_colour_brewer(palette = "Set1")+coord_cartesian(ylim = c(-0.81,-0.5))
p4.3=ggplot(data2[data2$type=="Neg",],aes(rank,bcl11beffect_1,col=type))+geom_point(size=0.1,alpha=0.3)+
  facet_wrap(~type,scales = "free")+
  ggpubr::stat_cor(label.x.npc = "middle")+#coord_cartesian(ylim=c(0,3))+
  stat_smooth(se=T,size=1.2)+
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x.npc = "middle")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),axis.text.x = element_text(angle = 0),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = "black"))+
  scale_colour_brewer(palette = "Set1")

pdf(paste0(dataname0,"bcl11beffect.pdf"))
print(p1)
print(p3)
print(p4.1)
print(p4.3)
dev.off()

pdf("bcl11b_neg.pdf")
p4.1
ggplot(data2[data2$type=="Neg",],aes(Pseudotime,bcl11beffect_1,col=type))+
  # geom_point(size=0.1,alpha=0.3)+
  facet_wrap(~type,scales = "free")+
  ggpubr::stat_cor(label.x.npc = "middle",label.y = -0.7)+
  stat_smooth(se=T,size=1.2)+
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x.npc = "middle",label.y = -0.8)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),axis.text.x = element_text(angle = 0),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = "black"))+
  scale_colour_brewer(palette = "Set1")+coord_cartesian(ylim = c(-0.81,-0.5))
ggplot(data2[data2$type=="Neg",],aes(Pseudotime,bcl11beffect_1,col=type))+
  # geom_point(size=0.1,alpha=0.3)+
  facet_wrap(~type,scales = "free")+
  ggpubr::stat_cor(label.x.npc = "middle")+
  stat_smooth(se=F,size=1.2)+
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x.npc = "middle")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),axis.text.x = element_text(angle = 0),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = "black"))+
  scale_colour_brewer(palette = "Set1")+coord_cartesian(ylim = c(-0.81,-0.5))
ggplot(data2[data2$type=="Neg",],aes(Pseudotime,bcl11beffect,col=type))+
  # geom_point(size=0.1,alpha=0.3)+
  facet_wrap(~type,scales = "free_y")+
  ggpubr::stat_cor(label.x = 0,label.y.npc = "middle")+#ylim(c(0,3))+
  stat_smooth(se=F,size=1.2)+
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x.npc = "middle",label.y.npc = "middle")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),axis.text.x = element_text(angle = 0),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = "black"))+
  scale_colour_brewer(palette = "Set1")+coord_cartesian(ylim = c(0.5,0.81))
ggplot(data2[data2$type=="Pos",],aes(Pseudotime,bcl11beffect,col=type))+
  facet_wrap(~type,scales = "free_y")+
  ggpubr::stat_cor(label.x = 0,label.y.npc = "middle")+#ylim(c(0,3))+
  stat_smooth(se=F,size=1.2)+
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.x.npc = "middle",label.y.npc = "middle")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),axis.text.x = element_text(angle = 0),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = "black"))+
  scale_colour_brewer(palette = "Set1")+coord_cartesian(ylim = c(0.775,0.79))
dev.off()

#------------------transcription factor enrichment along pseudotime------------------

hyperTest <- function(gene1, gene2, backgroundGenes=NULL, type="over"){
  q <- length(intersect(gene1, gene2))
  m <- length(gene2);
  background <- length(backgroundGenes)
  n <- background-m;
  k <- length(gene1);
  if(type == "over"){
    return(phyper(q-1, m, n, k, FALSE));	#return p value
  }
}

#TFs regulate targets along pseudotime
stage_genes=xlsx::read.xlsx("pseudotime_gene_function.xlsx",sheetName = "gene_cluster")

stage_gene_list=lapply(paste0("cluster",1:4),function(x){return(stage_genes$cell[stage_genes$cellcluster2==x])})
names(stage_gene_list)=paste0("s",1:4)

TF_target_list_to_mouse=readRDS("E:/databases/TF_target/20210615_TF_target_list_mouse.rds")
TF_target_list_to_mouse_stage=list()
for(gs in names(TF_target_list_to_mouse)){
  print(gs)
  tmpdb=TF_target_list_to_mouse[[gs]]
  backgroundgenes=unique(c(do.call("c",stage_gene_list),do.call("c",tmpdb)))
  
  for(stage in names(stage_gene_list)){
    genes=stage_gene_list[[stage]]
    a=sapply(tmpdb,function(x){hyperTest(gene1=genes,gene2 = x,backgroundGenes = backgroundgenes)})
    TF_target_list_to_mouse_stage[[gs]][[stage]]=names(a)[a<0.2]

  }
}


for(stage in paste0("s",1:4)){
  print(stage)
  allTF=unique(do.call("c",lapply(TF_target_list_to_mouse_stage,function(x){return(x[[stage]])})))
  m=matrix(0,nrow = length(allTF),ncol = length(TF_target_list_to_mouse_stage))
  rownames(m)=allTF
  colnames(m)=names(TF_target_list_to_mouse_stage)
  for(gs in names(TF_target_list_to_mouse_stage)){
    m[TF_target_list_to_mouse_stage[[gs]][[stage]],gs]=1
  }
  m=cbind(m,sum=apply(m,1,sum))
  m=m[order(m[,"sum"],decreasing = T),]
  xlsx::write.xlsx(m,file = paste0(dataname0,"stage_TF_intersection.xlsx"),sheetName = stage,append = T)
}

