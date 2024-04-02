### perform gene expression and pathway activity difference analysis based on TCGA-BRCA data

library(gridExtra)
library(lattice)
library(showtext)
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

load("E:/databases/pathway_genesets/02_kegg_API_download/kegg_api_geneSets20221019.RData")
load("E:/micedata/chipseq/20210113_all_m4_mut_chipseq_rabit_bcl11b/20210113_WT_VS_mut_bcl11b_pos_neg_targetmarkergene.RData")
load("E:/micedata/chipseq/bcl11b_chipseq_targetgene.RData")

names(hsa_kegg_list_des)=gsub(" - Homo sapiens \\(human\\)","",names(hsa_kegg_list_des))
forfigure=xlsx::read.xlsx("X:/analyze_data/FEdata/baihuiru/20211225_rmm4v1_m24v2/1virgin/20211228for figure.xlsx",sheetIndex = 2)
gs1=hsa_kegg_list_des[match(forfigure$Description,names(hsa_kegg_list_des),nomatch = 0)]
senesc_pathways <- fgsea::gmtPathways("X:/analyze_data/FEdata/publicdata/agingdata/code/senesccence_gene_sets.gmt")
bcl11btarget=toupper(bcl11b_chipseq_targetgene$SYMBOL)
markergenes=markergenes[markergenes$p_val<0.05,]
Bcl11bPosGene=toupper(intersect(bcl11b_chipseq_targetgene$SYMBOL,rownames(markergenes)[markergenes$avg_logFC<0]))
Bcl11bNegGene=toupper(intersect(bcl11b_chipseq_targetgene$SYMBOL,rownames(markergenes)[markergenes$avg_logFC>0]))
gs=c(gs1,senesc_pathways)
gs$bcl11btarget=bcl11btarget
gs$bcl11btargetpos=Bcl11bPosGene
gs$bcl11btargetneg=Bcl11bNegGene

gs=gs[c("bcl11btarget","bcl11btargetpos","bcl11btargetneg","Cellular senescence","COURTOIS_SENESCENCE_TRIGGERS","REACTOME_FORMATION_OF_SENESCENCE_ASSOCIATED_HETEROCHROMATIN_FOCI_SAHF","REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE","ECM-receptor interaction","Notch signaling pathway","NF-kappa B signaling pathway")]

pancancer_RNA=data.table::fread("X:/analyze_data/FEdata/publicdata/agingdata/code/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")%>%
  filter(!grepl("^\\?",gene_id))%>%
  mutate(gene_id=gsub("[|].*","",gene_id))%>%
  filter(!duplicated(gene_id))%>%
  column_to_rownames("gene_id")
pancancer_RNA0=pancancer_RNA
pancancer_RNA=pancancer_RNA0[,substr(colnames(pancancer_RNA0),14,15)%in%c("01","06","11")]

# a=colnames(pancancer_RNA)
# type=substr(a,14,15)
# t=pancancer_RNA[,type=="11"]
survival_data=openxlsx::read.xlsx("X:/analyze_data/FEdata/publicdata/agingdata/code/survival_data_all.xlsx")
pancancer_RNA=pancancer_RNA[,substr(colnames(pancancer_RNA),1,12)%in%survival_data$bcr_patient_barcode]


pdf(paste0(dataname0,"-TCGA-pathway",".pdf"),width = 10,height = 4)
plist=list()
for(type in sort(unique(survival_data$type))){
  print(type)

  brcaclinical=survival_data[survival_data$type==type,]
  BRCA=pancancer_RNA[,substr(colnames(pancancer_RNA),1,12)%in%brcaclinical$bcr_patient_barcode]
  patient=data.frame(table(substr(colnames(BRCA),1,12)))
  samples=sort(colnames(BRCA)[substr(colnames(BRCA),1,12)%in%patient$Var1[patient$Freq>1]])
  samples=samples[substr(samples,14,16)%in%c("01A","06A","11A","11B")]
  
  
  p=samples[substr(samples,14,16)%in%c("01A")]
  p=p[!duplicated(substr(p,1,12))]
  n=samples[substr(samples,14,16)%in%c("11A","11B")]
  n=n[!duplicated(substr(n,1,12))]
  # m=samples[substr(samples,14,16)%in%c("06A")]
  # m=m[!duplicated(substr(m,1,12))]
  
  p=p[substr(p,1,12)%in%substr(n,1,12)]
  n=n[substr(n,1,12)%in%substr(p,1,12)]
  if(length(n)<5){next()}
  # m=m[substr(m,1,12)%in%substr(p,1,12)]#3
  samples=c(n,p)
  # pam50=brcaclinical$PAM50[match(substr(samples,1,12),brcaclinical$patient)]
  # samples=samples[!is.na(pam50)]
  group=as.numeric(substr(samples,14,15))
  
  group=plyr::mapvalues(group,from = c(1,11),to=c("primary","normal"))
  
  # enrich.mat <- GSVA::gsva(expr = as.matrix(BRCA[,samples]), gset.idx.list = gs, kcdf = "Gaussian",method = "ssgsea", verbose = TRUE)
  BRCA[is.na(BRCA)]=0
  enrich.mat=list()
  for(i in names(gs)){
    g=gs[[i]]
    g=g[g%in%rownames(BRCA)]
    t=apply(BRCA[g,],2,mean)
    enrich.mat[[i]]=t
  }
  enrich.mat=t(do.call("cbind",enrich.mat))
  
  enrich.mat=enrich.mat[,samples]
  # enrich.mat[is.na(enrich.mat)]=0
  enrich.mat=t(scale(t(enrich.mat)))
  
  
  df=enrich.mat%>%data.frame(.)%>%tibble::rownames_to_column("pathway")%>%data.table::melt("pathway")%>%
    dplyr::rename(TCGA_barcode=variable,score=value)%>%
    dplyr::mutate(TCGA_barcode = gsub("\\.","-", TCGA_barcode),sampletype=group[match(TCGA_barcode,samples)])
  df$sampletype=factor(df$sampletype,levels = c("normal","primary","metastasis"))
  df$sampletype=droplevels(df$sampletype)
  # df$PAM50=brcaclinical$PAM50[match(substr(df$TCGA_barcode,1,12),brcaclinical$patient)]
  df$patient=substr(df$TCGA_barcode,1,12)
  df$location=plyr::mapvalues(substr(df$TCGA_barcode,14,15),from = c("01","11"),to=c("primary","normal"))
  sum(is.na(df))
  df=df[df$location%in%c("normal","primary"),]
  df$pathway=gsub("COURTOIS_SENESCENCE_TRIGGERS","Courtois senescence triggers",df$pathway)
  df$pathway=gsub("REACTOME_FORMATION_OF_SENESCENCE_ASSOCIATED_HETEROCHROMATIN_FOCI_SAHF","Formation of senescence associated heterochromatin foci",df$pathway)
  df$pathway=gsub("REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE","DNA damage telomere stress induced senescence",df$pathway)
  df$pathway=gsub("COURTOIS_SENESCENCE_TRIGGERS","Courtois senescence triggers",df$pathway)
  df$score[df$pathway=="bcl11btargetneg"]=-1*df$score[df$pathway=="bcl11btargetneg"]
  df$pathway=factor(df$pathway,levels = df$pathway[!duplicated(df$pathway)])
  
  # for(pathway in sort(unique(df$pathway))){
  #   print(paste(type,pathway))
  #   comp=c("normal","primary","metastasis")
  #   tmp=df[df$pathway%in%pathway,]
  #   p1=tmp%>%ggpaired(x="location",y="score",id="patient", fill = 'location',line.color = scales::alpha("gray", 2/3))+
  #     stat_compare_means(method = 't.test', hide.ns = T,paired = T)+labs(title=pathway)
  #   # p1.1=tmp%>%ggpaired(x="location",y="score",id="patient", fill = 'location',line.color = scales::alpha("gray", 2/3),facet.by = "PAM50")+
  #   #   stat_compare_means(method = 't.test', hide.ns = T,paired = T)+labs(title=pathway)
  #   
  #   print(p1)
  #   # print(p1.1)
  #   
  # }
  
  p=ggplot(df,aes(location,score,fill=location))+geom_boxplot(outlier.colour = NA)+facet_wrap(~pathway,ncol = 10,scales = "free_y")+
    ggpubr::stat_compare_means(method = "t.test",hide.ns = T,paired = T,comparisons = list(c("normal","primary")))+
    theme_classic()+scale_fill_brewer(palette = "Set1")+labs(y=paste("score in ",type))+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
  plist[[type]]=p
  print(p)
  
}
dev.off()
pdf(paste0(dataname0,"-TCGA-pathway",".pdf"),width = 10,height =50)
grid.arrange(grobs=plist,ncol=1)
dev.off()