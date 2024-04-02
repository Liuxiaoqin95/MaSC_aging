### pseudotime analysis based on filtered seuratData object

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

dataname0="20210831"
feature_sheet=data.frame(row.names = rownames(seuratData@assays$RNA@counts),gene_short_name=rownames(seuratData@assays$RNA@counts))
pd <- new("AnnotatedDataFrame", data = seuratData@meta.data)
fd <- new("AnnotatedDataFrame", data =feature_sheet)

HSMM=newCellDataSet(seuratData@assays$RNA@counts,
                    phenoData = pd,
                    featureData = fd,
                    lowerDetectionLimit=0.5,
                    expressionFamily=negbinomial.size())

HSMM <- detectGenes(HSMM, min_expr = 0.1)
HSMM <- estimateSizeFactors(HSMM)
HSMM=estimateDispersions(HSMM)
disp_table <- dispersionTable(HSMM)
HSMM_ordering_genes <- subset(disp_table, mean_expression >= 0.012 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <-setOrderingFilter(HSMM,ordering_genes = HSMM_ordering_genes)
HSMM <-reduceDimension(HSMM,max_components = 2, reduction_method = 'DDRTree',residualModelFormulaStr = "~Size_Factor+num_genes_expressed+percent.ercc+percent.mt+S.Score+G2M.Score+batcheffect")#

HSMM <-orderCells(HSMM)


p1=plot_cell_trajectory(HSMM, color_by = "State",cell_size = 1)+
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$State))), name = "State")

p2=plot_cell_trajectory(HSMM, color_by = "seurat_clusters",cell_size = 1)+ 
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$seurat_clusters))), name = "seurat_clusters")
p3=plot_cell_trajectory(HSMM, color_by="month",cell_size = 1,show_tree = T)+labs(title =which.max(table(pData(HSMM)[,c('month','State')])[1,]))+
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$month))), name = "month")
p4=plot_cell_trajectory(HSMM, color_by="month_version",cell_size = 1,show_tree = T)+ 
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$month_version))), name = "month_version")
p1.1=plot_cell_trajectory(HSMM, color_by = "State",cell_size = 1)+
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$State))), name = "State")+NoLegend()

p2.1=plot_cell_trajectory(HSMM, color_by = "seurat_clusters",cell_size = 1)+ 
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$seurat_clusters))), name = "seurat_clusters")+NoLegend()
p3.1=plot_cell_trajectory(HSMM, color_by="month",cell_size = 1,show_tree = T)+labs(title =which.max(table(pData(HSMM)[,c('month','State')])[1,]))+
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$month))), name = "month")+NoLegend()
p4.1=plot_cell_trajectory(HSMM, color_by="month_version",cell_size = 1,show_tree = T)+ 
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$month_version))), name = "month_version")+NoLegend()

pdf(paste0(dataname0,"mean_exp",mean_exp,"trajectory.pdf"))
print(p1)
print(p2)
print(p3)
print(p4)
print(p1.1)
print(p2.1)
print(p3.1)
print(p4.1)
plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 1,theta = 0)
plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 1,theta = 0)+NoLegend()

dev.off()


p1=plot_cell_trajectory(HSMM, color_by = "State",cell_size = 1,theta = 180)+
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$State))), name = "State")
p2=plot_cell_trajectory(HSMM, color_by = "seurat_clusters",cell_size = 1,theta = 180)+ 
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$seurat_clusters))), name = "seurat_clusters")
p3=plot_cell_trajectory(HSMM, color_by="month",cell_size = 1,show_tree = T,theta = 180)+labs(title =which.max(table(pData(HSMM)[,c('month','State')])[1,]))+
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$month))), name = "month")
p4=plot_cell_trajectory(HSMM, color_by="month_version",cell_size = 1,show_tree = T,theta = 180)+ 
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$month_version))), name = "month_version")
p1.1=plot_cell_trajectory(HSMM, color_by = "State",cell_size = 1,theta = 180)+
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$State))), name = "State")+NoLegend()
p2.1=plot_cell_trajectory(HSMM, color_by = "seurat_clusters",cell_size = 1,theta = 180)+ 
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$seurat_clusters))), name = "seurat_clusters")+NoLegend()
p3.1=plot_cell_trajectory(HSMM, color_by="month",cell_size = 1,show_tree = T,theta = 180)+labs(title =which.max(table(pData(HSMM)[,c('month','State')])[1,]))+
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$month))), name = "month")+NoLegend()
p4.1=plot_cell_trajectory(HSMM, color_by="month_version",cell_size = 1,show_tree = T,theta = 180)+ 
  scale_color_manual(values =rainbow(length(unique(pData(HSMM)$month_version))), name = "month_version")+NoLegend()

pdf(paste0(dataname0,"mean_exp",mean_exp,"trajectory_theta180.pdf"))
print(p1)
print(p2)
print(p3)
print(p4)
print(p1.1)
print(p2.1)
print(p3.1)
print(p4.1)
dev.off()

HSMM <-orderCells(HSMM,root_state = which.max(table(pData(HSMM)[,c('month','State')])[1,]))
p1=ggplot(pData(HSMM),aes(x=Pseudotime,col=Phase))+geom_density(size=1)+scale_color_manual(values = rainbow(length(unique(pData(HSMM)$Phase))))+theme_bw()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),panel.border = element_blank(),panel.grid = element_blank(),axis.line = element_line(size=1,colour = 'black'))

p2=ggplot(pData(HSMM),aes(x=Pseudotime,col=month))+geom_density(size=1)+scale_color_manual(values = rainbow(length(unique(pData(HSMM)$month))))+theme_bw()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),panel.border = element_blank(),panel.grid = element_blank(),axis.line = element_line(size=1,colour = 'black'))

p3=ggplot(pData(HSMM),aes(x=Pseudotime,col=month))+geom_density(col=NA,size=1)+
  scale_color_manual(values = rainbow(length(unique(pData(HSMM)$month))))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),panel.border = element_blank(),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = 'black'))+
  facet_wrap(~month,nrow = length(unique(pData(HSMM)$month)),scales = 'free_y')+
  geom_line(stat='density')+
  theme(strip.background = element_rect(fill = NA,size = 0.5),
        strip.text = element_blank())


p4=ggplot(pData(HSMM),aes(x=Pseudotime,col=month))+geom_density(col=NA,size=1)+
  scale_color_manual(values = rainbow(length(unique(pData(HSMM)$month))))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),panel.border = element_blank(),
                   panel.grid = element_blank(),axis.line = element_line(size=1,colour = 'black'))+
  facet_wrap(~month_version,nrow = length(unique(pData(HSMM)$month_version)),scales = 'free_y')+
  geom_line(stat='density')+
  theme(strip.background = element_rect(fill = NA,size = 0.5),
        strip.text = element_blank())

pdf(paste0(dataname0,"mean_exp",mean_exp,"density.pdf"))
print(p3)
print(p1)
print(p2)
print(p4)
dev.off()


HSMM=HSMM[,order(pData(HSMM)$Pseudotime)]
seuratData_v_rm0326_rmstromal0.1=seuratData
HSMM_v_rm0326_rmstromal0.1=HSMM


pData(HSMM)$combine=sapply(as.character(pData(HSMM)$month),function(x){if(x%in%c("m2","m4")){return("m2_m4")}else if(x%in%c(paste0("m",c(9:13)))){return("m9_m13")}else if(x%in%c(paste0("m",c(15:22)))){return("m15_m22")}else if(x%in%c(paste0("m",c(24:30)))){return("m24_30")}else{return(0)}})
pData(HSMM)$combine=factor(pData(HSMM)$combine,levels = c("m2_m4","m9_m13","m15_m22","m24_30"))
seuratData@meta.data$combine=pData(HSMM)$combine[match(colnames(seuratData),colnames(HSMM))]

##show gene trend along pseudotime

diff_test_res <- differentialGeneTest(HSMM,fullModelFormulaStr = "~sm.ns(Pseudotime)")
save(diff_test_res,file=paste0(dataname0,"HSMM_Pseudotime_diffgene.RData"))
xlsx::write.xlsx(diff_test_res,file = paste0(dataname0,"diff_test_res.xlsx"),sheetName = "diff_test_res")

##genes along with pseudotime

diff_test_res=diff_test_res[diff_test_res$use_for_ordering,]#7707
diff_test_res=diff_test_res[diff_test_res$pval<0.1,]#3640
diff_test_res=diff_test_res[diff_test_res$qval<0.05,]#1872

sig_gene_names <- row.names(diff_test_res)

annotation_col=data.frame(row.names =rownames(pData(HSMM)),month=pData(HSMM)$month,State=pData(HSMM)$State,Phase=pData(HSMM)$Phase)


annot_df_col=pData(HSMM)[,c("month","Phase","Pseudotime")]

column_color2=rainbow(length(unique(annot_df_col$month)))
names(column_color2)=levels(annot_df_col$month)
column_color3=rainbow(length(unique(annot_df_col$Phase)))
names(column_color3)=unique(annot_df_col$Phase)

anno_col=HeatmapAnnotation(df=annot_df_col,col=list(month=column_color2,Phase=column_color3,Pseudotime=circlize::colorRamp2(c(0, max(annot_df_col$Pseudotime)), c("white","darkgreen"))))

heatmap_matrix=plot_pseudotime_heatmap2_return_heatmap_matrix(HSMM[sig_gene_names,],cluster_rows = T,
                                                              add_annotation_col = annotation_col,
                                                              return_heatmap = T,hclust_method = 'ward.D2',
                                                              num_clusters=1)
bks <- seq(-3.1, 3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)


max_ind=apply(heatmap_matrix,1,function(x){a=order(x,decreasing = T)[1:20];xcol=ncol(heatmap_matrix);
s1=sum(sapply(a,function(y){y<0.15*xcol}))
s2=sum(sapply(a,function(y){y>=0.15*xcol&y<0.5*xcol}))
s3=sum(sapply(a,function(y){y>=0.5*xcol&y<0.85*xcol}))
s4=sum(sapply(a,function(y){y>=0.85*xcol}))
return(paste0("cluster",which.max(c(s1,s2,s3,s4))))
})
max_ind=factor(max_ind,levels = paste0("cluster",1:4))
heatmapresult=Heatmap(as.matrix(heatmap_matrix), name = "scaledata", 
                      cluster_rows = T, 
                      cluster_columns = F, 
                      top_annotation = anno_col, use_raster=T,
                      # column_split = annot_df_col$type,
                      row_split= max_ind,
                      # left_annotation = anno_row,
                      show_column_dend = F, show_row_dend = F, 
                      clustering_distance_rows='euclidean', show_row_names = F,show_column_names = F,
                      col=hmcols,
                      #column_dend_height = unit(40, "mm"),#row_names_gp = gpar(fontsize=10),
                      row_km = 1)
tiff(paste0(dataname0,"heatmap.tiff"))#,width = 1000,height = 800,res = 150
h1.1=draw(heatmapresult)
dev.off()



r.dend=column_dend(h1.1)
rc1.list=row_order(h1.1)
names(rc1.list)=paste0(names(rc1.list),".")

c.order=do.call("c",rc1.list)

allcell=data.frame(cell=rownames(heatmap_matrix)[c.order],cellcluster=names(c.order))
allcell$cellcluster2=gsub("[.].*","",allcell$cellcluster)
allcell$newcluster=allcell$cellcluster2

heatmap_matrix1=heatmap_matrix[allcell$cell,]
heatmapresult=Heatmap(as.matrix(heatmap_matrix1), name = "scaledata", cluster_rows = T, 
                      cluster_columns = F, 
                      top_annotation = anno_col, use_raster=T,
                      row_split= allcell$newcluster,
                      show_column_dend = F, show_row_dend = F, 
                      clustering_distance_rows='euclidean', show_row_names = F,show_column_names = F,
                      col=hmcols,
                      column_dend_height = unit(40, "mm"),#row_names_gp = gpar(fontsize=10),
                      row_km = 1)
tiff(paste0(dataname0,"heatmap2.tiff"),width = 1000,height = 800,res = 150)#
h1.1=draw(heatmapresult)
dev.off()

r.dend=column_dend(h1.1)
rc1.list=row_order(h1.1)
names(rc1.list)=paste0(names(rc1.list),".")

c.order=do.call("c",rc1.list)
allcell2=data.frame(cell=rownames(heatmap_matrix1)[c.order],cellcluster=names(c.order))
allcell2$cellcluster2=gsub("[.].*","",allcell2$cellcluster)
allcell2$cellcluster2=factor(allcell2$cellcluster2,levels = paste0("cluster",1:4))
heatmap_matrix2=heatmap_matrix1[allcell2$cell,]

heatmapresult=Heatmap(as.matrix(heatmap_matrix2), name = "scaledata", cluster_rows = F, 
                      cluster_columns = F, 
                      top_annotation = anno_col, use_raster=T,
                      # column_split = annot_df_col$type,
                      row_split= allcell2$cellcluster2,
                      # left_annotation = anno_row,
                      show_column_dend = F, show_row_dend = F, 
                      clustering_distance_rows='euclidean', show_row_names = F,show_column_names = F,
                      # col = c(colorRampPalette(c("grey60", "white"))(10), colorRampPalette(c(RColorBrewer::brewer.pal(9, "YlOrRd")))(60)),
                      col=hmcols,
                      column_dend_height = unit(40, "mm"),#row_names_gp = gpar(fontsize=10),
                      row_km = 1)
tiff(paste0(dataname0,"heatmap3.tiff"),width = 1000,height = 800,res = 150)#
h1.1=draw(heatmapresult)
dev.off()
xlsx::write.xlsx(allcell2,file = paste0(dataname0,"pseudotime_gene_function.xlsx"),sheetName = "gene_cluster",append = T)

