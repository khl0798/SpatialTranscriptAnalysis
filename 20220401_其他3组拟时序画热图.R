################20220401
#### C1亚群拟时序的结果
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1'
library(monocle)
library(ComplexHeatmap)
library(ggplot2)
load(file.path(indir,'C1_20220118_HSMM_myo.Rdata'))
outpath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1/9clusters_heatmap_chongxinxiugai_20220401'
sig_gene_names=read.delim(file.path(indir,'C1_20220118_gene_clusters9.xls'),sep="\t",as.is=T,header=T)

heatmap=plot_pseudotime_heatmap(HSMM_myo[rownames(sig_gene_names), ], 
                                num_clusters = 9, cores = 30, show_rownames = F,return_heatmap = TRUE)
row_cluster=cutree(heatmap$tree_row,9)
########################准备用complextheatmap画热图

show_rownames=FALSE
row_dist=heatmap$row_dist
hclust_method='ward.D2'
num_clusters=9
annotation_row=heatmap$annotation_row
annotation_col=heatmap$annotation_col
# annotation_colors=heatmap$annotation_colors
# col_gap_ind=heatmap$col_gap_ind
hmcols=heatmap$hmcols
heatmap_matrix=heatmap$heatmap_matrix
exp_rng <- range(heatmap_matrix)
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
library(dendsort)
row_cols=dscale(factor(1:9), hue_pal(l = 75))
#####颜色对应关系  5a 4b 2c 7d 8e 9f 3g 1h 6i
annotation_colors=list()
annotation_colors[['Cluster']]=c("#82B7FF", "#F590FF","#FF80DE","#00D65C","#00D4FF","#00DCBA",
                                 "#AEC500","#FF9289","#F0AD00")
names(annotation_colors[['Cluster']])=letters[1:9]
# row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
# row_dist[is.na(row_dist)] <- 1
library(plyr)
annotation=data.frame(gene=rownames(annotation_row),Cluster=as.character(annotation_row[,1]),
                      stringsAsFactors = F)
annotation[,2]=mapvalues(as.character(annotation[,2]),from=c(5,4,2,7,8,9,3,1,6),to=letters[1:9])
annotation=annotation[order(annotation[,2]),]
rownames(annotation)=annotation[,1]
annotation[,1]=NULL
annotation[,1]=factor(annotation[,1])
heatmap_matrix_new=heatmap_matrix[rownames(annotation),]
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix_new)))/2)
row_dist[is.na(row_dist)] <- 1
hclust_1 <- hclust(row_dist)

row_ha=rowAnnotation(anno=anno_block(gp=gpar(fontsize=80,fill=annotation_colors$Cluster),
                                     labels_gp = gpar(col = "black", fontsize = 20),
                                     labels=unique(annotation[,1])),width=unit(0.9,'cm'))
split=annotation_col
heatmap_legend_param = list(at = seq(-3,3))

dist_fun<-function(mat){
  y1=as.dist((1 - cor(Matrix::t(mat)))/2)
  y1[is.na(y1)]=1
  return(y1)
}
p2 <- Heatmap(heatmap_matrix_new, name = "Exp",col = hmcols,
              column_split = NULL, top_annotation = NULL,
              clustering_distance_rows =function(mat)dist_fun(mat),
              row_gap = unit(1.5, "mm"),
              heatmap_legend_param = list(
                title = "Exp", at = seq(-3,3), 
                labels = c("-3", "-2", "-1",'0','1','2','3'),
                grid_height = unit(2, "cm"),grid_width = unit(5, "mm"),
                labels_gp= gpar(fontsize = 15*0.8,fontface='bold'),title_gp=gpar(fontsize=15,fontface='bold')
              ),
              row_split = annotation,show_row_names = F,
              left_annotation = row_ha,row_title = NULL,
              column_title = NULL,cluster_row_slices = FALSE, cluster_columns=FALSE,
              show_column_dend=FALSE,show_column_names=FALSE,
              cluster_column_slices = FALSE)
png('C1拟时序图20220401.png',width=900,height=600)
draw(p2,padding=unit(c(1,2,3,4),'cm'))
dev.off()
png('C1拟时序图20220401V2.png',width=800,height=600)
draw(p2,padding=unit(c(1,2,3,4),'cm'))
dev.off()

p3 <- Heatmap(heatmap_matrix_new, name = "Exp",col = hmcols,
              column_split = NULL, top_annotation = NULL,
              clustering_distance_rows =function(mat)dist_fun(mat),
              row_gap = unit(1.5, "mm"),
              heatmap_legend_param = list(
                title = "Exp", at = seq(-3,3), 
                labels = c("-3", "-2", "-1",'0','1','2','3'),
                grid_height = unit(1.6, "cm"),grid_width = unit(5, "mm"),
                labels_gp= gpar(fontsize = 15*0.8,fontface='bold'),title_gp=gpar(fontsize=15,fontface='bold')
              ),
              row_split = annotation,show_row_names = F,
              left_annotation = row_ha,row_title = NULL,
              column_title = NULL,cluster_row_slices = FALSE, cluster_columns=FALSE,
              show_column_dend=FALSE,show_column_names=FALSE,
              cluster_column_slices = FALSE)
pdf('C1拟时序图20220401.pdf',width=9,height=6,onefile = FALSE)
draw(p3,padding=unit(c(1,1,1,3),'cm'))
dev.off()
pdf('C1拟时序图20220401V2.pdf',width=8,height=6,onefile = FALSE)
draw(p3,padding=unit(c(1,1,1,3),'cm'))
dev.off()
