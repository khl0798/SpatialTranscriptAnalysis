library(monocle)

HSMM_myo1 <- setOrderingFilter(HSMM_myo, genes355[,1])
HSMM_myo1 <- reduceDimension(HSMM_myo1, max_components = 5, method = "DDRTree")
HSMM_myo1 <- orderCells(HSMM_myo1,reverse=T)
cell_Pseudotime <- data.frame(pData(HSMM_myo1)$Pseudotime)
rownames(cell_Pseudotime) <- rownames(cell_metadata)
write.table(cell_Pseudotime, file = paste(name, "_subset_cell_Pseudotime.xls", sep = ""), sep = "\t", row.names = T, col.names = F)
Dims = t(monocle::reducedDimS(HSMM_myo1))
colnames(Dims) = c("Component 1", "Component 2")
write.table(Dims, file = paste(name, "_subset_cell_reducedDimS.xls", 
                                sep = ""), sep = "\t", row.names = T, col.names = NA)
save(HSMM_myo1, file = paste(name, "_HSMM_myo.Rdata", sep = ""))
Dims = t(monocle::reducedDimS(HSMM_myo1))[,1:2]
colnames(Dims) = c("Component 1", "Component 2")
write.table(Dims, file = paste(name, "_subset_cell_reducedDimS.xls", 
                               sep = ""), sep = "\t", row.names = T, col.names = NA)
##### pseudotime 图
name='test20220322'
pdf(paste(name, "_subset_monocle_Pseudotime.pdf", sep = ""), 
    width = 10, height = 6, onefile = F)
print(plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "Pseudotime", 
                           cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                            size = 15, colour = "black"), legend.position = "right", 
                                                  legend.text = element_text(face = "bold", size = 18, 
                                                                             colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                            size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                     size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                            size = 15), title = element_text(size = 15), 
                                                  strip.text = element_text(face = "bold", size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
        guides(shape = guide_legend(override.aes = list(size = 5))))

dev.off()
png(paste(name, "_subset_monocle_Pseudotime.png", sep = ""), width = 1000, height = 600)
print(plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "Pseudotime", 
                           cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                            size = 15, colour = "black"), legend.position = "right", 
                                                  legend.text = element_text(face = "bold", size = 18, 
                                                                             colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                            size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                     size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                            size = 15), title = element_text(size = 15), 
                                                  strip.text = element_text(face = "bold", size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
        guides(shape = guide_legend(override.aes = list(size = 5))))

dev.off()
###### cluster split 图
pdf(paste(name, "_subset_monocle_cluster2.pdf", sep = ""), width = 11, height = 7, onefile = F)
print(plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "clusters", 
                           cell_size = 0.5) + facet_wrap(~clusters, nrow = 4) + 
        theme(axis.title = element_text(face = "bold", size = 15, 
                                        colour = "black"), legend.position = "right", legend.text = element_text(face = "bold", 
                                                                                                                 size = 18, colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                                           size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                                                    size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                                           size = 15), title = element_text(size = 15), 
              strip.text = element_text(face = "bold", size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
        guides(color = guide_legend(override.aes = list(size = 5))))
# scale_color_manual(values=color_values))
dev.off()
HSMM_myo1 <- reduceDimension(HSMM_myo1, max_components = 6, method = "DDRTree")
HSMM_myo1 <- orderCells(HSMM_myo1,reverse=T)
png(paste(name, "_subset_monocle_cluster2.png", sep = ""), width = 1100, height = 700)
print(plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "clusters", 
                           cell_size = 0.5) + facet_wrap(~clusters, nrow = 4) + 
        theme(axis.title = element_text(face = "bold", size = 15, 
                                        colour = "black"), legend.position = "right", legend.text = element_text(face = "bold", 
                                                                                                                 size = 18, colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                                           size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                                                    size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                                           size = 15), title = element_text(size = 15), 
              strip.text = element_text(face = "bold", size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
        guides(color = guide_legend(override.aes = list(size = 5))))
# scale_color_manual(values=color_values))
dev.off()                                                                                                                                                                                                                       
#### cluster图
png(paste(name, "_subset_monocle_cluster.png", sep = ""), width = 1100, height = 700)
print(plot_cell_trajectory(HSMM_myo1,x=1,y=2, color_by = "clusters", 
                           show_branch_points = F, cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                                                    size = 15, colour = "black"), legend.position = "right", 
                                                                          legend.text = element_text(face = "bold", size = 18, 
                                                                                                     colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                    size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                             size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                    size = 15), title = element_text(size = 15), 
                                                                          strip.text = element_text(face = "bold",size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
        guides(color = guide_legend(override.aes = list(size = 5))))
# scale_color_manual(values=color_values))
dev.off()
pdf(paste(name, "_subset_monocle_cluster.pdf", sep = ""), width = 11, height = 7,onefile = F)
print(plot_cell_trajectory(HSMM_myo1,x=1,y=2, color_by = "clusters", 
                           show_branch_points = F, cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                                                    size = 15, colour = "black"), legend.position = "right", 
                                                                          legend.text = element_text(face = "bold", size = 18, 
                                                                                                     colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                    size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                             size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                    size = 15), title = element_text(size = 15), 
                                                                          strip.text = element_text(face = "bold",size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
        guides(color = guide_legend(override.aes = list(size = 5))))
# scale_color_manual(values=color_values))
dev.off()

diff_test_res2 <- differentialGeneTest(HSMM_myo1[ordering_genes1, ], 
                                       fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 10)
sig_gene_names <- subset(diff_test_res2, qval < 0.05)[, "gene_short_name"]
sig_gene_names = sig_gene_names[!duplicated(sig_gene_names)]
heatmap=plot_genes_branched_heatmap(HSMM_myo[sig_gene_names,],branch_point = 2,cluster_rows = T,cores=12,return_heatmap = T)
anno=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/anno_tf_go_kegg_ko.RDS')
row_cluster = cutree(heatmap$ph_res$tree_row,6)
row_cluster=data.frame(genes=names(row_cluster),genecluster=row_cluster,anno[names(row_cluster),])
write.table(row_cluster, file = paste(name, "_gene_clusters_branch1.xls", sep = ""), sep = "\t", row.names = T)
for(i in 1:6){
   tmp=row_cluster%>%filter(genecluster==i)
   print(head(tmp[,1:2]))
   write.table(tmp[,1:2],paste('genecluster',i,'.branch1.txt',sep=""),sep="\t",quote=F,row.names=F,col.names=F)
}
 png(file='branch_point2_heatmap.png',width=800,height=600)
 print(heatmap$ph_res)
#################### print(plot_genes_branched_heatmap(HSMM_myo[sig_gene_names,],cluster_rows = T,cores=12)) 此行代码应该暂时用不上
dev.off()
pdf(file='branch_point2_heatmap.pdf',width=8,height=6,onefile = F)
print(heatmap$ph_res)
#################### print(plot_genes_branched_heatmap(HSMM_myo[sig_gene_names,],cluster_rows = T,cores=12))
dev.off()
###########################20220328 重新画分支热图
show_rownames=FALSE
row_dist=heatmap$row_dist
hclust_method='ward.D2'
num_clusters=6
annotation_row=heatmap$annotation_row
annotation_col=heatmap$annotation_col
annotation_colors=heatmap$annotation_colors
col_gap_ind=heatmap$col_gap_ind
hmcols=heatmap$hmcols
heatmap_matrix=heatmap$heatmap_matrix
exp_rng <- range(heatmap_matrix)
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
library(dendsort)

# callback = function(hc, ...){dendsort(hc)}

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
library(plyr)
annotation=data.frame(gene=rownames(annotation_row),Cluster=as.character(annotation_row[,1]),stringsAsFactors = F)
annotation[,2]=mapvalues(as.character(annotation[,2]),from=c(6,2,4,3,5,1),to=letters[1:6])
annotation=annotation[order(annotation[,2]),]
rownames(annotation)=annotation[,1]
annotation[,1]=NULL
annotation[,1]=factor(annotation[,1])

heatmap_matrix_new=heatmap_matrix[rownames(annotation),]

annotation_colors[['Cluster']]=c("#FF9289","#FF81D7",'#96CA00','#E199FF','#62BCFF','#F5AA00')
names(annotation_colors[['Cluster']])=letters[1:6]
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix_new)))/2)
row_dist[is.na(row_dist)] <- 1
hclust_1 <- hclust(row_dist)
x=cutree(row_dist,6)
dend = reorder(as.dendrogram(hclust_1), wts=order(1:nrow(annotation)))
row_cluster <- as.hclust(dend)
num_clusters=6
annotation=outs$annotation
annotation_col=outs$annotation_col
annotation_colors=outs$annotation_colors
col_gap_ind=outs$col_gap_ind
bks=outs$bks
hmcols=outs$hmcols
hclust_method='ward.D2'
row_split=annotation_row

ha = HeatmapAnnotation(
  empty  = anno_empty(border = FALSE, height = unit(10, "mm")),
  Cluster = anno_block(gp = gpar(fontsize=10,fill=annotation_colors$`Cell Type`),
                       labels = "")
)
##unique(annotation_col[,1])
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
# p2=ComplexHeatmap::pheatmap(heatmap_matrix_new,name='Exp',col=hmcols,
#                   gaps_col=col_gap_ind,)
slices=diff(c(0, 101, 200))
annotation_names_col =TRUE
mat=heatmap_matrix_new
border_color =ifelse(nrow(mat) < 100 & ncol(mat) < 100,  "grey60", NA)
top_annotation=HeatmapAnnotation(df = annotation_col[, 
            ncol(annotation_col):1, drop = FALSE], col = annotation_colors[1], 
            show_legend = FALSE, show_annotation_name = FALSE, #annotation_names_col
            annotation_legend_param = list(labels_gp= gpar(fontsize = 15*0.8,fontface='bold'),
                        title_gp=gpar(fontsize=15,fontface='bold'),
                        grid_height = unit(2, "cm"), grid_width = unit(5, "mm")),
            simple_anno_size_adjust = TRUE,
            gp = gpar(col = border_color), annotation_name_gp = gpar(fontsize = 15, 
           fontface = "bold"), simple_anno_size = unit(20, "bigpts"), gap = unit(2, "bigpts"))
           # annotation_legend_param=list('Cell type'=list(fontsize=15)))

p2 <- Heatmap(heatmap_matrix_new, name = "Exp",col = hmcols,
              column_split = rep(seq_along(slices), times = slices), top_annotation = top_annotation,
              clustering_distance_rows =function(mat)dist_fun(mat),
              column_gap = unit(1.5, "mm"),
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
png('拟时序分支修改.png',width=900,height=600)
draw(p2,padding=unit(c(1,2,3,4),'cm'))
# print(p2)
dev.off()

top_annotation=HeatmapAnnotation(df = annotation_col[, 
                ncol(annotation_col):1, drop = FALSE], col = annotation_colors[1], 
               show_legend = TRUE, show_annotation_name = annotation_names_col,
               annotation_legend_param = list(labels_gp= gpar(fontsize = 15*0.8,fontface='bold'),
                title_gp=gpar(fontsize=15,fontface='bold'),
                grid_height = unit(1.4, "cm"), grid_width = unit(5, "mm")),
              simple_anno_size_adjust = TRUE,
              gp = gpar(col = border_color), annotation_name_gp = gpar(fontsize = 15, 
                fontface = "bold"), simple_anno_size = unit(20, "bigpts"), gap = unit(2, "bigpts"))
# annotation_legend_param=list('Cell type'=list(fontsize=15)))

p3 <- Heatmap(heatmap_matrix_new, name = "Exp",col = hmcols,
              column_split = rep(seq_along(slices), times = slices), top_annotation = top_annotation,
              clustering_distance_rows =function(mat)dist_fun(mat),
              column_gap = unit(1.5, "mm"),
              heatmap_legend_param = list(
                title = "Exp", at = seq(-3,3), 
                labels = c("-3", "-2", "-1",'0','1','2','3'),
                grid_height = unit(1.4, "cm"),grid_width = unit(5, "mm"),
                labels_gp= gpar(fontsize = 15*0.8,fontface='bold'),title_gp=gpar(fontsize=15,fontface='bold')
              ),
              row_split = annotation,show_row_names = F,
              left_annotation = row_ha,row_title = NULL,
              column_title = NULL,cluster_row_slices = FALSE, cluster_columns=FALSE,
              show_column_dend=FALSE,show_column_names=FALSE,
              cluster_column_slices = FALSE)

p4 <- Heatmap(heatmap_matrix_new, name = "Exp",col = hmcols,
              column_split = rep(seq_along(slices), times = slices), 
              top_annotation = top_annotation,
              clustering_distance_rows =function(mat)dist_fun(mat),
              column_gap = unit(1.5, "mm"),
              heatmap_legend_param=list(),
              show_heatmap_legend = F,
              # heatmap_legend_param = list(
              #   title = "Exp", at = seq(-3,3), 
              #   labels = c("-3", "-2", "-1",'0','1','2','3'),
              #   grid_height = unit(1.4, "cm"),grid_width = unit(5, "mm"),
              #   labels_gp= gpar(fontsize = 15*0.8,fontface='bold'),
              #   title_gp=gpar(fontsize=15,fontface='bold')
              # ),  show_heatmap_legend = FALSE,
              row_split = annotation,show_row_names = F,
              left_annotation = row_ha,row_title = NULL,
              column_title = NULL,cluster_row_slices = FALSE, cluster_columns=FALSE,
              show_column_dend=FALSE,show_column_names=FALSE,
              cluster_column_slices = FALSE)

# p4=p3+NoLegend()
pdf('拟时序分支修改test20220826V1.pdf',width=9,height=6,onefile=FALSE)
# print(p4)
draw(p4,padding=unit(c(1,1,1,3),'cm'))
dev.off()

###########20220826,拟时序分支热图重新出图
for(i in 9:12){
  print(i)
  pdf(paste('C12_11_14_8拟时序分支修改20220826width',i,'.pdf',sep=""),width=i,height=6,onefile=FALSE)
  draw(p3,padding=unit(c(1,1,1,3),'cm'))
  dev.off()
}


# grid::grid.rect(gp = grid::gpar("fill", col = NA))
# grid::grid.draw(ph_res$gtable)

############################# 重新转录因子分析tf20220322,共计3477个基因画分支热图
#服务器存储路径/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220310_monocle_12_11_8_14_adjust12_8_orientation/20220321
#本地存储路径E:\项目统计(正在做)\空间转录组\BHT201004\BHT201004_20210412售后分析\20220321_拟时序结果C12_11_14_8\
##########20220318,C12,C14,C11.0,C11.1,C8 思路：拟时间的分支图，6个genecluster，
#######转换每个基因对应的拟南芥的id，在和拟南芥的TF取交集计算
#######生成新的基因id，每个TF vs 其他的TF，然后遍历，分析explicit数据信息
############对6个genecluster中的1个tf基因vs 其他的tf基因分析
### genecluster path路径
# indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20211116_monocle_C12_14_11.0_11.1_8'
###########此处是C1的9个亚群的信息
# indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1/C1_gene_curveline'
# savePath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1/C1_gene_curveline'
#############此处是C1.0，1.4，11.0，11.1，8的tf分析
# indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_11_8_curveline'
# savePath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_11_8_curveline'
#############此处是C1.0，C1.2，C1.5，C14.0，C14.1，C12的tf分析
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17cluster_20220322_monocle_C12_14_11_8_tf_3477genes'
savePath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17cluster_20220322_monocle_C12_14_11_8_tf_3477genes'
anno=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/anno_tf_go_kegg_ko.RDS')
###杨树基因ID转换拟南芥的基因ID
########根据marker基因取交集
name='c14_12_11_8'
data_anno=read.delim(file.path(indir,'out_branch2_anno.xls'),sep="\t",as.is=T,header=T,check.names=F,row.names=1)
# rownames(data_anno)=data_anno[,1]
data_anno=data.frame(gene=rownames(data_anno),cluster=data_anno$genecluster,ID=anno[rownames(data_anno),'newGeneID'])

###根据拟南芥的ID过滤TF
#读取拟南芥TF
tf_info=read.delim('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/explicit/test_data/Ath_TF_list.txt',
                   sep="\t",as.is=T,header=T,check.names=F)
tf_id=tf_info[,2]
##########过滤
# data_anno$gene
filter_gene_tf=data_anno[data_anno[,'ID']%in%tf_id,]
filter_gene_tf=cbind(filter_gene_tf,anno[as.character(filter_gene_tf[,'gene']),])
filter_gene_tf$gene=as.character(filter_gene_tf$gene)
filter_gene_tf$ID=as.character(filter_gene_tf$ID)
###### 基因cluster信息修改
library(plyr)
filter_gene_tf$cluster=mapvalues(filter_gene_tf$cluster,from=1:6,to=letters[1:6])

write.table(filter_gene_tf,file.path(savePath,'filter_tf_noFCfilter.xlsx'),
            sep="\t",quote=F,row.names=F,col.names=T)

expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')

expr=subset(expr,idents=c("8","11.0","11.1","14.0","14.1","12"))

data_anno_tf_filter=filter_gene_tf
rna=expr[['Spatial']]@data
rna_filter=rna[rownames(rna)%in%data_anno_tf_filter$gene,]
data_anno_tf_filter$new_id=paste(data_anno_tf_filter[,'gene'],data_anno_tf_filter[,'cluster'],data_anno_tf_filter[,'ID'],sep="_")
###### 把cluster1-6转为a-f

rna1=rna_filter
rownames(data_anno_tf_filter)=data_anno_tf_filter[,'gene']
rownames(rna1)=data_anno_tf_filter[rownames(rna1),"new_id"]   ####258个杨树基因,216个拟南芥的转录因子
###match 拟南芥的基因ID信息 
# rna1=rna[rownames(rna)%in%c1[,2],]
# anno1=anno[anno[,1]%in%rownames(rna1),]
# rna1=rna1[match(anno1[,1],rownames(rna1)),]
# rna1=data.frame(rna1,geneID=anno1[match(rownames(rna1),anno1[,1]),'newGeneID'])
# rna1=rna1[!duplicated(rna1$geneID),]  ### 
out=c()
for(i in 1:nrow(rna1)){
  # tmp_rna=rna1[,-ncol(rna1)]
  tmp_tf_name=rownames(rna1)[i]
  print(i)
  tmp_target_name=rownames(rna1)[-i]
  tmp_tf_data=t(data.frame(rna1[i,]))
  rownames(tmp_tf_data)=tmp_tf_name
  tmp_target_data=rna1[-i,]
  rownames(tmp_target_data)=tmp_target_name
  tmp_out=explicit(tmp_target_data,tmp_tf_data,tmp_target_name,tmp_tf_name,path)
  # out[[tmp_tf_name]]=tmp_out
  colnames(tmp_out)=c('target','tf','beta','beta_pvalue')
  out=rbind(out,tmp_out)
}
library(data.table)
out_copy=out
out[,c('target_gene','target_cluster','target_ID')]=do.call(cbind,tstrsplit(out[,'target'], "_", fixed=TRUE))
out[,c('tf_gene','tf_cluster','tf_ID')]=do.call(cbind,tstrsplit(out[,'tf'], "_", fixed=TRUE))


target_anno=anno[out$target_gene,]
colnames(target_anno)=paste('target',colnames(target_anno),sep="_")
tf_anno=anno[out$tf_gene,]
colnames(tf_anno)=paste('tf',colnames(tf_anno),sep="_")

# results=data.frame(out,target_anno,tf_anno)
library(dplyr)
# library(tidyverse)
results=out
results_pvalue_filter=results%>%filter(beta_pvalue<0.05)
results_pvalue_filter_beta_pvalue=results_pvalue_filter%>%filter(abs(beta)>0.2)  #### 3034 个转录因子关系对
# write.table(results,'results.xlsx',sep="\t",quote=F,row.names=F,col.names=T)
write.table(results_pvalue_filter,file.path(savePath,'results_pvalue_filter.xlsx'),sep="\t",quote=F,row.names=F,col.names=T)
results_pvalue_filter_beta_pvalue=results_pvalue_filter_beta_pvalue%>%mutate(up_down=ifelse(beta>0,'up','down'))

write.table(results_pvalue_filter_beta_pvalue,file.path(savePath,paste(name,'_beta0.2_new.xlsx',sep="")),sep="\t",quote=F,row.names=F,col.names=T)
write.table(results_pvalue_filter_beta_pvalue,file.path(savePath,paste(name,'_beta0.2_new.txt',sep="")),sep="\t",quote=F,row.names=F,col.names=T)
results_pvalue_filter_beta_pvalue_cytoscape=results_pvalue_filter_beta_pvalue[,c('target_gene','target_cluster','tf_gene','tf_cluster','beta','beta_pvalue')]
results_pvalue_filter_beta_pvalue_cytoscape$target_gene=gsub(':',"|",anno[results_pvalue_filter_beta_pvalue_cytoscape$target_gene,'genesymbol'])
results_pvalue_filter_beta_pvalue_cytoscape$tf_gene=gsub(':',"|",anno[results_pvalue_filter_beta_pvalue_cytoscape$tf_gene,'genesymbol'])
results_pvalue_filter_beta_pvalue_cytoscape=results_pvalue_filter_beta_pvalue_cytoscape%>%mutate(up_down=ifelse(beta>0,'up','down'))
write.table(results_pvalue_filter_beta_pvalue_cytoscape,file.path(savePath,paste(name,'_beta0.2_new.txt',sep="")),sep="\t",quote=F,row.names=F,col.names=T)

target_anno=unlist(lapply(results_pvalue_filter_beta_pvalue_cytoscape$target_gene,function(t){unlist(strsplit(t,split="\\|"))[1]}))
target_anno=anno[target_anno,]
colnames(target_anno)=paste('target',colnames(target_anno),sep="_")
tf_anno=unlist(lapply(results_pvalue_filter_beta_pvalue_cytoscape$tf_gene,function(t){unlist(strsplit(t,split="\\|"))[1]}))
tf_anno=anno[tf_anno,]
colnames(tf_anno)=paste('tf',colnames(tf_anno),sep="_")
beta_results=data.frame(results_pvalue_filter_beta_pvalue_cytoscape,target_anno,tf_anno)
write.table(beta_results,'beta0.2_anno.xls',sep="\t",quote=F,row.names=F,col.names=T)
###
library(scales) ## 记录cluster的6个颜色
col2rgb=c("#FF9289","#FF81D7",'#96CA00','#E199FF','#62BCFF','#F5AA00')
##### 统计节点信息
get_node_info_up_down<-function(data,geneinfo,name){
  tmp=data[,c('target_gene','tf_gene','target_cluster','tf_cluster')]
  used_gene=data.frame(gene=tmp[,1],cluster=tmp[,c(3)])
  used_gene=rbind(used_gene,data.frame(gene=tmp[,2],cluster=tmp[,c(4)]))
  used_gene=used_gene[!duplicated(used_gene),]
  write.table(used_gene,paste(name,'.geneclusterinfo.txt',sep=""),sep="\t",quote=F,row.names=F,col.names=T)
}
get_node_info_up_down(data=results_pvalue_filter_beta_pvalue_cytoscape,geneinfo=gene_info,name=name)

###################统计node信息
cytoscaple=results_pvalue_filter_beta_pvalue[,1:10]
up_down=rep('up',nrow(cytoscaple))
up_down[cytoscaple$beta<0]='down'
cytoscaple=cytoscaple[,c(5,8,3:4,6:7,9:10)]
cytoscaple$up_down=up_down
colnames(cytoscaple)[1:2]=c('target','tf')
write.table(cytoscaple,file.path(savePath,'C12_14_11_8_beta0.15_cytoscape.txt'),sep="\t",quote=F,row.names=F,col.names=T)

setwd('E:\\项目统计(正在做)\\空间转录组\\BHT201004\\BHT201004_20210412售后分析\\20211124_2组数据网络图\\重新做网络图小图20211202\\C12_14_11.0_11.1_8')
genecluster=read.delim('branch_genecluster6_abcdef.txt',sep="\t",as.is=T)
select_info=read.delim('select_info.txt',sep="\t",as.is=T)
new_id=rep('noselect',nrow(genecluster))
new_id[genecluster$gene%in%select_info$locusName]='select'
genecluster$select_info=new_id
write.table(genecluster,'select_info_new.txt',sep="\t",quote=F,row.names=F,col.names=T)
select_info$cluster=genecluster[select_info$locusName,]$cluster
###################获得基因的注释信息
path='E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20220208_4组拟时序结果整理/tf'
setwd(path)
library(stringr)
x=list.files(getwd())
for(i in x){
  print(f)
  f=list.files(file.path(getwd(),i),pattern='node',full.names = T)
  if(file.exists(f)){
    tmp=read.delim(f,sep=",",as.is=T,header=T,check.names=F)
    gene_name=unlist(lapply(str_split(tmp$name,pattern='[|]'),function(t)t[1]))
    tmp=data.frame(tmp,anno[gene_name,])
    write.table(tmp,gsub('node','nodeanno',f),sep=",",quote=F,row.names=F,col.names=T)
  }
}
