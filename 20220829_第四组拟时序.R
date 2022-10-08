# 拟时序4的横坐标顺序：C1.0;  C1.3; C1.4;  C11.0, C11.1; C8
##存放路径
Idents(expr)=expr$cluster20220822_new
library(monocle)
library(VisiumSpatial)
cluster_list=c('C1.0','C1.3','C1.4','C11.0','C11.1','C8')
sub_expr=subset(expr,idents=cluster_list)
Idents(sub_expr)=factor(Idents(sub_expr),levels=cluster_list)
name='C1.01.31.4_11.011.1_8'
group='orig.ident'
species='human' ### 注释是不跑的，所以这个参数暂时没有什么用处
runmonoclesubcluster(sub_expr, cluster_list, name, group, species)

p1=plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime", 
                        cell_size = 0.5) + theme(axis.title = element_text(face = "bold", 
                                                                           size = 15, colour = "black"), legend.position = "right", 
                                                 legend.text = element_text(face = "bold", size = 13, 
                                                                            colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                           size = 15, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                    size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                           size = 15), title = element_text(size = 15), strip.text = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                                  size = 15))
pdf(paste(name, "_subset_monocle_Pseudotime.pdf", sep = ""), 
    width = 8, height = 6, onefile = F)
print(p1)
dev.off()
png(paste(name, "_subset_monocle_Pseudotime.png", sep = ""), 
    width = 800, height = 600)
print(p1)
dev.off()
p1=plot_cell_trajectory(HSMM_myo, color_by = "clusters", 
                        show_branch_points = F, cell_size = 0.5) + theme(axis.title = element_text(face = "bold", 
                                                                                                   size = 15, colour = "black"), legend.position = "right", 
                                                                         legend.text = element_text(face = "bold", size = 13, 
                                                                                                    colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                   size = 15, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                            size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                   size = 15), title = element_text(size = 15), strip.text = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                                                          size = 15))
pdf(paste(name, "_subset_monocle_cluster.pdf", sep = ""), 
    width = 8, height = 6, onefile = F)
print(p1)
dev.off()
png(paste(name, "_subset_monocle_cluster.png", sep = ""), 
    width = 800, height = 600)
print(p1)
dev.off()
p1=plot_cell_trajectory(HSMM_myo, color_by = "clusters", 
                        cell_size = 0.5) + facet_wrap(~clusters, nrow = 4) + 
  theme(axis.title = element_text(face = "bold", size = 15, 
                                  colour = "black"), legend.position = "right", legend.text = element_text(face = "bold", 
                                                                                                           size = 13, colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                                     size = 15, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                                              size = 10), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                                     size = 10), title = element_text(size = 15), strip.text = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                                                                            size = 15))

pdf(paste(name, "_subset_monocle_cluster_split.pdf", sep = ""), 
    width = 10, height = 7, onefile = F)
print(p1)
dev.off()
png(paste(name, "_subset_monocle_cluster_split.png", sep = ""), 
    width = 1000, height = 700)
print(p1)
dev.off()
p1=plot_complex_cell_trajectory(HSMM_myo, color_by = "clusters", 
                                 show_branch_points = F, cell_size = 0.5) + theme(axis.title = element_text(face = "bold", 
                                                                                                              size = 15, colour = "black"), legend.position = "right", 
                                                                                    legend.text = element_text(face = "bold", size = 13, 
                                                                                                               colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                              size = 15, colour = "black"), strip.text = element_text(face = "bold", 
                                                                                                                                                                                                                      size = 15))

pdf(paste(name, "_subset_monocle_cluster2.pdf", sep = ""), 
    width = 8, height = 6, onefile = F)
print(p1)
dev.off()
png(paste(name, "_subset_monocle_cluster2.png", sep = ""), 
    width = 800, height = 600)
print(p1)
dev.off()
VisiumSpatial:::plot_group_result(HSMM_myo, group, name)
gene_cluster=read.delim('C1.01.31.4_11.011.1_8_gene_clusters.xls',sep="\t",as.is=T,header=T,check.names=F)
for(i in unique(gene_cluster[,1])){
  tmp=rownames(gene_cluster)[gene_cluster[,1]==i]
  tmp_data=data.frame(gene=tmp,cluster=i,stringsAsFactors = F,check.names=F)
  write.table(tmp_data,paste(name,'.genecluster.',i,'.txt',sep=""),
              sep="\t",quote=F,row.names=F,col.names=F)
}
####################20220830第四组拟时序的作图，折线图和拟时间图
# cds <- NA
library(reshape2)
HSMM_myo1=HSMM_myo
cds=HSMM_myo1
cds_subset=HSMM_myo1

newdataB <- data.frame(Pseudotime = pData(cds)$Pseudotime)
cores=30
trend_formula = "~sm.ns(Pseudotime, df=3)"
BranchAB_exprs_test <- genSmoothCurves(cds[,], cores = cores, 
        trend_formula = trend_formula, relative_expr = T, new_data = newdataB)
BranchAB_exprs_test_copy=BranchAB_exprs_test
# BranchAB_exprs_test=BranchAB_exprs_test[sig_gene_names,]

BranchB_exprs <- BranchAB_exprs_test
BranchB_exprs <- log10(BranchB_exprs + 1)
heatmap_matrix <- BranchB_exprs
# heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1, sd) == 0, ] ### 暂时不过滤基因
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix), center = TRUE)) ##对基因做的归一化，不是按照spot
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0
scale_max=3
scale_min=-3
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
# heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
cds_pData <- pData(cds_subset)
cds_fData <- fData(cds_subset)
# colnames(heatmap_matrix)=cds@phenoData@data$original_cell_id
colnames(heatmap_matrix)=colnames(cds)
cds_exprs=melt(heatmap_matrix)
colnames(cds_exprs) <- c("f_id", "Cell", "expectation")
#######20220830,以下2步程序运行比较慢，修改后再次运行
cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
#####修改后的脚本
# x1=cbind(heatmap_matrix,cds_fData[rownames(heatmap_matrix),])
# x2=x1[,setdiff(colnames(x1),'Pseudotime')] ###删除掉pseudotime信息
# x3=x2%>%tibble::rownames_to_column(var='barcodes')
# x4=melt(x3,id.vars=c('orig.ident','clusters','Size_Factor','State','barcodes'))
# x5=data.frame(x4,Pseudotime=cds_pData[x4$barcodes,'Pseudotime'])
# colnames(x5)[6:7]=c('gene_short_name','expectation')
# cds_exprs=x5
########再添加拟时序信息
expr_plotdata=cds_exprs
library(scales)
my_color_palette <- hue_pal()(6)
# colnames(model_expectation) <- colnames(cds_subset)
# color_list=as.list(c('#FB61D7',"#00B6EB","#53B400","#C49A00","#f8766d","#00C094","#A58AFF"))
# color_list=as.list(c(my_color_palette[c(1,9,3,8,7,2,4,6,5)]))
###20220830 genecluster,4,2,1,3,6,5对应的是a,b,c,d,e,f
#factor_colors = dscale(factor(1:6), hue_pal(l = 75))

#color_list=as.list(c("#FF9289","#FF81D7",'#96CA00','#E199FF','#62BCFF','#F5AA00'))
color_list=as.list(c('#00BA38','#F8766D','#F564E3','#619CFF','#B79F00','#00BFC4'))
names(color_list)=letters[1:6]
expr_plotdata_copy=expr_plotdata
row_cluster=read.delim('C1.01.31.4_11.011.1_8_gene_clusters.xls',sep="\t",as.is=T,header=T,check.names=F)
row_cluster=data.frame(genes=rownames(row_cluster),cluster=row_cluster[,1],stringsAsFactors=F)
rownames(row_cluster)=row_cluster[,1]
row_cluster$genes=as.character(row_cluster$genes)
expr_plotdata=expr_plotdata[expr_plotdata$gene_short_name%in%row_cluster$genes,]
expr_plotdata$newgenecluster=row_cluster[expr_plotdata$gene_short_name,'cluster']
expr_plotdata$changegenecluster=plyr::mapvalues(expr_plotdata$newgenecluster,from=c(4,2,1,3,6,5),to=letters[1:6])


### #7CAE00 绿色。#f8766d。红色。#CD9600，土黄色。
levels=c("C1.0",  "C1.3",  "C1.4",  "C11.0", "C11.1", "C8")
expr_figure_line=list()
for(genecluster in letters[1:6]){
  tmp_figure=expr_plotdata[expr_plotdata$changegenecluster==genecluster,]
  tmp_figure$subcluster=tmp_figure$clusters
  # tmp_figure_outs=c()
  # for(i in as.character(levels(tmp_figure$subcluster))){
  #   tmp_figure_outs_1=tmp_figure[tmp_figure$subcluster==i,]
  #   # tmp_figure_outs_1$expectation=limitdata(tmp_figure_outs_1$expectation)
  #   tmp_figure_outs=rbind(tmp_figure_outs,tmp_figure_outs_1)
  # }
  tmp_figure_outs=tmp_figure
  tmp_max_value=aggregate(tmp_figure$expectation,by=list(tmp_figure$clusters,tmp_figure$changegenecluster),median)
  colnames(tmp_max_value)=c('subcluster','genecluster','value')
  tmp_max_value$subcluster=factor(as.character(tmp_max_value$subcluster),levels=levels)
  print(tmp_max_value)
  
  # p=ggplot(tmp_figure_outs,aes(x=subcluster,y=expectation,color=subcluster))+geom_line()+geom_point()
  # p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value))+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+labs(title = paste('genecluster',genecluster,sep=""))+
  # theme(plot.title = element_text(hjust = 0.5))
  # p_line=ggplot(tmp_max_value,aes(x=subcluster,y=value))+geom_line(size = 3,group=1,col="#3DFFC2")+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+labs(title = paste('genecluster',genecluster,sep=""))+
  #   theme(plot.title = element_text(hjust = 0.5))
  # color_list[[genecluster]]
  p=ggplot(tmp_max_value,aes(x=subcluster,y=value))
  p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value),size=2.5,col=color_list[[genecluster]])+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
    labs(title = paste(genecluster,sep="_"),x="",y='RelativeExp')+geom_point(size = 3,col=color_list[[genecluster]])+
    theme(plot.title = element_text(hjust = 0.5,size=30,face='bold',color='black'),text = element_text(size=25,color='black',face='bold'),
          axis.text.x = element_text(size=18,angle=90,hjust=1,color='black',face='bold'),
          axis.text.y=element_text(size=18,color='black',face='bold'))+
    scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2),labels=c('-2.0','-1.0','0','1.0','2.0'))
  expr_figure_line[[paste('genecluster',genecluster,sep="")]]=p
  # expr_figure_line_point_alpha[[paste('genecluster',genecluster,sep="")]]=p_line
}
library(gridExtra)
expr_line=grid.arrange(expr_figure_line[[1]],expr_figure_line[[2]],expr_figure_line[[3]],expr_figure_line[[4]],
                       expr_figure_line[[5]],expr_figure_line[[6]],ncol = 3)
for(width in c(16:20)){
  print(width)
  ggsave(paste('cluster6line_',gsub('\\-','',Sys.Date()),'nobar2width',width,'.png',sep=""),expr_line,width=width,height=11)
  ggsave(paste('cluster6line_',gsub('\\-','',Sys.Date()),'nobar2width',width,'.pdf',sep=""),expr_line,width=width,height=11)
}
###############带bar的折线图暂时用不到
# bar_compute<-function(data){
#   upper=median(data)+sd(data)
#   lower=median(data)-sd(data)
#   return(c(upper,lower))
# }
# 
# expr_figure_line=list()
# for(genecluster in 1:6){
#   tmp_figure=expr_plotdata[expr_plotdata$newgenecluster==genecluster,]
#   tmp_figure$subcluster=tmp_figure$clusters
#   # tmp_figure_outs=c()
#   # for(i in as.character(levels(tmp_figure$subcluster))){
#   #   tmp_figure_outs_1=tmp_figure[tmp_figure$subcluster==i,]
#   #   # tmp_figure_outs_1$expectation=limitdata(tmp_figure_outs_1$expectation)
#   #   tmp_figure_outs=rbind(tmp_figure_outs,tmp_figure_outs_1)
#   # }
#   tmp_figure_outs=tmp_figure
#   tmp_max_value=aggregate(tmp_figure$expectation,by=list(tmp_figure$clusters,tmp_figure$newgenecluster),function(x)cbind(median(x)[1],bar_compute(x)))
#   tmp_max_value=cbind(tmp_max_value[,1:2],tmp_max_value[[3]][,c(1,3,4)])
#   colnames(tmp_max_value)=c('subcluster','genecluster','value','max','min')
#   tmp_max_value$subcluster=factor(as.character(tmp_max_value$subcluster),levels=c("C12","C14.1",'C14.0',"C11.0","C11.1","C8"))
#   print(tmp_max_value)
#   
#   # p=ggplot(tmp_figure_outs,aes(x=subcluster,y=expectation,color=subcluster))+geom_line()+geom_point()
#   # p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value))+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+labs(title = paste('genecluster',genecluster,sep=""))+
#   # theme(plot.title = element_text(hjust = 0.5))
#   # p_line=ggplot(tmp_max_value,aes(x=subcluster,y=value))+geom_line(size = 3,group=1,col="#3DFFC2")+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+labs(title = paste('genecluster',genecluster,sep=""))+
#   #   theme(plot.title = element_text(hjust = 0.5))
#   # color_list[[genecluster]]
#   p=ggplot(tmp_max_value,aes(x=subcluster,y=value))
#   p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value),size=2.5,col='black')+
#     theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
#     labs(title = paste(genecluster,sep="_"),x="",y='RelativeExp')+geom_point(size = 3,col='black')+
#     theme(plot.title = element_text(hjust = 0.5,size=30,face='bold',color='black'),text = element_text(size=25,color='black',face='bold'),
#           axis.text.x = element_text(size=18,angle=90,hjust=1,face='bold',color='black'))+scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2),
#                                                                                                              labels=c('-2.0','-1.0','0','1.0','2.0'))
#   p=p+geom_errorbar(mapping=aes(x=subcluster, ymin=min, ymax=max),width=0.3)
#   
#   expr_figure_line[[paste('genecluster',genecluster,sep="")]]=p
#   # expr_figure_line_point_alpha[[paste('genecluster',genecluster,sep="")]]=p_line
# }
# library(gridExtra)
# expr_line=grid.arrange(expr_figure_line[[1]],expr_figure_line[[2]],expr_figure_line[[3]],expr_figure_line[[4]],
#                        expr_figure_line[[5]],expr_figure_line[[6]],
#                        ncol = 3)
# ggsave('cluster6line_20220314bar.png',expr_line,width=18,height=11)
# ggsave('cluster6line_20220314bar.pdf',expr_line,width=18,height=11)
####################基因折线图拟时序图等20220830
expr_plotdata=cds_exprs
files='genes.txt.txt'
gene_express=read.delim(file.path(getwd(),files),sep="\t",as.is=T,header=T,check.names=F)
c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%gene_express[,1],]
tmp_figure$clusters=factor(tmp_figure$clusters,levels=levels)
########再画基因折线
library(monocle)
HSMM_myo1=HSMM_myo
# C1_genes=read.delim(file.path('C12_8_11_14','genes.txt'),sep="\t",as.is=T,header=T,check.names=F)
# C1_genes=read.delim(file.path('C1_11_8_linecurve.txt'),sep="\t",as.is=T,header=T,check.names=F)
# expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/expr20220824.combin.data.RDS')

x2=pData(HSMM_myo1)
new_expr=as_matrix(expr[['Spatial']]@data)
new_expr=new_expr[,rownames(x2)]

for(i in 1:nrow(gene_express)){
  new_name=gene_express[i,2]
  x2=cbind(x2,new_expr[gene_express[i,1],])
}
colnames(x2)[6:ncol(x2)]=gene_express[1:nrow(gene_express),1]
# for(i in 1:nrow(new_tmp_data)){
#   x2=cbind(x2,new_expr[new_tmp_data[i,'gene'],])
# }
# colnames(x2)[70:ncol(x2)]=new_tmp_data[,'gene']

x2$clusters=factor(x2$clusters,levels=levels)
HSMM_myo1@phenoData@data=x2
c1_trajectorylist=trajectory_function(HSMM_myo1=HSMM_myo1,gene_express)

wrap_plots_combin(linelist=c1_linelist,trajectorylist=c1_trajectorylist,cluster_name='C1.01.31.4_11.011.1_8')
######################重新画图20220830
# /data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_11_8_heatmap_20220401
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8'
library(monocle)
library(ComplexHeatmap)
library(ggplot2)
library(pheatmap)
# load(file.path(indir,'1_14_12_HSMM_myo.Rdata'))
# outpath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_14_12_heatmap_20220401'
# sig_gene_names=read.delim(file.path(indir,'1_14_12_gene_clusters.xls'),sep="\t",as.is=T,header=T)
num_cluster=6
source('/data7/konghl/data/SpatialTranscriptome/software/plot_pseudotime_heatmap.R')
heatmap=plot_pseudotime_heatmap(HSMM_myo[rownames(row_cluster), ], 
                                num_clusters = num_cluster, cores = 30, show_rownames = F,return_heatmap = TRUE)
row_cluster=cutree(heatmap$ph_res$tree_row,num_cluster)
########################准备用complextheatmap画热图

show_rownames=FALSE
row_dist=heatmap$row_dist
hclust_method='ward.D2'
num_cluster=6
annotation_row=heatmap$annotation_row
annotation_col=heatmap$annotation_col

hmcols=heatmap$hmcols
heatmap_matrix=heatmap$heatmap_matrix
exp_rng <- range(heatmap_matrix)
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
#####颜色对应关系  2a 3b 4c 6d 5e 1f 
annotation_colors=list()

annotation_colors[['Cluster']]=c('#00BA38','#F8766D','#F564E3','#619CFF','#B79F00','#00BFC4')
names(annotation_colors[['Cluster']])=letters[1:num_cluster]

library(plyr)
annotation=data.frame(gene=rownames(annotation_row),Cluster=as.character(annotation_row[,1]),
                      stringsAsFactors = F)
annotation[,2]=mapvalues(as.character(annotation[,2]),from=c(4,2,1,3,6,5),to=letters[1:num_cluster])
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
                grid_height = unit(2, "cm"),grid_width = unit(5, "mm"),legend_height = unit(14, "cm"),
                labels_gp= gpar(fontsize = 15*0.8,fontface='bold'),title_gp=gpar(fontsize=15,fontface='bold')
              ),show_heatmap_legend =TRUE,
              row_split = annotation,show_row_names = F,
              left_annotation = row_ha,row_title = NULL,
              column_title = NULL,cluster_row_slices = FALSE, cluster_columns=FALSE,
              show_column_dend=FALSE,show_column_names=FALSE,
              cluster_column_slices = FALSE)
# for(width in c(9,10,11,12)){
width=9
png(paste('C1.01.31.4_11.011.1_8',gsub('\\-','',Sys.Date()),'heatmaplegend.png',sep=""),width=width*100,height=600)
draw(p2,padding=unit(c(1,2,3,4),'cm'))
dev.off()
pdf(paste('C1.01.31.4_11.011.1_8',gsub('\\-','',Sys.Date()),'heatmaplegend.pdf',sep=""),width=width,height=6)
draw(p2,padding=unit(c(1,2,3,4),'cm'))
dev.off()
# }
png(paste('C1.01.31.4_11.011.1_8',gsub('\\-','',Sys.Date()),'heatmaplegendV1.png',sep=""),width=width*100,height=600)
draw(p2,padding=unit(c(1,1,1,3),'cm'))
dev.off()
pdf(paste('C1.01.31.4_11.011.1_8',gsub('\\-','',Sys.Date()),'heatmaplegendV1.pdf',sep=""),width=width,height=6)
draw(p2,padding=unit(c(1,1,1,3),'cm'))
dev.off()

################################基因cluster的注释信息
anno=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/')
row_cluster_anno=data.frame(row_cluster,anno[row_cluster$genes,])
row_cluster_anno$cluster=mapvalues(row_cluster_anno$cluster,from=c(4,2,1,3,6,5),to=letters[1:6])
write.table(row_cluster_anno,'row_cluster_anno_20220830.xls',sep="\t",quote=F,row.names=F,col.names=T)
##################################拟时间label图信息调整20220830
plot_pseudotime_monocle<-function(name,HSMM_myo1,width=12,height=7,legend_height=2){
  p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "Pseudotime", 
                          cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                           size = 30, colour = "black"), legend.position = "right", 
                                                 legend.text = element_text(face = "bold", size = 25, 
                                                                            colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                           size = 30, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                    size = 25,colour='black'), axis.text.y = element_text(face = "bold", colour='black',
                                                                                                                                                                                                                                          size = 25), title = element_text(size = 30), 
                                                 legend.key.height = unit(2,"cm"))+
    guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
    guides(shape = guide_legend(override.aes = list(size = 6.5)))
  pdf(paste(name, "_subset_monocle_Pseudotime.pdf", sep = ""), 
      width = width, height = height, onefile = F)
  print(p1)
  dev.off()
  png(paste(name, "_subset_monocle_Pseudotime.png", sep = ""), width = width*100, height = height*100)
  print(p1)
  dev.off()
  p1=plot_cell_trajectory(HSMM_myo1,x=1,y=2, color_by = "clusters", 
                          show_branch_points = F, cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                                                   size = 30, colour = "black"), legend.position = "right", 
                                                                         legend.text = element_text(face = "bold", size = 25, 
                                                                                                    colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                   size = 30, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                            size = 25,colour='black'), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                                  size = 25,colour='black'), title = element_text(size = 30), 
                                                                         legend.key.height = unit(2,"cm"))+
    guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
    guides(color = guide_legend(override.aes = list(size = 6.5)))
  png(paste(name, "_subset_monocle_cluster.png", sep = ""), width = width*100, height = height*100)
  print(p1)
  # scale_color_manual(values=color_values))
  dev.off()
  pdf(paste(name, "_subset_monocle_cluster.pdf", sep = ""), width = width, height = height,onefile = F)
  print(p1)
  dev.off()
  p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "clusters", 
                          cell_size = 0.5) + facet_wrap(~clusters, nrow = 4) + 
    theme(axis.title = element_text(face = "bold", size = 30, 
                                    colour = "black"), legend.position = "right", legend.text = element_text(face = "bold", 
                                                                                                             size = 25, colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                                       size = 30, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                                                size = 25,colour='black'), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                      size = 25,colour='black'), title = element_text(size = 30), 
          strip.text = element_text(face = "bold", size = 18),legend.key.height = unit(legend_height,"cm"))+
    guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
    guides(color = guide_legend(override.aes = list(size = 6.5)))
  
  png(paste(name, "_subset_monocle_cluster_split.png", sep = ""), width = width*100, height = (height+2)*100)
  print(p1)
  # scale_color_manual(values=color_values))
  dev.off()     
  pdf(paste(name, "_subset_monocle_cluster_split.pdf", sep = ""), width = width, height = height+2, onefile = F)
  print(p1)
  dev.off()
}

# C1.0 PC-a   C1.1 PXy C1.2 PPh-a  C1.3 PC-b    C1.4 MC C1.5 PPh-b   C11.0 CZ 
# 63         29         63        120         83         37         90 
# C11.1 DXy  C14.0 PCL  C14.1 DPh
pData(HSMM_myo1)$clusters=as.character(pData(HSMM_myo1)$clusters)
library(plyr)#"14.0 PL1","14.1 Phloem"
new=c("C1.0 PC-a","C1.3 PC-b", "C1.4 MC", "C11.0 CZ", "C11.1 DXy", "C8 SXy")
old=c("C1.0",  "C1.3",  "C1.4",  "C11.0", "C11.1", "C8")
pData(HSMM_myo1)$clusters=mapvalues(pData(HSMM_myo1)$clusters,from=old,to=new)
###"C1.0"  "C1.2"  "C1.5"  "C12"   "C14.0" "C14.1"
pData(HSMM_myo1)$clusters=factor(pData(HSMM_myo1)$clusters,levels=new)
plot_pseudotime_monocle(name='C1.01.31.4_11.011.1_8',HSMM_myo1,legend_height = 3)
