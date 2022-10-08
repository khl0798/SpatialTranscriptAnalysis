####设置工作路径
outpath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_samples_cor'
####获取N2-1样本的表达值并检验信息
# n2_1_expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20210824/combin.subcluster.data.RDS')
n2_1_3_5_7=read.csv('N2-1-7-5-3.csv',header=T)
# n2_1_expr[['new_group']]=paste('N2_1_IN',n2_1_3_5_7[colnames(n2_1_expr),2],sep="_")
####获取C4-1，TW样本的表达值及检验信息
# c4_tw=readRDS('/data7/qisy/sc/10x/finished_proj/BHT201004/analysis/combin.data.RDS')
c4_3_5_7=read.csv("C4-1-TW.csv",header=T)
c4_3_5_7[,1]=gsub('-','-1_',c4_3_5_7[,1])
rownames(c4_3_5_7)=c4_3_5_7[,1]
c4=c4_3_5_7%>%filter(C4.1%in%c('3','5','7'))
tw=c4_3_5_7%>%filter(!C4.1%in%c('3','5','7'))
c4[,2]=paste('C4_1_IN',c4[,2],sep="_")
tw[,2]=paste('TW_IN',tw[,2],sep="_")
c4_tw_3_5_7=rbind(c4,tw)
rownames(c4_tw_3_5_7)=c4_tw_3_5_7[,1]
c4_tw=subset(c4_tw,cells=rownames(c4_tw_3_5_7))
c4_tw[['new_group']]=c4_tw_3_5_7[colnames(c4_tw),2]
####获取S2,S3样本的表达信息
s2_s3=readRDS('/data7/guohua/Visium/BHT201004_100g/analysis1/combin.data.RDS')
new_group=s2_s3@meta.data$orig.ident
names(new_group)=colnames(s2_s3)
new_group[grep('1_S2_1',new_group)]='IN3-R3'
new_group[grep('2_S2_2',new_group)]='IN3-R4'
new_group[grep('3_S3_1',new_group)]='IN5-R3'
new_group[grep('4_S3_2',new_group)]='IN5-R4'
s2_s3[['new_group1']]=new_group[colnames(s2_s3)]
######对c4_tw信息整理，分成4组
c4_tw_group=as.character(c4_tw[['new_group']][,1])
names(c4_tw_group)=colnames(c4_tw)
c4_tw_group[grep('C4_1_IN_3',c4_tw_group)]='IN3-R2'
c4_tw_group[grep('C4_1_IN_5',c4_tw_group)]='IN5-R2'
c4_tw_group[grep('C4_1_IN_7',c4_tw_group)]='IN7-R2'
c4_tw_group[grep('TW_IN_TW-7',c4_tw_group)]='IN7-R3'
c4_tw[['new_group1']]=c4_tw_group[colnames(c4_tw)]
#######对n2_1样本的信息整理，分成3组
n2_1_group=as.character(n2_1_expr[['new_group']][,1])
names(n2_1_group)=colnames(n2_1_expr)
n2_1_group[grep('N2_1_IN_3',n2_1_group)]='IN3-R1'
n2_1_group[grep('N2_1_IN_5',n2_1_group)]='IN5-R1'
n2_1_group[grep('N2_1_IN_7',n2_1_group)]='IN7-R1'
###"N2_1_IN_3" "N2_1_IN_7" "N2_1_IN_5"
n2_1_expr[['new_group1']]=n2_1_group[colnames(n2_1_expr)]

############ 3个RDS文件合并在一起试试把
n2_1_rna=n2_1_expr[['Spatial']]@data#####过滤掉不表达的基因
c4_tw_rna=c4_tw[['Spatial']]@data
s2_s3_rna=s2_s3[['Spatial']]@data

check_na <- function(data){
  tmp_sum=apply(data,1,function(t)sum(t==0))
  tmp_sum_check_na=data[tmp_sum!=ncol(data),]
  return(tmp_sum_check_na)
}
n2_1_rna_rmna=check_na(n2_1_rna)
c4_tw_rna_rmna=check_na(c4_tw_rna)
s2_s3_rna_rmna=check_na(s2_s3_rna)
######################求used genes，然后再计算mean值
merge_genes<-function(x,y){
  return(intersect(x,y))
}
#####交集基因28023个
used_genes=Reduce(merge_genes,list(rownames(n2_1_rna_rmna),rownames(c4_tw_rna_rmna),rownames(s2_s3_rna_rmna)))
used_rna=cbind(n2_1_rna_rmna[used_genes,],c4_tw_rna_rmna[used_genes,],s2_s3_rna_rmna[used_genes,])
#####求基因的mean值
used_group_info=c(n2_1_expr[["new_group1"]][colnames(n2_1_rna_rmna),1],
                  c4_tw[['new_group1']][colnames(c4_tw_rna_rmna),1],
                  s2_s3[['new_group1']][colnames(s2_s3_rna_rmna),1])

used_gene_means=aggregate(t(as.matrix(used_rna)),by=list(used_group_info),sum)
###
limitData <- function (data, min = NULL, max = NULL) 
{
  data2 <- data
  if (!is.null(min)) {
    data2[data2 < min] <- min
  }
  if (!is.null(max)) {
    data2[data2 > max] <- max
  }
  return(data2)
}
expr.data0.mean=used_gene_means
rownames(expr.data0.mean)=expr.data0.mean[,1]
expr.data0.mean=expr.data0.mean[,-1]
# expr.data0.mean=expr.data0.mean[-1,]
expr.data0.mean=t(expr.data0.mean)
scaletpm = apply(expr.data0.mean, 1, scale)
expr.data = t(scaletpm)
colnames(expr.data)=colnames(expr.data0.mean)
expr.data <- limitData(expr.data, min = -3, max = 3)
##### 准备好数据画热图
levels=c("IN3-R1", "IN3-R2", "IN3-R3","IN3-R4", "IN5-R1", "IN5-R2", "IN5-R3", 
         "IN5-R4", "IN7-R1", "IN7-R2", "IN7-R3")

split=factor(colnames(sct),levels=levels)
png(file='IN3_IN5_IN7.20220119counts.png',width=1000,height = 1500)
ha = HeatmapAnnotation(
  empty  = anno_empty(border = FALSE, height = unit(10, "mm")),
  Cluster = anno_block(gp = gpar(fontsize=10,fill=colors[1:11]),
                       labels = colnames(expr.data))
)

#split=factor(cell.cluster$newCluster,levels=unique(cell.cluster$newCluster))

# split=factor(rep(c('CD4minus_CD8A_plus',"CD4plus_CD8A_minus",'CD4minus_CD8A_minus'),c(7,7,10)),
#              levels=c('CD4minus_CD8A_plus','CD4plus_CD8A_minus','CD4minus_CD8A_minus'))
# row_ha=rowAnnotation(anno=anno_block(gp=gpar(fontsize=50,
#                                              fill=colors[33:45]),
#                                       labels=unique(find_used_genes_info$type)))
# annotation_legend_param = list(gp = list(title = "foo_top_anno")))
heatmap_legend_param = list(at = c(-3,-1.5,0, 1.5, 3), legend_height = unit(10, "cm"))
p1 <- Heatmap(counts.expr.data, name = "Exp",
              # col= colorRampPalette(c("#4393C3","white", "#D6604D"))(100),
              column_split = split, top_annotation = ha, 
              row_title = NULL,heatmap_legend_param=heatmap_legend_param,
              column_title = NULL,cluster_row_slices = FALSE, cluster_rows = TRUE,
              show_column_dend=FALSE,show_column_names=FALSE,show_row_names = FALSE,
              cluster_column_slices = FALSE)

# seekViewport("annotation_empty_1")
# loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
# seekViewport("annotation_empty_3")
# loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

# grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
#           just = c("left", "bottom"), gp = gpar(fill = "red"))
# grid.text("group 1", x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5)
# 
library(GetoptLong)  # for the function qq()
group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar(fontsize=15)) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  print(loc1)
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  print(loc2)
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5,
              gp = label_gp)
  }
}
print(p1)
p2=group_block_anno(1:4, "empty", gp = gpar(fill = "red"), label = "IN3")
p3=group_block_anno(5:8, "empty", gp = gpar(fill = "#76EE00"), label = "IN5")
p4=group_block_anno(9:11, "empty", gp = gpar(fill = "yellow"), label = "IN7")
# lgd1 = Legend(at = 1:13, legend_gp = gpar(fill = colors[33:45]), title = "term",nr = 1)
# pd = packLegend(lgd1, direction = "horizontal")
# pushViewport(viewport(width = 0.8, height = 0.8))
# grid.rect()
# draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
# draw(pd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
# popViewport()
dev.off()
######################画热图相关性
library(pheatmap)
cor_matrix=cor(expr.data0.mean,method='spearman')
png(file='cor_IN3_IN5_IN7_heatmaptest2.20220118.png')
pheatmap(cor_matrix,show_rownames=T,display_numbers = TRUE,
         show_colnames=T,cluster_cols=F,cluster_rows=F)
dev.off()      
#########################样本重新整合后做批次矫正
#tmp=merge(n2_1_expr,c4_tw)
#merge_expr=merge(tmp,s2_s3)
savePath=getwd()
sample.ident <- c()
expr.list <- list()
sample.Var <- c()
for (i in 1:length(sampleNames)) {
  sampleName <- sampleNames[i]
  print(sampleName)
  expr.list[[sampleName]] <- readRDS(paste0(savePath, "/", sampleName, ".expr.RDS"))
  
  sample.ident <- c(sample.ident, rep(sampleName, dim(expr.list[[sampleName]])[2]))
  sample.Var <- c(sample.Var, VariableFeatures(expr.list[[sampleName]]))
}
sample.ident <- as.factor(sample.ident)
message("[", Sys.time(), "] -----: combine raw matrix data")
suppressWarnings(expr <- merge(expr.list[[1]], expr.list[2:length(expr.list)]))
expr[["sample.ident"]] <- sample.ident
#######样本用harmony方法做批次矫正
expr2=FindVariableFeatures(expr2,selection.method='vst',nfeatures=22299)
expr2=ScaleData(expr2,verbose=FALSE,do.center=FALSE)
expr2=RunPCA(expr2)
#expr1=RunHarmony(expr1,group.by.vars='orig.ident',plot_convergence=TRUE,assay.use='Spatial')
### seurat MNN 方法
id=SplitObject(expr2,split.by='group')
suppressWarnings(expr.anchors <- FindIntegrationAnchors(object.list = id, 
                                                        dims = 1:30, anchor.features = 22299))
# expr1=SelectIntegrationFeatures(expr.anchors,nfeatures=2000)

exprMNN <- IntegrateData(anchorset = expr.anchors, dims = 1:30, verbose = F)
exprMNN <- ScaleData(exprMNN, verbose = FALSE,features = c(VariableFeatures(expr),genes))
DefaultAssay(exprMNN) <- "integrated"
exprMNN[["sample.ident"]] <- sample.ident
bool.plotHVG = F
saveRDS(expr.anchors@anchors, file = file.path(savePath, "anchors.RDS"))

harmony_embeddings<-Embeddings(expr1,'harmony')
# expr <- RunUMAP(expr, dims = 1:pc.use)
# expr <- RunTSNE(expr, dims = 1:pc.use,check_duplicates = FALSE)
expr1=RunTSNE(expr1,reduction='harmony',do.fast=TRUE,dims=1:30)
expr1=RunUMAP(expr1,reduction='harmony',dims=1:30)
library(patchwork)
p1=DimPlot(expr1,group.by = 'orig.ident',reduction = 'pca')+ggtitle('pca')
p2=DimPlot(expr1,group.by = 'orig.ident',reduction = 'harmony')+ggtitle('harmony')
p3=DimPlot(expr1,group.by = 'orig.ident',reduction = 'tsne')+ggtitle('tsne')
p4=DimPlot(expr1,group.by = 'orig.ident',reduction = 'umap')+ggtitle('umap')
# p1+p2+p3+p4
ggsave('test1.png',p1+p2+p3+p4,width=10,height=10)
####样本重新整合后求表达值，然后再做相关性分析
samplelist=c("1_S2_1", "1-N1-1", "2_S2_2", "2-N2-1", "3_S3_1", "3-C4-1", "4_S3_2", "4-TW-2")
###修改细胞的名字
new_names=colnames(expr1)
for(i in 1:8){
  print(c(paste('-1_',i,sep=""),paste('_',samplelist[i],sep="")))
  new_names=gsub(paste('-1_',i,sep=""),paste('_',samplelist[i],sep=""),new_names)
}
expr2=RenameCells(expr1,new.names=new_names)
####整理样本的分群信息
expr2=subset(expr2,orig.ident%in%c("1_S2_1", "2_S2_2", "2-N2-1", "3_S3_1", "3-C4-1", "4_S3_2", "4-TW-2"))
###样本分组信息
c4_tw=read.delim("C4-1-TW.csv",sep=",",as.is=T,header=T,check.names=F)
c4_tw_cellname=c4_tw[,1]
c4_tw_cellname=gsub('-1','_3-C4-1',c4_tw_cellname)
c4_tw_cellname=gsub('-2','_4-TW-2',c4_tw_cellname)
c4_tw[,1]=c4_tw_cellname
c4_tw[(grepl('3-C4-1',c4_tw[,1])&c4_tw[,2]=='3'),2]='IN3-R2'
c4_tw[(grepl('3-C4-1',c4_tw[,1])&c4_tw[,2]=='5'),2]='IN5-R2'
c4_tw[(grepl('3-C4-1',c4_tw[,1])&c4_tw[,2]=='7'),2]='IN7-R2'
c4_tw[(grepl('4-TW-2',c4_tw[,1])&c4_tw[,2]=='TW-7'),2]='IN7-R3'
colnames(c4_tw)=c('barcode','group')
rownames(c4_tw)=c4_tw[,1]
##
n2_1=read.delim('N2-1-7-5-3.csv',sep=",",as.is=T,header=T,check.names=F)
n2_1[,1]=gsub('-2','_2-N2-1',n2_1[,1])
n2_1[n2_1[,2]=='3',2]='IN3-R1'
n2_1[n2_1[,2]=='5',2]='IN5-R1'
n2_1[n2_1[,2]=='7',2]='IN7-R1'
colnames(n2_1)=c('barcode','group')
rownames(n2_1)=n2_1[,1]
##
total_info=data.frame(barcode=colnames(expr2),group=as.character(expr2[['orig.ident']][,1]),stringsAsFactors = F)
rownames(total_info)=total_info[,1]
total_info[total_info[,2]=='1_S2_1',2]='IN3-R3'
total_info[total_info[,2]=='2_S2_2',2]='IN3-R4'
total_info[total_info[,2]=='3_S3_1',2]='IN5-R3'
total_info[total_info[,2]=='4_S3_2',2]='IN5-R4'
total_info[rownames(n2_1),2]=n2_1[,2]
total_info[rownames(c4_tw),2]=c4_tw[,2]
expr2[['group']]=total_info[colnames(expr2),2]
rna=as.matrix(expr2[['Spatial']]@data)
#####使用共同的表达基因

expr2=subset(expr2,group='4-TW-2',invert=TRUE)

check_na <- function(data){
  tmp_sum=apply(data,1,function(t)sum(t!=0))
  tmp_sum_check_na=data[tmp_sum>=(ncol(data)*0.1),]
  return(rownames(tmp_sum_check_na))
}
expr2=subset(expr2,idents='4-TW-2',invert=TRUE)

Idents(expr2)=factor(expr2[['group']][,1],levels=c("IN3-R1", "IN3-R2", "IN3-R3", "IN3-R4", "IN5-R1", "IN5-R2", "IN5-R3", "IN5-R4",
                                                   "IN7-R1", "IN7-R2", "IN7-R3"))
expr2[['group']]=Idents(expr2)
id=SplitObject(expr2,split.by='group')
used_rna_list=list()
for(i in 1:length(id)){
  tmp_rna=id[[names(id)[i]]][['SCT']]@data
  used_rna_list[[names(id)[i]]]=check_na(tmp_rna)
}
######################求used genes，然后再计算mean值
merge_genes<-function(x,y){
  return(intersect(x,y))
}
#####交集基因28023个
used_genes=Reduce(merge_genes,used_rna_list)#####
# expr3=ScaleData(expr2,features=rownames(expr2))
# used_rna=as.matrix(expr2[['Spatial']]@data)[intersect(used_genes,rownames(expr2)),]

used_rna=as.matrix(expr2[['Spatial']]@data)[intersect(VariableFeatures(expr2),rownames(expr2)),]
s2_s3=readRDS('/data7/guohua/Visium/BHT201004_100g/analysis1/combin.data.RDS')
n2_1_expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20210824/combin.subcluster.data.RDS')
used_genes1=intersect(VariableFeatures(s2_s3),VariableFeatures(n2_1_expr))
used_genes1=c(used_genes1,VariableFeatures(expr2))
used_genes1=unique(used_genes1)
used_rna=as.matrix(expr2[['Spatial']]@data)[intersect(used_genes1,rownames(expr2)),]

#########limma包去除批次效应
library(limma)
# new_data <- removeBatchEffect(used_rna, batch = colnames(expr3))

#####求基因的mean值
used_gene_means=aggregate(t(as.matrix(used_rna)),by=list(Idents(expr2)),mean)
expr.data0.mean=used_gene_means
rownames(expr.data0.mean)=expr.data0.mean[,1]
expr.data0.mean=expr.data0.mean[,-1]
# expr.data0.mean=expr.data0.mean[-1,]
expr.data0.mean=t(expr.data0.mean)
# counts=trans(counts_sum)
# expr2.counts=CreateSeuratObject(counts = counts,min.cells=0,min.features=0)
# expr2.counts.nor=NormalizeData(expr2.counts,normalization.method = "LogNormalize")

# expr.data0.mean.norm=apply(expr.data0.mean,1,min_max_norm)
cor(expr.data0.mean,method='spearman')

# expr2 <- RunPCA(expr2, features = VariableFeatures(expr2))
# png(file = file.path(savePath, paste("pca.png",sep = "_")), width = 500, height = 400)
# print(DimPlot(expr2))
# dev.off()

######数据标准化
limitData <- function (data, min = NULL, max = NULL) 
{
  data2 <- data
  if (!is.null(min)) {
    data2[data2 < min] <- min
  }
  if (!is.null(max)) {
    data2[data2 > max] <- max
  }
  return(data2)
}
trans<-function(used_data){
  expr.data0.mean=used_data
  rownames(expr.data0.mean)=expr.data0.mean[,1]
  expr.data0.mean=expr.data0.mean[,-1]
  # expr.data0.mean=expr.data0.mean[-1,]
  expr.data0.mean=t(expr.data0.mean)
  return(expr.data0.mean)
}
expr.data0.mean=used_gene_means
rownames(expr.data0.mean)=expr.data0.mean[,1]
expr.data0.mean=expr.data0.mean[,-1]
# expr.data0.mean=expr.data0.mean[-1,]
expr.data0.mean=t(expr.data0.mean)
# counts=trans(counts_sum)
# expr2.counts=CreateSeuratObject(counts = counts,min.cells=0,min.features=0)
# expr2.counts.nor=NormalizeData(expr2.counts,normalization.method = "LogNormalize")

# expr.data0.mean.norm=apply(expr.data0.mean,1,min_max_norm)
cor(expr.data0.mean,method='spearman')
scale_data<-function(used_data){
  scaletpm = apply(expr.data0.mean, 1, scale)
  # scaletpm=t(expr.data0.mean)
  expr.data = t(scaletpm)
  colnames(expr.data)=colnames(expr.data0.mean)
  expr.data <- limitData(expr.data, min = -2.5, max = 2.5)
  return(expr.data)
}
counts.expr.data=scale_data(as.matrix(expr2.counts.nor@assays$RNA@data))
expr2 <- SCTransform(expr2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)


#####################################所有样本重新整合在一起看一下结果分析
runcellrangercombin<-function (datapath, datapath2, savePath, pc.use = 30){
  library("hdf5r")
  # agg <- read.csv(paste(datapath, "/outs/aggregation.csv", 
  #                       sep = ""), header = T, row.names = 1)
  # sampleNames <- rownames(agg)
  # exprf <- paste(datapath, "/outs/filtered_feature_bc_matrix.h5", 
  #                sep = "")
  # expr.data <- Seurat::Read10X_h5(filename = exprf)
  
  library(ggplot2)
  
  for (i in 1:length(sampleNames)) {
    id = paste("-", i, '^',sep = "")
    # expr.data.new <- expr.data[, grep(id, colnames(expr.data))]
    expr.data.new=readRDS(file.path(datapath2,paste(i,'expr.RDS',sep=".")))
    colnames(expr.data.new) <- gsub(id, "-1", colnames(expr.data.new))
    expr <- Seurat::CreateSeuratObject(counts = expr.data.new, 
                                       project = sampleNames[i], assay = "Spatial")
    expr$slice <- 1
    expr$region <- sampleNames[i]
    # imgpath <- paste(datapath, sampleNames[i], "/outs/spatial/", 
    #                  sep = "")
    # img <- Seurat::Read10X_Image(image.dir = imgpath)
    # Seurat::DefaultAssay(object = img) <- "Spatial"
    # img <- img[colnames(x = expr)]
    # expr[[sampleNames[i]]] <- img
    # plot1 <- VlnPlot(expr, features = "nCount_Spatial", pt.size = 0.1) + 
    #   NoLegend()
    # plot2 <- SpatialFeaturePlot(expr, features = "nCount_Spatial") + 
    #   theme(legend.position = "right")
    # pdf(paste(sampleNames[i], "_Total_UMI_in_spots.pdf", 
    #           sep = ""), width = 8, height = 4, onefile = F)
    # print(plot_grid(plot1, plot2))
    # dev.off()
    # plot1 <- VlnPlot(expr, features = "nFeature_Spatial", 
    #                  pt.size = 0.1) + NoLegend()
    # plot2 <- SpatialFeaturePlot(expr, features = "nFeature_Spatial") + 
    #   theme(legend.position = "right")
    # pdf(paste(sampleNames[i], "_Total_genes_in_spots.pdf", 
    #           sep = ""), width = 8, height = 4, onefile = F)
    # print(plot_grid(plot1, plot2))
    # dev.off()
    # saveRDS(expr, file = file.path(savePath, paste(sampleNames[i], 
    #                                                ".expr.RDS", sep = "")))
  }
  message("[", Sys.time(), "] -----: sample data combination")
  expr.list <- list()
  sample.ident <- c()
  sample.Var <- c()
  for (i in 1:length(sampleNames)) {
    sampleName <- sampleNames[i]
    print(sampleName)
    expr.list[[sampleName]] <- readRDS(paste0(savePath, "/", 
                                              sampleName, ".expr.RDS"))
    sample.ident <- c(sample.ident, rep(sampleName, dim(expr.list[[sampleName]])[2]))
    sample.Var <- c(sample.Var, VariableFeatures(expr.list[[sampleName]]))
  }
  sample.ident <- as.factor(sample.ident)
  message("[", Sys.time(), "] -----: combine raw matrix data")
  suppressWarnings(expr <- merge(expr.list[[1]], expr.list[2:length(expr.list)]))
  expr[["sample.ident"]] <- sample.ident
  
  # expr <- FindVariableFeatures(expr, nfeatures = 6000)
  # expr <- ScaleData(object = expr)
  # expr <- RunPCA(expr, verbose = FALSE)
  # expr <- FindNeighbors(expr, dims = 1:pc.use)
  # expr <- RunUMAP(expr, dims = 1:pc.use)
  # expr <- RunTSNE(expr, dims = 1:pc.use,check_duplicates = FALSE)
  # expr <- run_cluster(expr, sampleNames, savePath, resolution = 0.5, 
  #                     clusterStashName = "default")
  return(expr)
}
############################做样本间的相关性重复
expr.data=expr[['Spatial']]@data
data=read.delim('E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20220106拟时序测试/test1.txt',
                sep="\t",as.is=T,header=T,check.names=F)
data=data[,1:2]
# traverse <- function(a,i,innerl){
#   if(i < (ncol(df))){
#     alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
#     desc <- NULL
#     if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
#     else {
#       for(b in alevelinner) desc <- c(desc,traverse(b,i,innerl))
#       il <- NULL; if(innerl==TRUE) il <- a
#       # newickout <- as.list(desc)
#       # names(newickout)=alevelinner
#       (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
#       # newickout=paste(newickout,'+',sep="")
#     }
#   }
#   else { (newickout <- a) }
# }
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    ##if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    ##else {
    for(b in alevelinner) desc <- c(desc,traverse(b,i,innerl))
    il <- NULL; if(innerl==TRUE) il <- a
    # newickout=as.list(desc)
    # names(newickout)=alevelinner
    (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    ##}
  }else { (newickout <- a) }
  # return(newickout)
}
## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  # print(length(newick))
  # print(class(newick))
  # names(newick)=alevel
  # return(newick)
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}
df=data[1:200,]
outs=df2newick(df,innerlabel=TRUE)

# ddply(.data = res,
.variables = .(Cluster),
.fun = function(df, N) {
  if (length(df$Count) > N) {
    if (any(colnames(df) == "pvalue")) {
      idx <- order(df$pvalue, decreasing=FALSE)[1:N]
    } else {
      ## for groupGO
      idx <- order(df$Count, decreasing=T)[1:N]
    }
    return(df[idx,])
  } else {
    return(df)
  }
},
N=showCategory
)

################模仿数据画进化树的图
setwd("E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/计算相关性/CelltypeEvolution-master/CelltypeEvolution-master/evolution_tree")
HSMM_myo=readRDS('E:\\项目统计(正在做)\\空间转录组\\BHT201004\\BHT201004_20210412售后分析\\计算相关性\\new_HSMM_myo.RDS')
mean_data=read.delim('../../../mean.txt',sep="\t",as.is=T,header=T,row.names = 1)
library(dplyr)
str(mean_data %>% rowSums()%>%mutate(sums!=0))
ego <- enrichGO(gene = gene$ENTREZID, keyType = "ENTREZID", 
                OrgDb = OrgDb, ont = "ALL", pAdjustMethod = "BH", 
                pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
ego2 <- setReadable(ego, OrgDb = OrgDb, "ENTREZID")

HSMM_myo
#hc <- hclust(dist(mean_data_1[,1:8]), "ave")
hc=readRDS('../../../hc.RDS')
library(networkD3)
dendroNetwork(hc, height = 600)
#############################样本数据分析
setwd('E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/计算相关性')
C1=read.delim('C1_diff_anno_go_kegg_ko.xlsx',sep="\t",as.is=T,header=T,check.names=F)
all=read.delim('all_Conserved_go_kegg_ko.xls',sep="\t",as.is=T,header=T,check.names=F)
C1.3_up=C1%>%filter(cluster==3,avg_logFC>0)
C1.3_down=C1%>%filter(cluster==3,avg_logFC<0)
C14_up=all%>%filter(cluster==14,avg_logFC>0)
C14_down=all%>%filter(cluster==14,avg_logFC<0)
C1.3_14_up=intersect(C1.3_up$gene,C14_up$gene)
C1.3_14_down=intersect(C1.3_down$gene,C14_down$gene)
used_genes=c(C1.3_14_up,C1.3_14_down)
write.table(used_genes,'C1.3_14_up_down_inter.txt',sep="\t",quote=F,row.names=F,
            col.names=F)
############################空转分析样本的相关性, 20220118, N2-1,C4-1,S2-1,S2-2,S3-1,S3-2,TW
###先确定样本是否做了normalization分析
##############客户重新分析
expr2_new=subset(expr2,idents=c('IN3-R1','IN3-R2','IN3-R3','IN3-R4', 'IN5-R1','IN5-R2','IN7-R1','IN7-R2','IN5-R3','IN5-R4'))

group=as.character(Idents(expr2_new))
names(group)=colnames(expr2_new)
library(plyr)
group=mapvalues(group,from=c('IN5-R3'),to=c('IN7-R3'))
group=mapvalues(group,from=c('IN5-R4'),to=c('IN5-R3'))
expr2_levels=c("IN3-R1", "IN3-R2", "IN3-R3", "IN3-R4", "IN5-R1", "IN5-R2", "IN5-R3", "IN7-R1",
               "IN7-R2", "IN7-R3")
expr2_new[['new_group']]=factor(group[colnames(expr2_new)],levels=expr2_levels)

#expr2_new[['new_group']]=factor(group[colnames(expr2_new)],levels=sort(apply(expand.grid(c('IN3','IN5','IN7'),c('R1','R2','R3')), 1, paste, collapse="-")))
  # c('IN3','IN5','IN7'),c('R1','R2','R3'),sep="-"))
Idents(expr2_new)=expr2_new[['new_group']]
# expr2_new=FindVariableFeatures(expr2_new,nfeatures=2000)
#expr2_new=FindVariableFeatures(expr2_new,nfeatures=6000)
# used_genes1=intersect(VariableFeatures(s2_s3),VariableFeatures(n2_1_expr))
# 
# used_genes1=c(used_genes1,VariableFeatures(expr2_new))
# used_genes1=unique(used_genes1)
# 
# used_rna=as.matrix(expr2_new[['Spatial']]@data)[intersect(VariableFeatures(expr2_new),rownames(expr2_new)),]

#########limma包去除批次效应
library(limma)
# new_data <- removeBatchEffect(used_rna, batch = colnames(expr3))

#####求基因的mean值
used_gene_means=aggregate(t(as.matrix(used_rna)),by=list(Idents(expr2_new)),mean)
expr.data0.mean=used_gene_means
rownames(expr.data0.mean)=expr.data0.mean[,1]
expr.data0.mean=expr.data0.mean[,-1]
# expr.data0.mean=expr.data0.mean[-1,]
expr.data0.mean=t(expr.data0.mean)
# counts=trans(counts_sum)
# expr2.counts=CreateSeuratObject(counts = counts,min.cells=0,min.features=0)
# expr2.counts.nor=NormalizeData(expr2.counts,normalization.method = "LogNormalize")
# expr.data0.mean.norm=apply(expr.data0.mean,1,min_max_norm)
cor(expr.data0.mean,method='spearman')
###########画pca图
x=prcomp(t(expr.data0.mean))

pca_scores=x$x[,1:2]
pca_scores=data.frame(pca_scores,group=factor(rep(c('IN3','IN5','IN7'),each=3)))
p1=ggplot(pca_scores,aes(PC1,PC2,col=group,shape=group)) + 
  geom_point(size=3) + 
  geom_text(aes(label=rownames(pca_scores)),vjust = "outward") + 
  # geom_hline(yintercept = 0,lty=2,col="red") + 
  # geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="PCA_1",y="PCA_2",title = "PCA analysis")
ggsave('pca1_hvg.png',p1)
expr <- RunPCA(expr2_new,features=VariableFeatures(expr2_new))
png('pca_totalgene.png')
print(DimPlot(expr))
dev.off()
expr2_new_top6000=FindVariableFeatures(expr2_new,nfeatures=6000)
expr2_new_top6000 <- RunPCA(expr2_new_top6000,features=VariableFeatures(expr2_new_top6000))

get_mean<-function(data_expr,genes,filename,group_info){
  used_rna=as.matrix(data_expr[['Spatial']]@data)[intersect(genes,rownames(data_expr)),]
  used_gene_means=aggregate(t(as.matrix(used_rna)),by=list(Idents(data_expr)),mean)
  expr.data0.mean=used_gene_means
  rownames(expr.data0.mean)=expr.data0.mean[,1]
  expr.data0.mean=expr.data0.mean[,-1]
  expr.data0.mean=t(expr.data0.mean)
  # expr.data0.mean.scale=apply(expr.data0.mean,1,scale)
  # rownames(expr.data0.mean.scale)=colnames(expr.data0.mean)
  # x=prcomp(t(expr.data0.mean))
  # pca_scores=x$x[,1:2]
  # pca_scores=data.frame(pca_scores,group=factor(rep(group_info,each=3)))
  # p1=ggplot(pca_scores,aes(PC1,PC2,col=group,shape=group)) + 
  #   geom_point(size=3) + 
  #   geom_text(aes(label=rownames(pca_scores)),vjust = "outward") + 
  #   # geom_hline(yintercept = 0,lty=2,col="red") + 
  #   # geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  #   theme_bw() + theme(legend.position = "none") + 
  #   theme(plot.title = element_text(hjust = 0.5)) + 
  #   labs(x="PCA_1",y="PCA_2",title = paste("PCA",filename, "analysis"))
  # ggsave(paste('pca1_',filename,'.png',sep=""),p1)
  # ggsave(paste('pca1_',filename,'.pdf',sep=""),p1)
  return(expr.data0.mean)
}

x=get_mean(data_expr=expr2_new_top6000,genes=used_genes_checkna,filename='rmpercent10genes',group_info=c('IN3','IN5','IN7'))
###########
IN3_expr=subset(expr2_new,indent=c('IN3-R1','IN3-R2','IN3-R3'))

x=prcomp(t(average_mean))
######################

# rna=as.data.frame(expr2_new_copy[['Spatial']]@data)
# rna_split=split(rna,f=expr2_new_copy[['orig.ident']][,1])

expr2_new_split=SplitObject(expr2_new_copy,split.by='orig.ident')
rna_list=list()
for(e in 1:length(expr2_new_split)){
  rna_tmp=expr2_new_split[[e]][['Spatial']]@data
  used_genes_tmp=check_na(rna_tmp)
  rna_list[[names(expr2_new_split)[e]]]=used_genes_tmp
}
used_genes_checkna=Reduce(merge_genes,rna_list)

expr2_new=expr2_new_copy[used_genes_checkna,]
expr2_new=FindVariableFeatures(expr2_new,nfeatures=6000)
expr2_new_in3=subset(expr2_new,ident=c('IN3-R1','IN3-R2','IN3-R3'))
expr2_new_in3=FindVariableFeatures(expr2_new_in3,nfeatures=2000)

samples_info=c()
for(i in samplelist){
  f=read.delim(paste(i,'_qc/cellManifest-all.txt',sep=""),sep="\t",as.is=T,header=T,check.names=F)
  samples_info=rbind(samples_info,c(i,nrow(f),median(f$nGene)))
}
###############SAM IN1 IN2样本的重复性严重
n1_1=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/1-N1-1.expr.RDS')
n1_2=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT211036-2/analysis/N1-2.expr.RDS')
n1_3=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT211036-2/analysis/N1-3.expr.RDS')
check_spotnumber<-function(expr_data,sam_loupe_info){
  tmp_loupe_info=read.delim(sam_loupe_info,sep=",",as.is=T,header=T,check.names=F)
  used_spot=intersect(colnames(expr_data),gsub('-2','-1',tmp_loupe_info[,1]))
  print(c(ncol(expr_data),nrow(tmp_loupe_info),length(used_spot)))
  rownames(tmp_loupe_info)=tmp_loupe_info[,1]
  return(tmp_loupe_info)
}
n1_3_loupe_info=check_spotnumber(expr_data=n1_3,sam_loupe_info='SAM-N1-3.csv')
n1_1_loupe_info=check_spotnumber(expr_data=n1_1,sam_loupe_info='SAM-N1-1.csv')
n1_2_loupe_info=check_spotnumber(expr_data=n1_2,sam_loupe_info='SAM-N1-2.csv')
###### 
rownames(n1_3_loupe_info)=gsub('-2','-1',n1_3_loupe_info[,1])
n1_1[['new_group']]=n1_1_loupe_info[colnames(n1_1),2]
n1_2[['new_group']]=n1_2_loupe_info[colnames(n1_2),2]
n1_3[['new_group']]=n1_3_loupe_info[colnames(n1_3),2]
n1_sam=merge(n1_1,merge(n1_2,n1_3))
n1_sam_split=SplitObject(n1_sam,split.by='new_group')
sam_rna_list=list()
for(e in 1:length(n1_sam_split)){
  rna_tmp=n1_sam_split[[e]][['Spatial']]@data
  used_genes_tmp=check_na(rna_tmp)
  sam_rna_list[[names(n1_sam_split)[e]]]=used_genes_tmp
}
sam_used_genes_checkna=Reduce(merge_genes,sam_rna_list)
# sam_in_used_genes=intersect(sam_used_genes_checkna,used_genes_checkna)
sam_used_expr=n1_sam[sam_used_genes_checkna,]

###############
sam_in_expr=merge(sam_used_expr,expr2_new)
Idents(sam_in_expr)=sam_in_expr[['new_group']]

# sam_in_expr=merge(sam_used_expr,expr2_new[sam_in_used_genes,])
# levels=names(table(sam_in_expr@meta.data$new_group))
# levels=levels[c(16:18,1:15)]

sam_in_expr@meta.data$new_group=factor(sam_in_expr@meta.data$new_group,levels=levels)
sam_in_expr=FindVariableFeatures(sam_in_expr,nfeatures=6000)
x=get_mean(data_expr=sam_in_expr,genes=VariableFeatures(sam_in_expr),
           filename='rmpercent10genes_hvg6000',group_info=c('SAM','IN1','IN2','IN3','IN5','IN7'))
###########重新计算used genes,20220216, 重新计算数据的值
select_spot<-function(filename,used_expr1,used_expr2,invert=F,idents,rm_id=1){
  # select_samr1=read.delim('sam-r1-v2.csv',sep=",",as.is=T,header=T,check.names=F)
  select_samr1=read.delim(filename,sep=",",as.is=T,header=T,check.names=F)
  if(rm_id==1){
    rm_samr1=setdiff(colnames(used_expr1)[Idents(used_expr1)==idents],select_samr1[,1])
    n1_sam_rmsamr1top=subset(used_expr1,cells=rm_samr1,invert=TRUE)
    sam_in_expr_new=merge(n1_sam_rmsamr1top,used_expr2)
  }
  if(rm_id==2){
    rm_samr1=setdiff(colnames(used_expr2)[Idents(used_expr2)==idents],select_samr1[,1])
    n1_sam_rmsamr1top=subset(used_expr2,cells=rm_samr1,invert=TRUE)
    sam_in_expr_new=merge(used_expr1,n1_sam_rmsamr1top)
  }
  return(sam_in_expr_new)
}
# select_spot(filename='sam-r1-v2.csv',used_expr1=n1_sam,used_expr2=expr2_new,invert=F,idents='SAM-R1',rm_id=1)
used_expr2=expr2_new
select_samr1=read.delim('in3r3v8.csv',sep=",",as.is=T,header=T,check.names=F)
select_samr1[,1]=gsub('-1','_1_S2_1',select_samr1[,1])
rm_samr1=setdiff(colnames(used_expr2)[Idents(used_expr2)=='IN3-R3'],select_samr1[,1])
n1_sam_rmsamr1top=subset(used_expr2,cells=rm_samr1,invert=TRUE)
expr2_rmin3r3=n1_sam_rmsamr1top

#######选择基因，重新分析，重新考虑在重复样本中的分析，
# select_samr1=read.delim('sam-r1-v2.csv',sep=",",as.is=T,header=T,check.names=F)
# n1_sam_rmsam=setdiff(colnames(n1_sam)[Idents(n1_sam)=='SAM-R1'],select_samr1[,1])
# n1_sam_rmsamr1top=subset(n1_sam,cells=n1_sam_rmsam,invert=TRUE)
# n1_sam_rmsamr1=n1_sam_rmsamr1top
sam_in_expr_new=merge(n1_sam_rmsamr1,expr2_rmin3r3)
n1_sam_split=SplitObject(sam_in_expr_new,split.by='new_group')
sam_rna_list=list()
for(e in 1:length(n1_sam_split)){
  rna_tmp=n1_sam_split[[e]][['Spatial']]@data
  used_genes_tmp=check_na(rna_tmp)
  sam_rna_list[[names(n1_sam_split)[e]]]=used_genes_tmp
}
# select_names=sort(names(sam_rna_list))[c(1:6,10:18)]
##################################

sam_used_genes_checkna=Reduce(merge_genes,sam_rna_list)
# check_in3 <- function()
# check_in3=Reduce(merge_genes,sam_rna_list[sort(names(sam_rna_list))[7:9]])
# sam_used_genes_checkna=unique(c(sam_used_genes_checkna,check_in3))
sam_in_expr_new_checkna=sam_in_expr_new[sam_used_genes_checkna,]
# n1_sam_rmsamr1top=subset(used_expr1,cells=rm_samr1,invert=TRUE)
sam_in_expr_new_checkna=FindVariableFeatures(sam_in_expr_new_checkna,nfeatures=6000)
# genes=unique(c(VariableFeatures(sam_in_expr_new_checkna),sam_used_genes_checkna))
levels=c("SAM-R1","SAM-R2", "SAM-R3","IN1-R1", "IN1-R2", "IN1-R3", "IN2-R1", "IN2-R2", "IN2-R3", "IN3-R1", "IN3-R2",
         "IN3-R3", "IN5-R1", "IN5-R2", "IN5-R3", "IN7-R1", "IN7-R2", "IN7-R3")
x=get_mean(data_expr=sam_in_expr_new_checkna,genes=VariableFeatures(sam_in_expr_new_checkna),
           filename='rmsamr1top_rmpercent10genes_hvg6000',group_info=c('SAM','IN1','IN2','IN3','IN5','IN7'))#######先看一下sam人

#######################
# sam_in_expr=FindVariableFeatures(sam_in_expr_new_checkna,nfeatures=6000)
# x=get_mean(data_expr=x=get_mean(data_expr=sam_in_expr,genes=VariableFeatures(sam_in_expr),
#                                 filename='rmsamr1top_rmpercent10genes_hvg6000',group_info=c('SAM','IN1','IN2','IN3','IN5','IN7'))#######先看一下sam人
#            ,genes=VariableFeatures(sam_in_expr),
#            filename='rmsamr1top_rmpercent10genes_hvg6000',group_info=c('SAM','IN1','IN2','IN3','IN5','IN7'))#######先看一下sam人
otus_hvg6000=apply(x,2,scale)
otus_pca=apply(otus_hvg6000,1,scale)
rownames(otus_pca)=colnames(x)
#############################画热图

otus_pca_convert=t(otus_pca)
otus_pca_convert=otus_pca_convert[,levels]
otus_pca_convert[otus_pca_convert> 3]=3
otus_pca_convert[otus_pca_convert< -3]=-3

split=factor(colnames(otus_pca_convert),levels=levels)
png(file='SAM_IN1_IN2_IN3_IN5_IN7.20220216counts.png',width=1000,height = 1500)
pdf(file='SAM_IN1_IN2_IN3_IN5_IN7.20220216counts.pdf',width=16,height = 20)
ha = HeatmapAnnotation(
  empty  = anno_empty(border = FALSE, height = unit(10, "mm")),
  Cluster = anno_block(gp = gpar(fontsize=10,fill=c_color[1:18]),
                       labels = colnames(otus_pca_convert))
)

#split=factor(cell.cluster$newCluster,levels=unique(cell.cluster$newCluster))

# split=factor(rep(c('CD4minus_CD8A_plus',"CD4plus_CD8A_minus",'CD4minus_CD8A_minus'),c(7,7,10)),
#              levels=c('CD4minus_CD8A_plus','CD4plus_CD8A_minus','CD4minus_CD8A_minus'))
# row_ha=rowAnnotation(anno=anno_block(gp=gpar(fontsize=50,
#                                              fill=colors[33:45]),
#                                       labels=unique(find_used_genes_info$type)))
# annotation_legend_param = list(gp = list(title = "foo_top_anno")))

heatmap_legend_param = list(at = c(-3,-1.5,0, 1.5, 3), legend_height = unit(10, "cm"))
p1 <- Heatmap(otus_pca_convert, name = "Exp",
              # col= colorRampPalette(c("#4393C3","white", "#D6604D"))(100),
              column_split = split, top_annotation = ha, 
              row_title = NULL,heatmap_legend_param=heatmap_legend_param,
              column_title = NULL,cluster_row_slices = FALSE, cluster_rows = TRUE,
              show_column_dend=FALSE,show_column_names=FALSE,show_row_names = FALSE,
              cluster_column_slices = FALSE)

# seekViewport("annotation_empty_1")
# loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
# seekViewport("annotation_empty_3")
# loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

# grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
#           just = c("left", "bottom"), gp = gpar(fill = "red"))
# grid.text("group 1", x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5)
# 
library(GetoptLong)  # for the function qq()
group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar(fontsize=15)) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  print(loc1)
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  print(loc2)
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5,
              gp = label_gp)
  }
}
print(p1)
p2=group_block_anno(1:3, "empty", gp = gpar(fill = "red"), label = "SAM")
p3=group_block_anno(4:6, "empty", gp = gpar(fill = "#76EE00"), label = "IN1")
p4=group_block_anno(7:9, "empty", gp = gpar(fill = "yellow"), label = "IN2")
p5=group_block_anno(10:12, "empty", gp = gpar(fill = "green"), label = "IN3")
p6=group_block_anno(13:15, "empty", gp = gpar(fill = "#00BFC4"), label = "IN5")
p7=group_block_anno(16:18, "empty", gp = gpar(fill = "#F564E3"), label = "IN7")

# lgd1 = Legend(at = 1:13, legend_gp = gpar(fill = colors[33:45]), title = "term",nr = 1)
# pd = packLegend(lgd1, direction = "horizontal")
# pushViewport(viewport(width = 0.8, height = 0.8))
# grid.rect()
# draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
# draw(pd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
# popViewport()
dev.off()
