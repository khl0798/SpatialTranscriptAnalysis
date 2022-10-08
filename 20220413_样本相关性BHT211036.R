#########对数据进行整理
library(stringr)
# sampleNames=c('BHT211036-1-T1','BHT211036-1-T2','BHT211036-1-T3','BHT211036-1-T4','BHT211036-2-N1-2',
#               'BHT211036-2-N1-3','BHT211036-2-T1-2','BHT211036-2-T3-3','BHT211036-3-T1-1',
#               'BHT211036-3-T2-3','BHT211036-3-T3-3','BHT211036-3-T4-1')
sampleNames=c( "BHT211036-1-T1",   "BHT211036-1-T3",   "BHT211036-2-N1-2", "BHT211036-2-T1-2",
               "BHT211036-3-T1-1", "BHT211036-3-T3-3", "BHT211036-1-T2",   "BHT211036-1-T4" , 
              "BHT211036-2-N1-3", "BHT211036-2-T3-3", "BHT211036-3-T2-3", "BHT211036-3-T4-1")
####### revise spot name 
spot_type=c("BHT211036-1_T3-2.csv", "BHT211036-1_T4-1.csv", "BHT211036-2-T3-3.csv",
            "BHT211036-3-T3-3.csv", "BHT211036-3-T4-1.csv")
spot_type_list=list()
spot_type_list[[1]]=c('-3','-BHT211036-1-T3')
spot_type_list[[2]]=c('-4','-BHT211036-1-T4')
spot_type_list[[3]]=c('-4','-BHT211036-2-T3-3')
spot_type_list[[4]]=c('-3','-BHT211036-3-T3-3')
spot_type_list[[5]]=c('-4','-BHT211036-3-T4-1')
new_list=list()
for(i in 1:length(spot_type)){
  tmp=read.delim(spot_type[i],sep=",",as.is=T,header=T,check.names=F)
  tmp=tmp[grep('OW|TW',tmp[,2]),]
  print(dim(tmp))
  tmp[,1]=gsub(spot_type_list[[i]][1],spot_type_list[[i]][2],tmp[,1])
  colnames(tmp)=c('barcode','group')
  new_list[[i]]=tmp
}
new_list_outs=dplyr::bind_rows(new_list)
expr_sub=subset(expr,cells=new_list_outs[,1])
rownames(new_list_outs)=new_list_outs[,1]
expr_sub[['stim_name']]=new_list_outs[colnames(expr_sub),2]                            
###### 修改barcodes的名字
# new.names=colnames(expr)
# for(i in 1:length(sampleNames)){
#   new.names=gsub(paste('-',sampleNames[i],sep=""),paste('-',i,sep=""),new.names)
# }
# expr[['newbarcode']]=new.names
# changecellsname<-function(expr_data){
#   newcells=unlist(lapply(str_split(colnames(expr_data),pattern='-|_'),function(t)t[1]))
#   newcells=paste(newcells,expr_data[['new_group']][,1],sep="_")
#   expr_data=RenameCells(expr_data,new.names=newcells)
#   return(expr_data)
# }
# n1_sam_changename=changecellsname(n1_sam)
# expr2_new_changename=changecellsname(expr2_new)

####### 整体的19个样本的cluster信息
# total_expr=merge(n1_sam_changename,expr2_new_changename)
# ########根据圈图移除不要的spot信息
# select_spot<-function(filename,used_expr1,invert=F,idents){
#   # select_samr1=read.delim('sam-r1-v2.csv',sep=",",as.is=T,header=T,check.names=F)
#   select_samr1=read.delim(filename,sep=",",as.is=T,header=T,check.names=F)
#   select_samr1=unlist(lapply(str_split(select_samr1[,1],pattern='-'),function(t)t[1]))
#   select_samr1=paste(select_samr1,idents,sep="_")
#   rm_samr1=setdiff(colnames(used_expr1)[Idents(used_expr1)==idents],select_samr1)
#   print(paste(idents,length(rm_samr1)))
#   n1_sam_rmsamr1top=subset(used_expr1, cells=rm_samr1, invert=TRUE)
#   return(n1_sam_rmsamr1top)
# }
# total_expr_rmsampart=select_spot(filename='../20220214/sam-r1-v2.csv',used_expr1=total_expr,
#                                  invert=F,idents='SAM-R1')
# 
# total_expr_rmsampart=subset(total_expr_rmsampart,idents=c('IN3-R4'),invert=TRUE)
# 
# total_expr_rmsampart=select_spot(filename='../20220214/IN3R3V4.csv',used_expr1=total_expr_rmsampart,
#                                  invert=F,idents='IN3-R3')
check_na <- function(data){
  tmp_sum=apply(data,1,function(t)sum(t!=0))
  tmp_sum_check_na=data[tmp_sum>=(ncol(data)*0.1),]
  return(rownames(tmp_sum_check_na))
}
#########
# expr_sub_data=data.frame(t(as.matrix(expr_sub$Spatial@data)),stim_name=expr_sub$stim_name)
total_sam_split=SplitObject(expr_sub,split.by='stim_name')
total_rna_list=list()
for(e in 1:length(total_sam_split)){
  rna_tmp=total_sam_split[[e]][['Spatial']]@data
  used_genes_tmp=check_na(rna_tmp)
  total_rna_list[[names(total_sam_split)[e]]]=used_genes_tmp
}
merge_genes<-function(x,y){
  return(intersect(x,y))
}
levels=names(table(expr_sub@meta.data$stim_name))
levels=levels[c(1:6,10:18)]
used_genes_checkna=Reduce(merge_genes,total_rna_list)
# in3_genes_checkna=Reduce(merge_genes,total_rna_list[c('IN3-R1','IN3-R2','IN3-R3')])
# used_genes_checkna=unique(c(used_genes_checkna,in3_genes_checkna))
total_used_expr=expr_sub[used_genes_checkna,]
# total_used_expr=FindVariableFeatures(total_used_expr,nfeatures=6000)
# used_expr=
get_mean<-function(data_expr,genes,filename){
  used_rna=as.matrix(data_expr[['Spatial']]@data)[intersect(genes,rownames(data_expr)),]
  used_gene_means=aggregate(t(as.matrix(used_rna)),by=list(Idents(data_expr)),mean)
  expr.data0.mean=used_gene_means
  rownames(expr.data0.mean)=expr.data0.mean[,1]
  expr.data0.mean=expr.data0.mean[,-1]
  expr.data0.mean=t(expr.data0.mean)
  return(expr.data0.mean)
}
Idents(expr_sub)=expr_sub$stim_name
x1=get_mean(data_expr=expr_sub,genes=used_genes_checkna,filename='rep10')
# x=get_mean(data_expr=total_used_expr,genes=unique(c(VariableFeatures(total_used_expr),in3_genes_checkna)),
#            filename='rmpercent10genes_hvg6000')
# x1=get_mean(data_expr=total_used_expr,genes=VariableFeatures(total_used_expr),
#             filename='rmpercent10genes_hvg6000')
otus_hvg6000=apply(x1,2,scale)
otus_pca=apply(otus_hvg6000,1,scale)
rownames(otus_pca)=colnames(x1)
colnames(outs_pca)=rownames(x1)
library(ade4)
otus_pca.dudi <- dudi.pca(otus_pca,center=TRUE,scale=F,scannf=F,nf=3)
ratio <- inertia.dudi(otus_pca.dudi)
#2d-pca
#pca_pc1 <- as.numeric(sprintf("%0.2f",(ratio$TOT)[1,3]))
pca_pc1 <- floor(otus_pca.dudi$eig[1]*10000/sum(otus_pca.dudi$eig))/100
#pca_pc2 <- as.numeric(sprintf("%0.2f",((ratio$TOT)[2,3]-(ratio$TOT)[1,3])))
pca_pc2 <- floor(otus_pca.dudi$eig[2]*10000/sum(otus_pca.dudi$eig))/100
pca_pc3 <- floor(otus_pca.dudi$eig[3]*10000/sum(otus_pca.dudi$eig))/100
otu_axe_pc123 <- otus_pca.dudi$li[,1:3]
otu_axe_pc123$Sample <- rownames(otu_axe_pc123)
otu_axe_pc123$Group=unlist(strsplit(otu_axe_pc123$Sample,split="-R2|-R1|-R3"))
pca_data=otu_axe_pc123
# pca_data <- merge(otu_axe_pc123,info[,c(1,2)], by="Sample",sort =T)
color_cc <- c_color[1:length(unique(pca_data$Group))]
shape_cc <- c(15, 16, 17, 7, 19, 2, 11, 21, 14, 13, 22, 23, 24, 25, 18)[1:length(unique(pca_data$Group))]
color_cc=c("#FF0000","#6A5ACD","#1E90FF","#FF6F00FF","#DA70D6","#5F9EA0")
p <- ggplot(pca_data,aes(x=Axis1,y=Axis2,label=Sample))
if (length(levels(pca_data$Group)) > 6) {
  p1 <- p + geom_point(size=4)+ aes(color=Group)
  p2 <- p1 + scale_colour_manual(values = color_cc) 
} else {p1 <- p + geom_point(size=4)+ aes(color=Group,shape=Group)
p1 <- p1+xlim(c(-120,140))+ylim(c(-100,50))

p2 <- p1 + scale_colour_manual(values = color_cc) + scale_shape_manual(values = shape_cc)
}
p3 <- p2 + xlab(paste("PC1(",pca_pc1,"%)",sep=''))+ylab(paste("PC2(",pca_pc2,"%)",sep=''))+ geom_text(hjust=0.3,vjust=-1.5,size=4.5, show.legend = F)+geom_vline(xintercept=0,color="grey50",linetype=2)+geom_hline(yintercept=0,color="grey50",linetype=2)+  ggtitle("PCA")
p4 <- p3 + theme_bw() + theme(panel.grid.major=element_line(colour=NA)) + theme(plot.title = element_text(hjust = 0.5, vjust=0.5,size=10))
p4=p4+ggforce::geom_mark_ellipse(aes(color=Group,fill=Group),show.legend=FALSE,label.colour='white',con.colour = "white")
p4=p4+ theme(legend.title = element_text(size=30), #change legend title font size
             legend.text = element_text(size=26),axis.title=element_text(size=16,color='black'),
             axis.text=element_text(size=12,color='black'),
             plot.title = element_text(size = 32, color='black'))+guides(shape = guide_legend(override.aes = list(size = 6)))

ggsave(paste("pca_pc1_pc2_",gsub('-','',Sys.Date()),".hvgV3ellipse.png",sep=""),p4,width=10,height=10)
ggsave(paste("pca_pc1_pc2_",gsub('-','',Sys.Date()),".hvgV3ellipse.pdf",sep=""),p4,width=10,height=10)


otu_axe_pc123 <- otus_pca.dudi$li[,1:3]
sites <- data.frame("Sample"=rownames(otu_axe_pc123), otu_axe_pc123)
write.table(sites,"pca_sites.xls",sep = "\t", quote = F, col.names = T, row.names = F)
ggsave(paste("pca_pc1_pc2_",gsub('-','',Sys.Date()),".nohvg.png",sep=""),p4,width=8,height=7.5)
ggsave(paste("pca_pc1_pc2_",gsub('-','',Sys.Date()),".nohvg.pdf",sep=""),p4,width=8,height=7.5)
###################读取sctransfrom的数据进行分析


#################################整理交集基因
path='E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20220224_重新整理基因表达/'
setwd(path)
library(readxl)
data=read.delim('42950genes.txt',sep="\t",as.is=T,header=T,check.names=F)
genes=readRDS('31917genes_1-n1-1_2-n2-1.RDS')
rownames(data)=data[,1]
data1=data[genes,]
write.table(data1,'31917genes.txt',sep="\t",quote=F,row.names=F,col.names=T)
############################################画热图

otus_pca_convert=t(otus_pca)

levels=names(table(expr_sub@meta.data$stim_name))
# otus_pca_convert=otus_pca_convert[,levels]
# levels=levels[c(16:18,1:15)]
otus_pca_convert=otus_pca_convert[,levels]
otus_pca_convert[otus_pca_convert> 2]=2
otus_pca_convert[otus_pca_convert< -2]=-2
library(ggsci)
library(scales)
c_color=pal_d3('category20')(20)
c_color=c_color[c(2:5,7:20)]

split=factor(colnames(otus_pca_convert),levels=levels)
# png(file='SAM_IN1_IN2_IN3_IN5_IN7.20220415counts_10samples_usedgenes_3378_V1.png',width=1000,height = 1500)
pdf(file='SAM_IN1_IN2_IN3_IN5_IN7.20220415counts_10samples_usedgenes_3378_V1.pdf',width=16,height = 20)
ha = HeatmapAnnotation(
  empty  = anno_empty(border = FALSE, height = unit(10, "mm")),
  Cluster = anno_block(gp = gpar(fontsize=10,fill=c_color[1:18],alpha=0.4),
                       labels = colnames(otus_pca_convert))
)

heatmap_legend_param = list(at = c(-2,-1,0, 1, 2), legend_height = unit(10, "cm"))
p1 <- Heatmap(otus_pca_convert, name = "Exp",column_gap=unit(2,'mm'),
              # col= colorRampPalette(c("#4393C3","white", "#D6604D"))(100),
              column_split = split, top_annotation = ha, 
              row_title = NULL,heatmap_legend_param=heatmap_legend_param,
              column_title = NULL,cluster_row_slices = FALSE, cluster_rows = TRUE,
              show_column_dend=FALSE,show_column_names=FALSE,show_row_names = FALSE,
              cluster_column_slices = FALSE)
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
p2=group_block_anno(1:3, "empty", gp = gpar(fill = c_color[11],alpha=0.4), label = "OW-24H-IN7")
p3=group_block_anno(4:5, "empty", gp = gpar(fill = c_color[12],alpha=0.4), label = "OW-72H-IN7")
p4=group_block_anno(6:8, "empty", gp = gpar(fill = c_color[13],alpha=0.4), label = "TW-24H-IN7")
p5=group_block_anno(9:10, "empty", gp = gpar(fill = c_color[14],alpha=0.4), label = "TW-72H-IN7")
# p6=group_block_anno(13:15, "empty", gp = gpar(fill = color_cc[4],alpha=0.4), label = "IN5")
# p7=group_block_anno(16:18, "empty", gp = gpar(fill = color_cc[5],alpha=0.4), label = "IN7")

dev.off()
##########################