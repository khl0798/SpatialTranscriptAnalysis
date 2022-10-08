library(Seurat)
library("VisiumSpatial")
library(plyr)
library(reshape2)
library(pheatmap)
library(ggplot2)

expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.data.RDS')
c1_3_8_9=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')
savePath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220421_chongxinhuatu'
namelist=list.files(savePath,pattern = '.txt')
cluster1=c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C13','C12','C14.1','C14.0','C11.0','C11.1','C8')
###################################直接读取注释文件
anno=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/anno_tf_go_kegg_ko.RDS')
new_list=list()
gene_id_list=list()
total_expr_cluster=as.character(Idents(expr))
names(total_expr_cluster)=colnames(expr)
total_expr_cluster[colnames(c1_3_8_9)]=paste('C',as.character(Idents(c1_3_8_9)),sep="")
expr[['clusters']]=total_expr_cluster[colnames(expr)]

expr.data0 <- expr[["Spatial"]]@data
expr.data=t(as_matrix(expr.data0))
data = data.frame(expr.data, cluster =expr$clusters)

cluster_sample_gene_mean1 <- aggregate(x = data[, 1:(ncol(data) - 1)], by = list(data$cluster), FUN = mean)
rownames(cluster_sample_gene_mean1) <- cluster_sample_gene_mean1$Group.1
cluster_sample_gene_mean1 <- t(cluster_sample_gene_mean1[, -1])

cluster_sample_gene_mean_scale1=apply(cluster_sample_gene_mean1,1,scale)
cluster_sample_gene_mean_scale1=t(cluster_sample_gene_mean_scale1)
colnames(cluster_sample_gene_mean_scale1)=colnames(cluster_sample_gene_mean1)
new_list=list()
gene_id_list=list()
for(f in namelist){###批量画气泡图
  print(f)
  cellcycle=read.delim(file.path(savePath,f),sep="\t",as.is=T,header=T,check.names = F)
  colnames(cellcycle)=c('gene','geneID')
  cellcycle[,1]=gsub(' ','',cellcycle[,1])
  cellcycle=cellcycle[which(!duplicated(cellcycle[,1])),]
  rownames(cellcycle)=cellcycle[,1]
  name=gsub('.txt','',f)
  print(name)
  used_genes=intersect(cellcycle[,'gene'],rownames(expr))
  cellcycle=cellcycle[used_genes,]
  tmp_gene_id=cellcycle[,'geneID']
  genelist=cellcycle[,'gene']
  
  cluster_sample_gene_mean_scale=cluster_sample_gene_mean_scale1[cellcycle[,'gene'],cluster1]
  ####
  # x1=pheatmap(cluster_sample_gene_mean_scale)
  # x2=x1$tree_row$labels[x1$tree_row$order]
  x2=genelist
  new_list[[name]]=x2
  cluster_sample_gene_mean_scale=cluster_sample_gene_mean_scale[x2,]
  
  cluster_gene_mean=melt(cluster_sample_gene_mean_scale)
  colnames(cluster_gene_mean)=c('gene','cluster','express')
  
  pct.1=expr.data[,x2]
  idents=as.character(Idents(expr))
  names(idents)=colnames(expr)
  pct.1=data.frame(pct.1,cluster=idents[rownames(pct.1)])
  pct.1_outs=aggregate(pct.1[,1:(ncol(pct.1)-1)],by=list(pct.1$cluster),
                       function(x){sum(x!=0)/length(x)})
  pct.1_outs_data=melt(pct.1_outs)
  
  colnames(pct.1_outs_data)=c('cluster','gene','pct')
  
  pct_cluster_gene_mean=merge(pct.1_outs_data,cluster_gene_mean,sort=FALSE,by=c('cluster','gene'))
  
  new_data=pct_cluster_gene_mean
  gene_id=cellcycle[match(new_data[,'gene'],cellcycle[,1],),2]
  new_data_gene=cellcycle[match(x2,cellcycle[,1]),2]
  
  gene_id=factor(gene_id,levels=rev(new_data_gene),labels=rev(new_data_gene))
  
  new_data[,'gene']=gene_id
  new_data[,'cluster']=factor(new_data[,'cluster'],levels=cluster1)
  p<-ggplot(data=new_data,aes(y=gene,x=cluster))+
    geom_point(aes(color=express,size=pct)) +
    scale_size_continuous(range = c(1,5.5))+labs(x="",y="")+
    scale_colour_gradient(low = "grey", high = "Firebrick4")+
    theme_bw()+guides(size = guide_legend(order = 1), fill = guide_legend(order = 2))
  p1 <- p+theme(
    axis.text=element_text(color='black',face='plain'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=18,colour="black",face = 'bold',),
    axis.text.x=element_text(size=18,colour="black"),
    legend.title = element_text(size=18), #change legend title font size
    legend.text = element_text(size=15)
  ) + theme(axis.text.x = element_text(angle = 270,hjust=0)) #hjust = -0.1,vjust = 0.8  hjust=1,vjust=0.2
  p2 <- p+theme(
    axis.text=element_text(color='black',face='plain'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=18,colour="black",family='Arial',face = 'bold',),
    axis.text.x=element_text(size=18,colour="black"),
    legend.title = element_text(size=18), #change legend title font size
    legend.text = element_text(size=15)
  ) + theme(axis.text.x = element_text(angle = 270,hjust=0)) #hjust = -0.1,vjust = 0.8  hjust=1,vjust=0.2
  height=(20*nrow(cellcycle))/50
  if(height<5){
    if(nrow(cellcycle)==5){
      p1=p1+theme(plot.margin=unit(c(6,0.5,6,0.5),'lines'))
      p2=p2+theme(plot.margin=unit(c(6,0.5,6,0.5),'lines'))
      height=5
    }
    if(nrow(cellcycle)==7){
      p1=p1+theme(plot.margin=unit(c(9.5,0.5,8.5,0.6),'lines'))
      p2=p2+theme(plot.margin=unit(c(9.5,0.5,8.5,0.6),'lines'))
      height=7
    }
  }
  
  pdf(file.path(savePath, paste(gsub('-','',Sys.Date()),name,"_17clusters_scale",nrow(cellcycle),"genes_nobold_pointsize18_Firebrick4.pdf", sep = "")),
      width = 15, height = height, onefile = F) #
  print(p1)
  dev.off()
  png(file.path(savePath, paste(gsub('-','',Sys.Date()),name,"_17clusters_scale",nrow(cellcycle),"genes_nobold_pointsize18_Firebrick4.png", sep = "")),
      width = 15 * 100, height =(height) * 100 ) #
  print(p2)
  
  dev.off()
  nn2=8
  nn=length(x2)
  ##### 18.823
  height=ceiling(nn/5)*3
  
  pdf(file.path(savePath, paste(name, "_marker_umap.pdf", sep = "")), width = 20, height =height ,onefile = FALSE)
  p1=FeaturePlot(expr, features = x2, cols = c("grey",  "red"), reduction = "umap",
                 combine = FALSE,order=TRUE,min.cutoff = 'q5',max.cutoff = 'q95')
  for(i in 1:length(p1)){
    p1[[i]]$labels$title = cellcycle[x2,2][i]
  }
  print(CombinePlots(p1,ncol=5))
  dev.off()
  
  
  png(file.path(savePath, paste(name, "_marker_umap.png", sep = "")), width = 20 * 100, height = height * 100)
  print(CombinePlots(p1,ncol=5))
  dev.off()
  
  
  outs=anno[x2,]
  write.table(outs,file.path(savePath,paste(paste(name,"_",length(x2),"genes",sep=""),'.xls',sep="")),
              sep="\t",quote=F,row.names=T,col.names=T)
  gene_id_list[[name]]=cellcycle[match(x2,cellcycle[,1]),]
}
#######################20220422重新画umap图，只在1，11，14，这3个亚群中展示
library(Seurat)
library(ggplot2)
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')
path='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220421_chongxinhuatu/'
namelist=list.files(path,pattern = 'txt$')
savePath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220421_chongxinhuatu/umap'

##### 18.823
for(s in 1:length(namelist)){
  name=namelist[s]
  gene_express=read.delim(file.path(path,namelist[s]),sep="\t",as.is=T,header=T,check.names=F)
  gene_express=gene_express[!duplicated(gene_express[,1]),]
  colnames(gene_express)=c('gene','geneID')
  rownames(gene_express)=gene_express[,1]
  # cluster=c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C11.0','C11.1','C14.0','C14.1')
  used_cells=colnames(expr)[Idents(expr)%in%cluster1]
  noused_cells=setdiff(colnames(expr),used_cells)
  expr1=expr
  expr.data=expr1@assays$Spatial@data
  expr.data[gene_express[,1],noused_cells]=0
  expr1@assays$Spatial@data=expr.data
  nn2=8
  nn=length(gene_express[,1])
  height=ceiling(nn/5)*3
  pdf(file.path(savePath, paste(name, "_marker_umap.pdf", sep = "")), width = 20, height =height ,onefile = FALSE)
  p1=FeaturePlot(expr1, features = gene_express[,1], cols = c("grey",  "red"), reduction = "umap",
                 combine = FALSE,order=TRUE,min.cutoff = 'q5',max.cutoff = 'q95')
  for(i in 1:length(p1)){
    p1[[i]]$labels$title = gene_express[i,2]
  }
  print(CombinePlots(p1,ncol=5))
  dev.off()
  png(file.path(savePath, paste(name, "_marker_umap.png", sep = "")), width = 20*100, height =height*100 )
  print(CombinePlots(p1,ncol=5))
  dev.off()
}
