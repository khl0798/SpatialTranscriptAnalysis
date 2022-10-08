########################################基因的小提琴图等展示信息20220822，6张txt文档画小提琴图
library(Seurat)
library(cowplot)
####/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220530_vlnplot
# expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')
library(ggplot2)
#Idents(expr)=expr$clusters
#####20220822小提琴图做图： C1.0;  C1.1; C1.2;  C1.3; C1.4; C1.5;   C12; C14.1;C14.0, C11.0, C11.1; C8				
#####20220824 小提琴图做图： C1.0;  C1.1; C1.2;  C1.3; C1.4; C1.5;   C12; C14.1;C14.0, C11.0, C11.1; C8
##### （其中只在C1.0，C1.4，C14.0, C11.0做小提琴图）

cluster20220822=expr$new20220815
library(plyr)
cluster20220822_new=mapvalues(cluster20220822,from=c('C1.0 PC-a', 'C1.1 PXy', 'C1.2 PPh-a',
                                                     "C1.3 PC-b", "C1.4 MC","C1.5 PPh-b", "C11.0 CZ", 
                                                     "C11.1 DXy", 'C12','C14.0 PCL', 'C14.1 DPh','C8'),
                              to=c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C11.0','C11.1','C12','C14.0','C14.1','C8'))
cluster20220822_new[!cluster20220822_new%in%c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C11.0','C11.1','C12','C14.0','C14.1','C8')]='null'
expr[['cluster20220822_new']]=cluster20220822_new
Idents(expr)=expr$cluster20220822_new
# Idents(expr)=factor(paste('C',as.character(expr$default),sep=""),levels=paste('C',0:16,sep=""))
# gene1=read.delim('17clustersVlnplot20220601.txt',sep="\t",as.is=T,header=T,check.names=F)
#cluster1=c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C12','C14.1','C14.0','C11.0','C11.1','C8')
# cluster1=paste('C',0:16,sep="")
#color_list=readRDS('../17clusters_20220527_qipaotuqiepiantu/color_list.20220528.RDS')
#17clusters颜色信息是什么呢在哪里呢？20220530color_list.RDS
sub_expr=subset(expr,idents=c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C11.0','C11.1','C12','C14.0','C14.1','C8'))
# cluster_color_list
########################先获取所有的基因的小提琴的data信息，后面可以直接获取信息等等
tmp1=FetchData(sub_expr,vars = rownames(sub_expr))
library(dplyr)
# tmp=FetchData(sub_expr,vars = gene1[,1])
library(reshape2)
x=tmp1%>%tibble::rownames_to_column('barcodes')%>%melt('barcodes')
x$cluster=Idents(expr)[x[,1]]
x$cluster=factor(as.character(x$cluster),levels=c("C1.0",  "C1.1", "C1.2",  "C1.3", "C1.4", "C1.5",
                                                  "C12", "C14.1","C14.0", "C11.0", "C11.1", "C8"))
# levels=c("C1.0",  "C1.1", "C1.2",  "C1.3", "C1.4", "C1.5",
#          "C12", "C14.1","C14.0", "C11.0", "C11.1", "C8")
#################20220825修改levels
levels=c("C1.5","C1.2", 'C1.0',  "C1.1",  "C1.3", "C1.4",    "C12", "C14.1","C14.0", "C11.0", "C11.1", "C8")

#小提琴图的颜色高亮展示每个基因在C14.1, C14.0, C11.0, C11.1 中的表达情况，颜色就用各自cluster的颜
colnames(x)[2]='gene'
x=data.frame(x,new='new',vlnplot_value=0)
x$cluster=as.character(x$cluster)
x$new=as.character(x$new)
cluster9_cor_new=c("#99CC00FF","#FFCD00FF","#6EE2FFFF",'#e7de2c',"#7bf1a8","#188cbe",
                   "#FFBFE5" , "#d1b1cb","#2196f3","#38b000", "#fcbc5d","#F79D1E")
names(cluster9_cor_new)=c("C1.0",  "C1.1", "C1.2",  "C1.3", "C1.4", "C1.5",
                          "C12", "C14.1","C14.0", "C11.0", "C11.1", "C8")
# for(i in 1:length(cluster_color_list)){
#   select_id_tmp=x$gene==names(cluster_color_list)[i]&x$cluster%in%cluster_color_list[[i]]
#   x$vlnplot_value[select_id_tmp]=x$value[select_id_tmp]
#   x$new[select_id_tmp]=x$cluster[select_id_tmp]
# }
x$gene=as.character(x$gene)

namelist=list.files(getwd(),pattern='*txt$')

for(name in namelist){
  print(name)
  gene1=read.delim(name,sep="\t",as.is=T,header=T,check.names=F)
  gene1=gene1[!duplicated(gene1[,1]),]
  rownames(gene1)=gene1[,1]
  head(gene1)
  # if(name!='vlnplot20220824_sheet1V1.txt'){
  # genename=do.call(rbind,strsplit(gene1[,2],split="\\|"))
  # gene1$genename=genename[,2]
  # }else{
  #   gene1$genename=gene1[,2]
  # }
  #gene1[,2]=paste(gene1[,1],gene1[,2],sep="|")
  #gene_cluster_color=read.delim('gene3_20220525_showcluster.txt',sep="\t",as.is=T,header=T,check.names=F)
  # cluster_color_list=split(gene1[,3],f=gene1[,1])
  # cluster_color_list=lapply(cluster_color_list,function(t)unlist(stringr::str_split(t,pattern=",")))
  ####C12, C14.1, C14.0, C11.0, C11.1, C8
  #sub_expr=subset(expr,idents=cluster1)
  #Idents(sub_expr)=factor(Idents(sub_expr),levels=cluster1)

  # for(i in 1:length(cluster_color_list)){
  #   select_id_tmp=x$gene==names(cluster_color_list)[i]&x$cluster%in%cluster_color_list[[i]]
  #   x$vlnplot_value[select_id_tmp]=x$value[select_id_tmp]
  #   x$new[select_id_tmp]=x$cluster[select_id_tmp]
  # }
  ########选择基因数据信息
  vlnplot_tmp=x[x$gene%in%gene1[,1],]
  vlnplot_tmp$gene=as.character(vlnplot_tmp$gene)
  
  vlnplot_tmp$newgene=factor(gene1[vlnplot_tmp$gene,3],levels=unique(gene1$genename))
  vlnplot_tmp$cluster=factor(vlnplot_tmp$cluster,levels=levels)
  vlnplot_tmp$vlnplot_cluster=vlnplot_tmp$cluster
  # vlnplot_tmp$vlnplot_cluster[!vlnplot_tmp$vlnplot_cluster%in%c("C1.0","C1.4","C14.0","C11.0")]='new'
  
  ######读取数据
  ######vlnplot_value列控制小提琴图的生成，值为0即不会生成小提琴，
  #小提琴图的颜色高亮展示每个基因在C14.1, C14.0, C11.0, C11.1 中的表达情况，颜色就用各自cluster的颜色
  ##E77F0F，#38B000，#D1B1CB，#2196F3 C11.1 11.0 14.1 14.0
  #vlnplot_data_tmp$new=factor(vlnplot_data_tmp$new,levels=c('C14.1', 'C14.0', 'C11.0', 'C11.1','new'))
  # levels=c("C3","C5","C10","C16","C1","C13","C9","C0","C4","C2","C6","C12","C14","C11","C8","C7","C15")
  # vlnplot_data_tmp$cluster=factor(vlnplot_data_tmp$cluster,levels=levels)# fill=c(new)
  vlnplot_tmp_value=vlnplot_tmp
  # vlnplot_tmp_value=vlnplot_tmp_value[vlnplot_tmp_value$cluster%in%c("C1.0","C1.4","C14.0", "C11.0"),]# 20220824修改
  #######20220825修改C1.2; C1.0;  C1.1; C1.4;    C12; C14.1;C14.0, C11.0, C11.1; C8
  vlnplot_tmp_value=vlnplot_tmp_value[vlnplot_tmp_value$cluster%in%c("C1.2","C1.0","C1.1","C1.4","C12","C14.1","C14.0","C11.0", 
                                                                      "C11.1","C8"),] ## 20220825修改
  tmp_levels=setdiff(levels,c('C1.3','C1.5'))
  vlnplot_tmp$cluster=factor(vlnplot_tmp$cluster,levels=levels,order=T)
  # vlnplot_tmp_value$vlnplot_cluster=factor(vlnplot_tmp_value$vlnplot_cluster,levels=tmp_levels,order=T)
  # vlnplot_tmp_value=vlnplot_tmp
  vlnplot_tmp$vlnplot_value=vlnplot_tmp$value
  vlnplot_tmp$vlnplot_value[!vlnplot_tmp$cluster%in%c("C1.2","C1.0","C1.1","C1.4","C12","C14.1","C14.0","C11.0", "C11.1","C8")]=0
  vlnplot_tmp$new=as.character(vlnplot_tmp$cluster)
  vlnplot_tmp$new[!vlnplot_tmp$cluster%in%c("C1.2","C1.0","C1.1","C1.4","C12","C14.1","C14.0","C11.0", 
                                                 "C11.1","C8")]='new'
  vlnplot_tmp$new=factor(vlnplot_tmp$new,levels=c(tmp_levels,'new'))
  # vlnplot_tmp_value$cluster=as.character(vlnplot_tmp_value$cluster)
  # vlnplot_tmp_value=vlnplot_tmp_value[vlnplot_tmp_value$cluster%in%c("C1.2","C1.0","C1.1","C1.4","C12","C14.1","C14.0","C11.0", 
  #                                                                    "C11.1","C8"),] ## 20220825修改
  # vlnplot_tmp_value$vlnplot_cluster=as.character(vlnplot_tmp_value$cluster)
  # # vlnplot_tmp$vlnplot_cluster[!vlnplot_tmp$vlnplot_cluster%in%tmp_levels]=NA
  # vlnplot_tmp_value$vlnplot_cluster=factor(vlnplot_tmp_value$vlnplot_cluster,levels=tmp_levels)
  vlnplot_tmp$cluster=as.character(vlnplot_tmp$cluster)
  new_data=c()
  for(i in levels){
    tmp=vlnplot_tmp[vlnplot_tmp$cluster==i,]
    new_data=rbind(new_data,tmp)
  }
  new_data_value=c()
  for(i in tmp_levels){
    tmp=vlnplot_tmp_value[vlnplot_tmp_value$cluster==i,]
    new_data_value=rbind(new_data_value,tmp)
  }
  new_data$cluster=factor(new_data$cluster,levels=levels,order=F)
  new_data_value$vlnplot_cluster=factor(new_data_value$vlnplot_cluster,levels=tmp_levels,order=F)
  new_data$vlnplot_cluster='new'
  new_data$vlnplot_cluster[as.character(new_data$cluster)%in%tmp_levels]=as.character(new_data$cluster[as.character(new_data$cluster)%in%tmp_levels])
  new_data$vlnplot_cluster=factor(new_data$vlnplot_cluster,levels=c(tmp_levels,'new'))
  new_data$vlnplot_value[!as.character(new_data$cluster)%in%tmp_levels]=NA
  a <- ggplot() + #ggnewscale::new_scale()+ #stroke = 0.3
    geom_violin(data=new_data,mapping=aes(x=cluster,y=vlnplot_value,fill=vlnplot_cluster),
                scale = "width", trim = TRUE,width=1) +
    geom_jitter(data=new_data, mapping=aes(x=factor(cluster),y=value),
                width = 0.2,size=0.5,stroke=0.3)+theme_bw()+ #ggnewscale::new_scale("fill")+
     #ggnewscale::new_scale('fill')+#adjust = 1,
    # geom_jitter(data=new_data, mapping=aes(x=cluster,y=value),width = 0.2,size=0.5,stroke=0.3)+theme_bw()+ #ggnewscale::new_scale("fill")+
    scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(rows = vars(newgene), scales = "free", switch = "y") +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          strip.background.y = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          axis.title.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          axis.title.y.right=element_text(size=34,face='bold',colour='black'),
          axis.text.y.left = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold",size=25),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=25,colour='black',face='bold'),
          # axis.text.y=element_text(size=14,face='bold',colour='black'),
          strip.text.y.left = element_text(angle = 0))+ #hjust=1,忽略了这个调整
    ggtitle("") + xlab("cluster") + ylab("Expression Level")+
    scale_fill_manual(values=cluster9_cor_new[tmp_levels])
  ###,hjust=1
  a1 <- a+theme(axis.text.y = element_text(vjust = 0,size=25,colour='black',face='bold'))
  ####30个基因长度是27，
  height=dim(gene1)[1]*30/27
  ####3带基因名字的是width=20
  ggsave(paste(name,'_lnplot_',gsub('\\-','',Sys.Date()),'.png',sep=""),a1,width=16,height=height,limitsize = F)
  ggsave(paste(name,'_lnplot_',gsub('\\-','',Sys.Date()),'.pdf',sep=""),a1,width=16,height=height,limitsize = F)
}
#####################################做基因的umap图20220824
# UMAP做图： C1.0;  C1.1; C1.2;  C1.3; C1.4; C1.5;   C12; C14.1;C14.0, C11.0, C11.1; C8（每排五个基因）
# UMAP做图： C1.0;  C1.4; C14.0, C11.0（每排五个基因）
namelist=list.files(getwd(),pattern='*txt$')
cluster1=c("C1.0",  "C1.1", "C1.2",  "C1.3", "C1.4", "C1.5",
           "C12", "C14.1","C14.0", "C11.0", "C11.1", "C8")
# cluster1=c('C1.0','C1.4','C14.0','C11.0')
# Idents(expr)=expr$new 
path='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220824_R1_vlnplot'
savePath=path
for(s in 1:length(namelist)){
  name=namelist[s]
  gene_express=read.delim(file.path(path,namelist[s]),sep="\t",as.is=T,header=T,check.names=F)
  if(name=='vlnplot20220824_sheet1V1.txt'){
    gene_express$new=paste(gene_express[,1],gene_express[,2],sep="|")
    gene_express=gene_express[,c(1,3)]
    # colnames(gene_express)=c('gene','geneID')
  }
  if(ncol(gene_express)>2){
    gene_express=gene_express[,c(1,2)]
  }
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
  # if(s %in% c(1,2,3)){
  width=20
  ncol=5
  # }else if(s==4){
  #   width=15
  #   ncol=4
  # }
  library(patchwork)
  p=list()
  p1=FeaturePlot(expr1, features = gene_express[,1], cols = c("grey",  "red"), reduction = "umap",
                 combine = FALSE,order=TRUE,min.cutoff = 'q5',max.cutoff = 'q95')
  for(i in 1:length(p1)){
    p1[[i]]$labels$title = gene_express[i,2]
  }
  # p2=FeaturePlot(expr1, features = gene_express[6,1], cols = c("grey",  "red"), reduction = "umap",
  #                combine = FALSE,order=TRUE)
  # p2[[1]]$labels$title=gene_express[6,2]
  # p3=FeaturePlot(expr1, features = gene_express[7:nrow(gene_express),1], cols = c("grey",  "red"), reduction = "umap",
  #                combine = FALSE,order=TRUE,min.cutoff = 'q5',max.cutoff = 'q95')
  # for(i in 1:length(p3)){
  #   p3[[i]]$labels$title = gene_express[i+6,2]
  # }
  # p=c(p1,p2,p3)
  pdf(file.path(savePath, paste(name, "_marker_umap.pdf", sep = "")), width = width, height =height ,onefile = FALSE)
  print(CombinePlots(p1,ncol=5))
  dev.off()
  png(file.path(savePath, paste(name, "_marker_umap.png", sep = "")), width = width*100, height =height*100 )
  print(CombinePlots(p1,ncol=5))
  dev.off()
}

###################################拟时序图和折线图20220824
# "拟时序1的折线图横坐标顺序：C1.5；C1.2; C1.0； C1.3；C1.1; C1.4				（做拟时序图和折线图，具体参考右图）	"
# 拟时序2的横坐标顺序：C12; C14.1;C14.0, C11.0, C11.1; C8
# 
# 拟时序3的横坐标顺序： C1.0;  C1.2;  C1.5; C14.0, C14.1;,C12; 
# 
# 拟时序4的横坐标顺序：C1.0;  C1.3; C1.4;  C11.0, C11.1; C8
#########拟时序1得折线图顺序 20220824
library(reshape2)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(VisiumSpatial)
#C1拟时序存放路径
#/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1'
load(file.path(indir,'C1_gene_curveline/20220127_curveline.rdata'))
namelist=list.files(getwd(),pattern='*txt$')
s=1
gene_express=c()
for(s in 1:3){
  tmp_gene_express=read.delim(file.path(getwd(),namelist[s]),sep="\t",as.is=T,header=T,check.names=F)
  if(namelist[s]=='vlnplot20220824_sheet1V1.txt'){
    tmp_gene_express$new=paste(tmp_gene_express[,1],tmp_gene_express[,2],sep="|")
    tmp_gene_express=tmp_gene_express[,c(1,3)]
    # colnames(gene_express)=c('gene','geneID')
  }
  if(ncol(tmp_gene_express)>2){
    tmp_gene_express=tmp_gene_express[,c(1,2)]
  }
  tmp_gene_express=tmp_gene_express[!duplicated(tmp_gene_express[,1]),]
  colnames(tmp_gene_express)=c('gene','geneID')
  rownames(tmp_gene_express)=tmp_gene_express[,1]
  gene_express=rbind(gene_express,tmp_gene_express)
}
tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%gene_express[,1],]
tmp_figure$clusters=factor(tmp_figure$clusters,levels=c("C1.5","C1.2", "C1.0", "C1.3","C1.1", "C1.4"))
########再画基因折线图得figure_list，然后把折线图和拟时序图都放在一起，一张图中，上下得排布
####定义C1基因得折线图
levels=c("C1.5","C1.2", "C1.0", "C1.3","C1.1", "C1.4")
geneline<-function(expr_plotdata,gene_express,cluster_levels){
  expr_figure_line=list()
  for(i in 1:nrow(gene_express)){
    #tmp_figure=expr_plotdata[expr_plotdata$newgenecluster==genecluster,]
    # tmp_genes=gene_info[gene_info[,'genecluster']==genecluster,1]
    # tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%tmp_genes,]
    tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%gene_express[i,1],]
    tmp_figure$clusters=factor(tmp_figure$clusters,levels=cluster_levels)
    tmp_figure$subcluster=tmp_figure$clusters
    # tmp_figure_outs=c()
    # for(i in as.character(levels(tmp_figure$subcluster))){
    #   tmp_figure_outs_1=tmp_figure[tmp_figure$subcluster==i,]
    #   # tmp_figure_outs_1$expectation=limitdata(tmp_figure_outs_1$expectation)
    #   tmp_figure_outs=rbind(tmp_figure_outs,tmp_figure_outs_1)
    # }
    # tmp_figure_outs=tmp_figure
    tmp_max_value=aggregate(tmp_figure$expectation,by=list(tmp_figure$clusters),median)
    colnames(tmp_max_value)=c('subcluster','value')
    # tmp_max_value$subcluster=factor(as.character(tmp_max_value$subcluster),levels=c("C12","C14","C11.0","C11.1","C8"))
    print(tmp_max_value)
    tmp_max_value[tmp_max_value$value > 2,'value']=2
    tmp_max_value[tmp_max_value$value < (-2),'value']=-2
    
    # p=ggplot(tmp_figure_outs,aes(x=subcluster,y=expectation,color=subcluster))+geom_line()+geom_point()
    # p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value))+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+labs(title = paste('genecluster',genecluster,sep=""))+
    # theme(plot.title = element_text(hjust = 0.5))
    # p_line=ggplot(tmp_max_value,aes(x=subcluster,y=value))+geom_line(size = 3,group=1,col=color_list[[genecluster]])+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+labs(title = paste('genecluster',genecluster,sep=""))+
    #   theme(plot.title = element_text(hjust = 0.5))
    p=ggplot(tmp_max_value,aes(x=subcluster,y=value))
    p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value),size=2.5,col='black')+
      theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
      labs(title = gene_express[i,2],x="",y='RelativeExp')+geom_point(size = 3,col='black')+
      theme(plot.title = element_text(hjust = 0.5,size=38,face='bold'),
            text = element_text(size=26,face='bold',colour="black"),
            axis.text.x = element_text(size=26,angle=90,hjust=1,face='bold',colour="black"),
            axis.text.y = element_text(size=26,face='bold',colour="black"))+scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2),
                                                                                               labels=c('-2.0','-1.0','0','1.0','2.0'))
    expr_figure_line[[gene_express[i,1]]]=p
    # expr_figure_line_point_alpha[[paste('genecluster',genecluster,sep="")]]=p_line
  }
  return(expr_figure_line)
}
c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
library(monocle)
HSMM_myo1=HSMM_myo
# C1_genes=read.delim(file.path('C12_8_11_14','genes.txt'),sep="\t",as.is=T,header=T,check.names=F)
# C1_genes=read.delim(file.path('C1_11_8_linecurve.txt'),sep="\t",as.is=T,header=T,check.names=F)
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')

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

x2$clusters=factor(x2$clusters,levels=c('C1.5','C1.2','C1.0','C1.3','C1.1','C1.4'))
HSMM_myo1@phenoData@data=x2
c1_trajectorylist=trajectory_function(HSMM_myo1=HSMM_myo1,gene_express)


#######整合结果
wrap_plots_combin(linelist=c1_linelist,trajectorylist=c1_trajectorylist,cluster_name='C1')
  
wrap_plots_combin<-function(linelist,trajectorylist,cluster_name){
  num=length(linelist)
  fi=trajectorylist[[1]]
  fi_pdf=trajectorylist[[2]]
  id=floor(num/5)
  other=num-5*id
  library(patchwork)
  new_figure_list=list()
  for(i in 1:id){
  new_figure_list=c(new_figure_list,fi[((i-1)*5+1):(5*i)],linelist[((i-1)*5+1):(5*i)])
  }
  figure_list=wrap_plots(new_figure_list,ncol=5)
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes1_',5*id,'.png',sep=""),figure_list,
         width=35,height=5*id*2,
         dpi=300,limitsize = F)######## 4行5列对应宽35，高20
  if(5*id < num){
    new_figure_list=list()
    new_figure_list=c(fi[(5*id):num],linelist[(5*id):num])
    figure_list=wrap_plots(new_figure_list,ncol=5)
    ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes',5*id,'_',num,'.png',sep=""),figure_list,
           width=35,height=20,dpi=300,limitsize = F)
  }
  new_figure_list=list()
  for(i in 1:id){
    new_figure_list=c(new_figure_list,fi[((i-1)*5+1):(5*i)],linelist[((i-1)*5+1):(5*i)])
  }
  figure_list=wrap_plots(new_figure_list,ncol=5)
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes1_',5*id,'.pdf',sep=""),figure_list,
         width=35,height=5*id*2,
         dpi=300,limitsize = F)######## 4行5列对应宽35，高20
  if(5*id<num){
  new_figure_list=list()
  new_figure_list=c(fi[(5*id):num],linelist[(5*id):num])
  # new_figure_list=c(fi[11:12],expr_line[11:12])
  figure_list=wrap_plots(new_figure_list,ncol=5)
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes',5*id,'_',num,'.pdf',sep=""),figure_list,
         width=35,height=20,dpi=300,limitsize = F)
  }
} 


##############C1基因得拟时序

trajectory_function<-function(HSMM_myo1,gene_express){
  fi=list()
  for(i in 1:nrow(gene_express)){
    p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = gene_express[i,1], 
                            cell_size = 1.2) + 
      theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
      # labs(title = C1_genes[i,2],x="",y='RelativeExp')+
      theme(plot.title = element_text(hjust = 0.5,size=38,face='bold'),
            text = element_text(size=26,face='bold',colour="black"),
            axis.text.x = element_text(size=26,hjust=1,face='bold',colour="black"),
            axis.text.y = element_text(size=26,face='bold',colour="black"))+
      # theme(axis.title = element_text(face = "bold", size = 26, colour = "black"), legend.position = "right", 
      #       legend.text = element_text(face = "bold", size = 26, 
      #                                  colour = "black"), legend.title = element_text(face = "bold", 
      #         size = 26, colour = "black"), axis.text.x = element_text(face = "bold", colour='black',
      #         size = 26), axis.text.y = element_text(face = "bold", size = 26,colour='black'), 
      #       title = element_text(size = 38,face='bold',family='Arial'), 
      #       strip.text = element_text(face = "bold", size = 22),legend.key.height = unit(1.3,"cm"))+
      labs(title=gene_express[i,2],family = "Arial",face='bold',col='black')+
      theme(plot.title = element_text(hjust = 0.5,vjust = 0.3)) + #也就加上这一行
      guides(fill=guide_legend(label.hjust=1.5,label.vjust=1.5))+
      guides(shape = guide_legend(override.aes = list(size = 3)))+
      #low="#CDC9C9",mid='#FFFF00', high = "#FF0000"
      scale_colour_gradientn(colours=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ),
                             name=NULL)
    # scale_fill_gradientn(colors=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ))
    # scale_color_viridis_c()
    fi[[gene_express[i,1]]]=p1
  }
  fi_pdf=list()
  for(i in 1:nrow(gene_express)){
    p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = gene_express[i,1], 
                            cell_size = 1.2) + 
      theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
      # labs(title = C1_genes[i,2],x="",y='RelativeExp')+
      theme(plot.title = element_text(hjust = 0.5,size=38,face='bold'),
            text = element_text(size=26,face='bold',colour="black"),
            axis.text.x = element_text(size=26,hjust=1,face='bold',colour="black"),
            axis.text.y = element_text(size=26,face='bold',colour="black"))+
      # theme(axis.title = element_text(face = "bold", size = 26, colour = "black"), legend.position = "right",
      #       legend.text = element_text(face = "bold", size = 26,
      #                                  colour = "black"), legend.title = element_text(face = "bold",
      #   size = 26, colour = "black"), axis.text.x = element_text(face = "bold", colour='black',
      #     size = 26), axis.text.y = element_text(face = "bold", size = 26,colour='black'),
      #       title = element_text(size = 38,face='bold'),
      #       strip.text = element_text(face = "bold", size = 22),legend.key.height = unit(1.3,"cm"))+
      labs(title=gene_express[i,2],face='bold',col='black')+
      # theme(plot.title = element_text(hjust = 0.5,vjust = 0.3)) + #也就加上这一行
      guides(fill=guide_legend(label.hjust=1.5,label.vjust=1.5))+
      guides(shape = guide_legend(override.aes = list(size = 3)))+
      #low="#CDC9C9",mid='#FFFF00', high = "#FF0000"
      scale_colour_gradientn(colours=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ),
                             name=NULL)
    
    # scale_fill_gradientn(colors=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ))
    # scale_color_viridis_c()
    fi_pdf[[gene_express[i,1]]]=p1
  }
  return(list(fi,fi_pdf))
}

########################## 拟时序2的横坐标顺序：C12; C14.1;C14.0, C11.0, C11.1; C8 20220824
#最新路径/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220311_genemonocle_C1_C12_8_11_14/
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220311_genemonocle_C1_C12_8_11_14/'
expr_plotdata=readRDS(file.path(indir,'expr_plotdata_copy.RDS'))
HSMM_myo1=readRDS(file.path(indir,'HSMM_myo1.RDS'))
levels=c('C12','C14.1','C14.0','C11.0','C11.1','C8')
c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
# HSMM_myo=HSMM_myo1
c1_traje=get_hsmm_myo_info(HSMM_myo=HSMM_myo1,gene_express=gene_express,levels=levels)
wrap_plots_combin(linelist=c1_linelist,trajectorylist=c1_traje,cluster_name='C12_14_11_8')
########################################拟时序3 拟时序3的横坐标顺序： C1.0;  C1.2;  C1.5; C14.0, C14.1;,C12; 
##最新路径/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_14_12_curveline
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_14_12_curveline'
expr_plotdata=readRDS(file.path(indir,'c1_14_12_cds_exprs.RDS'))
load(file.path(dirname(indir),'1_14_12_HSMM_myo.Rdata'))
levels=c("C1.0",  "C1.2",  "C1.5", "C14.0", "C14.1","C12")
c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
# HSMM_myo=HSMM_myo1
c1_traje=get_hsmm_myo_info(HSMM_myo=HSMM_myo,gene_express=gene_express,levels=levels)
wrap_plots_combin(linelist=c1_linelist,trajectorylist=c1_traje,cluster_name='C1_14_12')
###################################拟时序4的横坐标顺序：C1.0;  C1.3; C1.4;  C11.0, C11.1; C8

get_hsmm_myo_info<-function(HSMM_myo,gene_express,levels){
  HSMM_myo1=HSMM_myo
  # C1_genes=read.delim(file.path('C12_8_11_14','genes.txt'),sep="\t",as.is=T,header=T,check.names=F)
  # C1_genes=read.delim(file.path('C1_11_8_linecurve.txt'),sep="\t",as.is=T,header=T,check.names=F)
  expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')
  
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
  return(c1_trajectorylist)
}

