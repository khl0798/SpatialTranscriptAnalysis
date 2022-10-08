library(reshape2)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(VisiumSpatial)
library(monocle)
########20220826 画拟时序1和4的折线图和拟时间图等
# 拟时序1的折线图横坐标顺序：C1.5；C1.2; C1.0； C1.3；C1.1; C1.4	

#############结果存放路径/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220826_R1_C1andC12_14_11_8_genecurve_pseudotime
savePath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220826_R1_C1andC12_14_11_8_genecurve_pseudotime'


#C1拟时序存放路径
#/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1'
load(file.path(indir,'C1_gene_curveline/20220127_curveline.rdata'))
namelist=list.files(getwd(),pattern='*txt$')
s=1
gene_express=c()
for(s in seq(1,length(namelist))){
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
    tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%gene_express[i,1],]
    tmp_figure$clusters=factor(tmp_figure$clusters,levels=cluster_levels)
    tmp_figure$subcluster=tmp_figure$clusters
    tmp_max_value=aggregate(tmp_figure$expectation,by=list(tmp_figure$clusters),median)
    colnames(tmp_max_value)=c('subcluster','value')
    print(tmp_max_value)
    tmp_max_value[tmp_max_value$value > 2,'value']=2
    tmp_max_value[tmp_max_value$value < (-2),'value']=-2
    
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
  }
  return(expr_figure_line)
}
c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
library(monocle)
HSMM_myo1=HSMM_myo
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')

x2=pData(HSMM_myo1)
new_expr=as_matrix(expr[['Spatial']]@data)
new_expr=new_expr[,rownames(x2)]

for(i in 1:nrow(gene_express)){
  new_name=gene_express[i,2]
  x2=cbind(x2,new_expr[gene_express[i,1],])
}
colnames(x2)[6:ncol(x2)]=gene_express[1:nrow(gene_express),1]

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
    width_nn=num-(5*id)
    new_figure_list=list()
    new_figure_list=c(fi[(5*id+1):num],linelist[(5*id+1):num])
    figure_list=wrap_plots(new_figure_list,ncol=4,nrow=2)
    ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes',5*id+1,'_',num,'.png',sep=""),figure_list,
           width=7*width_nn-0.8,height=10,dpi=300,limitsize = F)
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
    new_figure_list=c(fi[(5*id+1):num],linelist[(5*id+1):num])
    figure_list=wrap_plots(new_figure_list,ncol=width_nn)
    ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes',5*id+1,'_',num,'.pdf',sep=""),figure_list,
           width=7*width_nn-0.8,height=10,dpi=300,limitsize = F)
  }
} 


##############C1基因得拟时序

trajectory_function<-function(HSMM_myo1,gene_express){
  fi=list()
  for(i in 1:nrow(gene_express)){
    p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = gene_express[i,1], 
                            cell_size = 1.2) + 
      theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
      theme(plot.title = element_text(hjust = 0.5,size=38,face='bold'),
            text = element_text(size=26,face='bold',colour="black"),
            axis.text.x = element_text(size=26,face='bold',colour="black"),
            axis.text.y = element_text(size=26,face='bold',colour="black"))+
      labs(title=gene_express[i,2],family = "Arial",face='bold',col='black')+
      theme(plot.title = element_text(hjust = 0.5,vjust = 0.3)) + #也就加上这一行
      guides(fill=guide_legend(label.hjust=1.5,label.vjust=1.5))+
      guides(shape = guide_legend(override.aes = list(size = 3)))+
      scale_colour_gradientn(colours=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ),
                             name=NULL)
    fi[[gene_express[i,1]]]=p1
  }
  fi_pdf=list()
  for(i in 1:nrow(gene_express)){
    p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = gene_express[i,1], 
                            cell_size = 1.2) + 
      theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
      theme(plot.title = element_text(hjust = 0.5,size=38,face='bold'),
            text = element_text(size=26,face='bold',colour="black"),
            axis.text.x = element_text(size=26,face='bold',colour="black"),  #hjust=1,
            axis.text.y = element_text(size=26,face='bold',colour="black"))+
      labs(title=gene_express[i,2],face='bold',col='black')+
      guides(fill=guide_legend(label.hjust=1.5,label.vjust=1.5))+
      guides(shape = guide_legend(override.aes = list(size = 3)))+
      scale_colour_gradientn(colours=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ),
                             name=NULL)
    fi_pdf[[gene_express[i,1]]]=p1
  }
  return(list(fi,fi_pdf))
}
######################20220826,C1_11_8拟时序画图 # 拟时序4的横坐标顺序：C1.0;  C1.4;  C11.0, C11.1; C8	
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_11_8_curveline'
load(file.path(indir,'20220201.rdata'))
namelist=list.files(getwd(),pattern='*txt$')
s=1
levels=c('C1.0','C1.4','C11.0','C11.1','C8')
gene_express=c()
for(s in seq(1,length(namelist))){
  tmp_gene_express=read.delim(file.path(getwd(),namelist[s]),sep="\t",as.is=T,header=T,check.names=F)
  if(namelist[s]=='vlnplot20220824_sheet1V1.txt'){
    tmp_gene_express$new=paste(tmp_gene_express[,1],tmp_gene_express[,2],sep="|")
    tmp_gene_express=tmp_gene_express[,c(1,3)]
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
tmp_figure$clusters=factor(tmp_figure$clusters,levels=levels)
c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
library(monocle)
HSMM_myo1=HSMM_myo
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')

x2=pData(HSMM_myo1)
new_expr=as_matrix(expr[['Spatial']]@data)
new_expr=new_expr[,rownames(x2)]

for(i in 1:nrow(gene_express)){
  new_name=gene_express[i,2]
  x2=cbind(x2,new_expr[gene_express[i,1],])
}
colnames(x2)[6:ncol(x2)]=gene_express[1:nrow(gene_express),1]
x2$clusters=factor(x2$clusters,levels=levels)
HSMM_myo1@phenoData@data=x2
c1_trajectorylist=trajectory_function(HSMM_myo1=HSMM_myo1,gene_express)


#######整合结果
wrap_plots_combin(linelist=c1_linelist,trajectorylist=c1_trajectorylist,cluster_name='C1_11_8')
####################20220826,重新更新切片图等
###存放路径是/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220826_R1_C1andC12_14_11_8_genecurve_pseudotime/
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')
Idents(expr)=expr$default
c1_new=read.csv('C1_6clusters-20220812.csv',sep=",",as.is=T,header=T,check.names=F)
new_idents=as.character(Idents(expr))
names(new_idents)=colnames(expr)
c1.3cells=c1_new%>%filter(C1_6clusters=='C1.3 PC-b')
c1.3cells=c1.3cells[,1]
new_idents[c1.3cells]='1'
Idents(expr)=factor(new_idents[colnames(expr)],levels=0:16)
library(scales)
color=hue_pal()(17)
for(pt.size.factor in c(1.2,1.3)){
  print(pt.size.factor)
  p1=SpatialDimPlot(expr,label = F,image.alpha = 1,ncol = 2,cols=color,crop=F,pt.size.factor = pt.size.factor)
  ggsave(paste('17clusters_img_20220826pt',pt.size.factor,'.png',sep=""),p1,width = 8 + 2 * 2, height = 6)
  ggsave(paste('17clusters_img_20220826pt',pt.size.factor,'.pdf',sep=""),p1,width = 8 + 2 * 2, height = 6)
}

############20220907，切片图
  for(slice in unique(expr$orig.ident)){
    new_labels=as.character(sort(unique(expr@meta.data[expr$orig.ident==slice,idname])))
    print(c(slice,new_labels))
    new_color=color[new_labels]
    names(new_color)=new_labels
    p1=SpatialDimPlot(expr,crop=FALSE,label=FALSE,images=slice,
                      image.alpha = 0.7,
                      pt.size.factor = 1.2,ncol = 1,repel=TRUE)+
      theme(legend.text = element_text(size=8,face='bold'),
            legend.title = element_text(face = "bold", size = 12, colour = "black"),
            title = element_text(size = 12,colour='black',face='bold'),
            plot.title = element_text(hjust = 0.5))+
      guides(fill= guide_legend(override.aes = list(size = 10)))+labs(fill='Cluster',title=new_names[slice])+
      scale_fill_manual(values=new_color,labels=new_labels)
    figure_list[[slice]]=p1
  }
  p1=CombinePlots(figure_list,ncol=2)
  if(length(unique(expr$orig.ident))==3){
    height=6
  }else{
    height=12
  }
  ggsave(paste(gsub('-','',Sys.Date()),name,'img.png',sep="."),p1,dpi=500,width=17,height=12)
  ggsave(paste(gsub('-','',Sys.Date()),name,'img.pdf',sep="."),p1,dpi=500,width=17,height=12)
  new_names=c('N1-1','N2-1')
  names(new_names)=c('1-N1-1','2-N2-1')
  figure_list1=list()
  for(slice in unique(expr@meta.data$orig.ident)){
    subsetcells=colnames(expr)[(expr$region==slice)]
    expr_sub=subset(expr,cells=subsetcells)
    new_labels=as.character(sort(unique(expr@meta.data[expr$orig.ident==slice,idname])))
    print(c(slice,new_labels))
    new_color=color[new_labels]
    names(new_color)=new_labels
    expr_sub=subset(expr_sub,idents=setdiff(new_labels,'null'))
    p1=SpatialDimPlot(expr_sub,crop=FALSE,label=FALSE,images=paste('X',gsub('-','.',slice),sep=""),
                      image.alpha = 0.7,
                      pt.size.factor = 1.2,ncol = 1,repel=TRUE)+
      theme(legend.text = element_text(size=8,face='bold'),
            legend.title = element_text(face = "bold", size = 12, colour = "black"),
            title = element_text(size = 12,colour='black',face='bold'),
            plot.title = element_text(hjust = 0.5))+
      guides(fill= guide_legend(override.aes = list(size = 10)))+labs(fill='Cluster',title=new_names[slice])+
      scale_fill_manual(values=new_color[setdiff(new_labels,'null')],labels=setdiff(new_labels,'null'))
    figure_list1[[slice]]=p1
  }
  p1=CombinePlots(figure_list1,ncol=2)
  if(length(unique(expr$orig.ident))==3){
    height=6
  }else{
    height=12
  }
  ggsave(paste(gsub('-','',Sys.Date()),name,'imgV1.png',sep="."),p1,dpi=500,width=17,height=12)
  ggsave(paste(gsub('-','',Sys.Date()),name,'imgV1.pdf',sep="."),p1,dpi=500,width=17,height=12)
}