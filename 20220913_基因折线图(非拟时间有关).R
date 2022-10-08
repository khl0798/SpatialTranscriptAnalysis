########################20220909基因折线图非拟时间有关
##首先是C11.1， C8 与C 14.1， C12 做差异分析
library('VisiumSpatial')
library('Seurat')
library("ggplot2")
library('dplyr')
# expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')

# total_expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.data.RDS')

# expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')

# expdata=t(as.matrix(expr[['Spatial']]@data))
species=NULL
# Idents(expr)=expr$total_new_info
Idents(expr)=expr$clusters
savePath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220513_linecurve_4groups/20220516'
namelist=list.files(pattern='.txt')##折线图横坐标的顺序:C1.0， C1.4,  C14.0, C11.0
Idents(expr)=expr$cluster20220822_new
expr1=subset(expr,idents=c("C1.1","C1.2","C1.5","C1.3","C1.0","C1.4","C14.1","C8","C11.0",
                           "C12","C11.1","C14.0"))
data=t(as.matrix(expr1[['Spatial']]@data))
data=log(data+1)
data[data> 3]=3
data[data< (-3)]=-3
data=data.frame(data,cluster=paste(as.character(Idents(expr1)),sep=""))
data_mean=aggregate(data[,1:(ncol(data)-1)],by=list(data$cluster),mean)
#######bar图需要计算
# bar_compute<-function(data,upper=TRUE){
#   # data=log(data+1)
#   if(upper){
#     upper=mean(data)+sd(data)/2
#     return(upper)
#   }else{
#     lower=mean(data)-sd(data)/2
#     return(lower)
#   }
# }
# tmp_upper=aggregate(data[,1:(ncol(data)-1)],by=list(data$cluster),function(t)bar_compute(t,upper=TRUE))
# tmp_lower=aggregate(data[,1:(ncol(data)-1)],by=list(data$cluster),function(t)bar_compute(t,upper=FALSE))
trans<-function(tmp_data){
  rownames(tmp_data)=tmp_data[,1]
  tmp_data=t(tmp_data[,-1])
  # tmp_data=log(tmp_data+1)
  return(tmp_data)
}
# tmp_upper=trans(tmp_upper)
# tmp_lower=trans(tmp_lower)
data_mean=trans(data_mean)
namelist=list.files(getwd(),pattern='*txt$')
# levels=c('C1.0','C1.3','C1.4','C11.0','C11.1','C8')
levels_list=list()
levels_list[['C1_11_8']]=c('C1.0','C1.3','C1.4','C11.0','C11.1','C8')
levels_list[['C1_14_12']]=c('C1.0','C1.2','C1.5','C14.0','C14.1','C12')
levels_list[['C12_14_11_8']]=c('C12','C14.1','C14.0','C11.0','C11.1','C8')
for(levels_name in names(levels_list)){
  levels=levels_list[[levels_name]]
  for(i in 1:length(namelist)){
    print(i)
    tmp_gene_express=read.delim(file.path(getwd(),namelist[i]),sep="\t",as.is=T,header=T,check.names=F)
    tmp_gene_express=tmp_gene_express[!duplicated(tmp_gene_express[,1]),]
    # print(paste('gene number is',dim(tmp_gene_express)[1],'and total gene in ',sum(tmp_gene_express[,1]%in%expr_plotdata$gene_short_name)))
    colnames(tmp_gene_express)=c('gene','geneID')
    rownames(tmp_gene_express)=tmp_gene_express[,1]
    gene_express=tmp_gene_express
    
    anno=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/anno_tf_go_kegg_ko.RDS')
    
    library(reshape2)
    x1_data=melt(data_mean)
    colnames(x1_data)=c('genecluster','subcluster','value')
    x1_data[,1:2]=apply(x1_data[,1:2],2,as.character)
    x1_data_max_min=x1_data[x1_data$genecluster%in%gene_express[,1],]
    x1_data_max_min=x1_data_max_min[x1_data_max_min$subcluster%in%levels,]
    x1_data_max_min$genecluster=factor(x1_data_max_min$genecluster,levels=gene_express[,1])
    x1_data_max_min$subcluster=factor(x1_data_max_min$subcluster,levels=levels)
    
    library(dplyr)
   
    genes_split=list()
    genes_split[[namelist[i]]]=gene_express[,1]
    for(name in names(genes_split)){
      print(name)
      expr_figure_line=list()
      # expr_figure_line_bar=list()
      for(s in genes_split[[name]]){
        tmp_max_value=data.frame(x1_data_max_min)%>%filter(genecluster==s)
        # print(tmp_max_value)
        p=ggplot(tmp_max_value,aes(x=subcluster,y=value))
        # p1=ggplot(tmp_max_value,aes(x=subcluster,y=value))
        p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value),size=2.5)+
          theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
          labs(title = gsub(':',"|",gene_express[s,2]),x="",y='Log Exp')+geom_point(size = 3)+
          theme(plot.title = element_text(hjust = 0.5,size=30,face='bold'),
                text = element_text(size=22,face='bold',colour="black"),
                axis.text.x = element_text(size=22,angle=90,hjust=1,face='bold',colour="black"),
                axis.text.y = element_text(size=22,face='bold',colour="black"))
        # scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2),labels=c('-2.0','-1.0','0','1.0','2.0'))
        expr_figure_line[[paste('genecluster',s,sep="")]]=p
      }
      library(patchwork)
      expr_line=wrap_plots(expr_figure_line,ncol=5)
      id=floor(dim(gene_express)[1]/5)
      if(id*5<dim(gene_express)[1]){
        id=id+1
      }
      cluster_name=paste(paste(levels,collapse="_"),namelist[i],sep="_")
      ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'VnobarV1.png',sep=""),expr_line,width=35+3,height=5*id)
      ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'VnobarV1.pdf',sep=""),expr_line,width=35,height=5*id,limitsize = F)
    }
  }
}