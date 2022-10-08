#####基因的折线图
## 王亮，我突然想到，如果用C11.1， C8 与C 14.1， C12 里有显著差异表达的基因，
## 再分别与C11.0，C14.0的marker gene取交集，是不是就可以得到在 C11.0， C11.1里高表达的基因，
## 同时满足这些基因肯定是在C14.0， C14.1 低表达？

##首先是C11.1， C8 与C 14.1， C12 做差异分析
library('VisiumSpatial')
library('Seurat')
library("ggplot2")
library('dplyr')
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')

# total_expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.data.RDS')

# expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')

# expdata=t(as.matrix(expr[['Spatial']]@data))
species=NULL
Idents(expr)=expr$total_new_info
savePath='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220513_linecurve_4groups/20220516'
namelist=list.files(pattern='.txt')##折线图横坐标的顺序:C1.0， C1.4,  C14.0, C11.0
for(i in 1:length(namelist)){
  expr1=subset(expr,idents=c('C1.0','C1.4','C14.0','C11.0'))
  data=t(as.matrix(expr1[['Spatial']]@data))
  markergene=read.delim(namelist[i],sep="\t",as.is=T,header=T,check.names=F)
  markergene=markergene[!duplicated(markergene[,1]),]
  rownames(markergene)=markergene[,1]
  data=data[,markergene[,1]]
  data=log(data+1)
  data[data> 3]=3
  data[data< (-3)]=-3
  data=data.frame(data,cluster=paste(as.character(Idents(expr1)),sep=""))
  data_mean=aggregate(data[,1:(ncol(data)-1)],by=list(data$cluster),mean)
  bar_compute<-function(data,upper=TRUE){
    # data=log(data+1)
    if(upper){
      upper=mean(data)+sd(data)/2
      return(upper)
    }else{
      lower=mean(data)-sd(data)/2
      return(lower)
    }
  }
  tmp_upper=aggregate(data[,1:(ncol(data)-1)],by=list(data$cluster),function(t)bar_compute(t,upper=TRUE))
  tmp_lower=aggregate(data[,1:(ncol(data)-1)],by=list(data$cluster),function(t)bar_compute(t,upper=FALSE))
  trans<-function(tmp_data){
    rownames(tmp_data)=tmp_data[,1]
    tmp_data=t(tmp_data[,-1])
    # tmp_data=log(tmp_data+1)
    return(tmp_data)
  }
  tmp_upper=trans(tmp_upper)
  tmp_lower=trans(tmp_lower)
  data_mean=trans(data_mean)
  
  anno=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/anno_tf_go_kegg_ko.RDS')
  
  library(reshape2)
  x1_data=melt(data_mean)
  x1_max=melt(tmp_upper)
  x1_min=melt(tmp_lower)
  colnames(x1_data)=c('genecluster','subcluster','value')
  colnames(x1_max)=c('genecluster','subcluster','max')
  colnames(x1_min)=c('genecluster','subcluster','min')
  x1_data[,1:2]=apply(x1_data[,1:2],2,as.character)
  x1_max[,1:2]=apply(x1_max[,1:2],2,as.character)
  x1_min[,1:2]=apply(x1_min[,1:2],2,as.character)
  
  cluster1=c('C1.0','C1.4','C14.0','C11.0')
  x1_data_max=merge(x1_data,x1_max,by=c('genecluster','subcluster'))
  x1_data_max_min=merge(x1_data_max,x1_min,by=c('genecluster','subcluster'))
  
  x1_data_max_min$subcluster=factor(x1_data_max_min$subcluster,levels=paste(cluster1,sep=""))
  # x1_data$genecluster=factor(x1_data$subcluster,levels=c('C12','C14.1','C14.0','C11.0','C11.1','C8'))
  
  # x1_data$value=log(x1_data$value)
  
  library(dplyr)
  # genes_split=list()
  # genes_split[['final']]=markergene[,1]
  ###col=color_list[[genecluster]]
  # num=floor(nrow(markergene)/100)
  # genes_split=split(used_genes[1:(num*100)],f = rep(1:num,each=100))
  # genes_split[[as.character(num+1)]]=used_genes[((num*100)+1):length(used_genes)]
  genes_split=list()
  genes_split[[namelist[i]]]=markergene[,1]
  for(name in names(genes_split)){
    print(name)
    expr_figure_line=list()
    expr_figure_line_bar=list()
    for(s in genes_split[[name]]){
      tmp_max_value=data.frame(x1_data_max_min)%>%filter(genecluster==s)
      # print(tmp_max_value)
      p=ggplot(tmp_max_value,aes(x=subcluster,y=value))
      p1=ggplot(tmp_max_value,aes(x=subcluster,y=value))
      p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value),size=2.5)+
        theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
        labs(title = gsub(':',"|",markergene[s,2]),x="",y='Log Exp')+geom_point(size = 3)+
        theme(plot.title = element_text(hjust = 0.5,size=25,face='bold'),
              text = element_text(size=22,face='bold',colour="black"),
              axis.text.x = element_text(size=22,angle=90,hjust=1,face='bold',colour="black"),
              axis.text.y = element_text(size=22,face='bold',colour="black"))
      p1=p1+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value),size=2.5)+
        theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
        labs(title = gsub(':',"|",markergene[s,2]),x="",y='Log Exp')+geom_point(size = 3)+
        theme(plot.title = element_text(hjust = 0.5,size=25,face='bold'),
              text = element_text(size=22,face='bold',colour="black"),
              axis.text.x = element_text(size=22,angle=90,hjust=1,face='bold',colour="black"),
              axis.text.y = element_text(size=22,face='bold',colour="black"))+
        geom_errorbar(mapping=aes(x=subcluster, ymin=min, ymax=max),width=0.3)
      # scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2),labels=c('-2.0','-1.0','0','1.0','2.0'))
      expr_figure_line[[paste('genecluster',s,sep="")]]=p
      expr_figure_line_bar[[paste('genecluster',s,sep="")]]=p1
    }
    # require(gridExtra)
    library(patchwork)
    expr_line=wrap_plots(expr_figure_line,ncol=5)
    height=(88*round(nrow(markergene)/5))/20
    ggsave(paste('gene_',name,'_',gsub('-','',Sys.Date()),'VnobarV1.png',sep=""),expr_line,width=30,height=height)
    ggsave(paste('gene_',name,'_',gsub('-','',Sys.Date()),'VnobarV1.pdf',sep=""),expr_line,width=30,height=height,limitsize = F)
    expr_line_bar=wrap_plots(expr_figure_line_bar,ncol=5)
    ggsave(paste('gene_',name,'_',gsub('-','',Sys.Date()),'VbarV1.png',sep=""),expr_line_bar,width=30,height=height)
    ggsave(paste('gene_',name,'_',gsub('-','',Sys.Date()),'VbarV1.pdf',sep=""),expr_line_bar,width=30,height=height,limitsize = F)
  }
}

# for(i in 1:length(namelist)){
#   markergene=read.delim(namelist[i],sep="\t",as.is=T,header=T,check.names=F)
#   markergene=markergene[!duplicated(markergene[,1]),]
#   expr1=subset(expr,idents=c('C1.0','C1.4','C14.0','C11.0'))
#   expr1=expr1[markergene[,1],]
#   rownames(markergene)=markergene[,1]
#   used_genes=markergene[,1]
#   # data=as.matrix(expr1[used_genes,][['Spatial']]@data)
#   #### 以下通过AverageExpression计算均值可能有点问题，expm1值的计算是exp(x)-1,可能不太对。
#   # data_mean=AverageExpression(expr1,features = markergene[,1],slot = 'data')
#   # data_mean=data.frame(gene=rownames(data_mean$Spatial),data_mean$Spatial,check.names = F)
#   # colnames(data_mean)=paste('C',colnames(data_mean),sep="")
#   # data_mean[,2:ncol(data_mean)]=log(data_mean[,2:ncol(data_mean)])
#   ###### 重新计算画图的
#   
#   data=t(as.matrix(expr1[['Spatial']]@data))
#   data=data.frame(data,cluster=paste(as.character(Idents(expr1)),sep=""))
#   data_mean=aggregate(data[,1:(ncol(data)-1)],by=list(data$cluster),mean)
#   rownames(data_mean)=data_mean[,1]
#   data_mean=t(data_mean[,-1])
#   data_mean=log(data_mean+1)  ##对数据取log再做归一化
#   # data_mean_scale=data_mean
#   data_mean_scale=t(apply(data_mean,1,scale))
#   colnames(data_mean_scale)=colnames(data_mean)
#   rownames(data_mean_scale)=rownames(data_mean)
#   anno=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/anno_tf_go_kegg_ko.RDS')
#   
#   library(reshape2)
#   x1_data=melt(data_mean_scale)
#   colnames(x1_data)=c('genecluster','subcluster','value')
#   x1_data$subcluster=factor(x1_data$subcluster,levels=c('C1.0','C1.4','C14.0','C11.0'))
#   # x1_data$value=log(x1_data$value)
#   
#   library(dplyr)
#   
#   num=floor(nrow(markergene)/100)
#   genes_split=list()
#   genes_split[[namelist[i]]]=used_genes
#   # genes_split=split(used_genes[1:(num*100)],f = rep(1:num,each=100))
#   # genes_split[[as.character(num+1)]]=used_genes[((num*100)+1):length(used_genes)]
#   
#   for(name in names(genes_split)){
#     print(name)
#     expr_figure_line=list()
#     for(i in genes_split[[name]]){
#       tmp_max_value=data.frame(x1_data)%>%filter(genecluster==i)
#       # print(tmp_max_value)
#       p=ggplot(tmp_max_value,aes(x=subcluster,y=value))
#       p=p+geom_line(data=tmp_max_value,group = 1,mapping=aes(x=subcluster,y=value),size=2.5)+
#         theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
#         labs(title = gsub(':',"|",anno[i,'genesymbol']),x="",y='Log Exp')+geom_point(size = 3)+
#         theme(plot.title = element_text(hjust = 0.5,size=25,face='bold'),
#               text = element_text(size=22,face='bold',colour="black"),
#               axis.text.x = element_text(size=22,angle=90,hjust=1,face='bold',colour="black"),
#               axis.text.y = element_text(size=22,face='bold',colour="black"))
#       # scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2),labels=c('-2.0','-1.0','0','1.0','2.0'))
#       expr_figure_line[[paste('genecluster',i,sep="")]]=p
#     }
#     # require(gridExtra)
#     library(patchwork)
#     expr_line=wrap_plots(expr_figure_line,ncol=5)
#     height=(88*round(length(expr_figure_line)/5))/20
#     ggsave(paste('gene_',name,'_',gsub('-','',Sys.Date()),'V1.png',sep=""),expr_line,width=30,height=height,limitsize = F)
#     ggsave(paste('gene_',name,'_',gsub('-','',Sys.Date()),'V1.pdf',sep=""),expr_line,width=30,height=height,limitsize = F)
#   }
# }