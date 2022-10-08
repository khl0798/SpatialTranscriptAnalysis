##########################################20210924, 和bulk数据相关性比较
## /data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20210924
rna_bulk=read.delim('gene_expression.xls',sep="\t",as.is=T,header=T,check.names=F)
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20210824/combin.subcluster.data.RDS')
rna_sp=as.matrix(expr[['Spatial']]@data)
#### 先过滤掉在所有表达值为0的基因
check_sum=apply(rna_sp,1,function(x)sum(x==0))
rna_sp_rmgenes=rna_sp[check_sum!=ncol(rna_sp),]
rna_sp_scale_mean=apply(rna_sp_rmgenes,1,mean)  #### 31907
#### 读取bulk RNA基因的表达值
rna_bulk_data=rna_bulk[,c(1,5:10)]
## 过滤掉在所有表达值为0的基因   
rna_bulk_check_na=apply(rna_bulk_data[,2:7],1,function(x)sum(x==0))
rna_bulk_data_rm_0=rna_bulk_data[rna_bulk_check_na!=6,]
rna_bulk_data_rm_0_mean=apply(rna_bulk_data_rm_0[,2:7],1,mean)
names(rna_bulk_data_rm_0_mean)=rna_bulk_data_rm_0[,1]   ### 31386
####
used_genes=intersect(names(rna_sp_scale_mean),names(rna_bulk_data_rm_0_mean))  #### 29230
####
####
rna_sp_umi=apply(rna_sp_rmgenes,1,mean)
rna_sp_umi_mean=scale(rna_sp_umi)

cor(rna_sp_umi_mean[used_genes,1],rna_bulk_data_rm_0_mean[used_genes],method = 'pearson')
sp_bulk_rna=data.frame(sp=rna_sp_umi_mean[used_genes,1],bulk=rna_bulk_data_rm_0_mean[used_genes])
sp_bulk_rna[,2]=scale(sp_bulk_rna[,2])
library(ggplot2)
p=ggplot(sp_bulk_rna,aes(sp,bulk))+geom_point()+geom_abline()+
  annotate(geom="text", x=0, y=30, label="correlation:0.825",color="black",
           hjust=0)+theme_bw()+
  labs(x='SpatialRNA(UMI)',y='bulk mRNA-seq(FPKM)')+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave('sp_bulk_corrrelation.png',p,width=6,height=6)
ggsave('sp_bulk_corrrelation.pdf',p,width=6,height=6)
#########单个样本的correlation
### 提取N1样本
rna_sp_rmgenes_sample=data.frame(t(rna_sp_rmgenes),as.character(colnames(expr)),stringsAsFactors = F)
library(dplyr)
colnames(rna_sp_rmgenes_sample)[ncol(rna_sp_rmgenes_sample)]='samples'
rna_sp_rmgenes_sample=t(rna_sp_rmgenes_sample)
####取交集的基因
expr[['stim']]=expr[['orig.ident']]
expr1=subset(expr,stim=='1-N1-1')
rna1=as.matrix(expr1[['Spatial']]@counts)
rna1_check_na=apply(rna1,1,function(x)sum(x==0))
rna1_sp_rm_0=rna1[rna1_check_na!=ncol(rna1),]

name='1-N1-1'
n1_sp=rna_sp_rmgenes[,colnames(expr)[expr[['orig.ident']][,1]==name]]
n1_sp_check_na=apply(n1_sp,1,function(x)sum(x==0))
n1_sp_rm_0=n1_sp[n1_sp_check_na!=ncol(n1_sp),]
n1_sp_mean=apply(n1_sp_rm_0,1,mean)  #### 
n1_sp_scale_mean=scale(n1_sp_mean)  #30484   30296
print(paste('sp genes',nrow(n1_sp_scale_mean),sep=" "))
n1_rna=rna_bulk_data_rm_0[,c(1,grep(substr(name,3,4),colnames(rna_bulk_data_rm_0)))]
n1_bulk_check_na=apply(n1_rna[,2:4],1,function(x)sum(x==0))
n1_bulk_data_rm_0=n1_rna[n1_bulk_check_na!=3,]  ### 30017
n1_bulk_data_rm_0_mean=apply(n1_bulk_data_rm_0[,2:4],1,mean)
n1_bulk_data_rm_0_mean=scale(n1_bulk_data_rm_0_mean)
rownames(n1_bulk_data_rm_0_mean)=n1_bulk_data_rm_0[,1]   ### 30017   30362
print(paste('bulk genes',nrow(n1_bulk_data_rm_0_mean),sep=" "))

used_genes=intersect(rownames(n1_sp_scale_mean),rownames(n1_bulk_data_rm_0_mean))  #### 27736  28005
print(paste('union genes',length(used_genes),sep=" ")) 

sp_bulk_rna=data.frame(sp=n1_sp_scale_mean[used_genes,1],bulk=n1_bulk_data_rm_0_mean[used_genes,1])
corr=cor(sp_bulk_rna[,1],sp_bulk_rna[,2],method='pearson')
label=paste('correlation:',round(corr,3),sep="")
print(label)
p=ggplot(sp_bulk_rna,aes(sp,bulk))+geom_point()+geom_abline()+
  annotate(geom="text", x=0, y=floor(max(max(sp_bulk_rna[,1]),max(sp_bulk_rna[,2]))), label=label,color="black",
           hjust=0)+theme_bw()+
  labs(x='SpatialRNA(UMI)',y='bulk mRNA-seq(FPKM)')+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
##   	scale_x_log10() +    	stat_smooth(method=lm) +
ggsave(paste(name,'_sp_bulk_corrrelation.png',sep=""),p,width=6,height=6)
ggsave(paste(name,'_sp_bulk_corrrelation.pdf',sep=""),p,width=6,height=6)
## 
library(eulerr)
vd=euler(c(VisiumRNA=2677,bulkRNA=2156,"VisiumRNA&bulkRNA"=29230))
p=plot(vd,fills=list(fill=c('pink','skyblue'),alpha=0.6),
       labels=list(col='black'),edges=FALSE,quantities = TRUE,
       adjust_labels = TRUE,legend = T)
ggsave('E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20210924/BC210523-1/BC210523-1venn.png',p,
       width=8.5,height=6)
ggsave('E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20210924/BC210523-1/BC210523-1venn.pdf',p,
       width=8.5,height=6)
vd=euler(c(VisiumRNA=2748,bulkRNA=2281,"VisiumRNA&bulkRNA"=27736))
p=plot(vd,fills=list(fill=c('pink','skyblue'),alpha=0.6),
       labels=list(col='black'),edges=FALSE,quantities = TRUE,
       adjust_labels = TRUE,legend = T)
ggsave('E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20210924/BC210523-1/1-N1-1_venn.png',p,
       width=8.5,height=6)
ggsave('E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20210924/BC210523-1/1-N1-1_venn.pdf',p,
       width=8.5,height=6)
vd=euler(c(VisiumRNA=2291,bulkRNA=2357,"VisiumRNA&bulkRNA"=28005))
p=plot(vd,fills=list(fill=c('pink','skyblue'),alpha=0.6),
       labels=list(col='black'),edges=FALSE,quantities = TRUE,
       adjust_labels = TRUE,legend = T)
ggsave('E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20210924/BC210523-1/2-N2-1_venn.png',p,
       width=8.5,height=6)
ggsave('E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20210924/BC210523-1/2-N2-1_venn.pdf',p,
       width=8.5,height=6)