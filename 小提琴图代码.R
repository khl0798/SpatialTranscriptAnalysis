library(Seurat)
library(cowplot)
####/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220530_vlnplot
# expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')
library(ggplot2)
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')
levels=c("C1.5","C1.2", 'C1.0',  "C1.1",  "C1.3", "C1.4",    "C12", "C14.1","C14.0", "C11.0", "C11.1", "C8")
sub_expr=subset(expr,idents=levels)
# cluster_color_list
########################先获取所有的基因的小提琴的data信息，后面可以直接获取信息等等
tmp1=FetchData(sub_expr,vars = rownames(sub_expr))### vars是基因
library(dplyr)
# tmp=FetchData(sub_expr,vars = gene1[,1])
library(reshape2)
x=tmp1%>%tibble::rownames_to_column('barcodes')%>%melt('barcodes')
x$cluster=Idents(expr)[x[,1]]

# x$cluster=factor(as.character(x$cluster),levels=levels)

#################
#小提琴图的颜色高亮展示每个基因在C14.1, C14.0, C11.0, C11.1 中的表达情况，颜色就用各自cluster的颜
colnames(x)[2]='gene'
x=data.frame(x,new='new',vlnplot_value=0)  ###x 添加2列，new列（实际中没有用到这一列的信息），vlnplot_value列是决定是否画小提琴的列
x$cluster=as.character(x$cluster)
x$new=as.character(x$new)
cluster9_cor_new=c("#99CC00FF","#FFCD00FF","#6EE2FFFF",'#e7de2c',"#7bf1a8","#188cbe",
                   "#FFBFE5" , "#d1b1cb","#2196f3","#38b000", "#fcbc5d","#F79D1E")
names(cluster9_cor_new)=levels
x$gene=as.character(x$gene)

gene1=read.delim('genefeatures.txt',sep="\t",as.is=T,header=T,check.names=F) ###读取需要画图的基因文件
gene1=gene1[!duplicated(gene1[,1]),]
rownames(gene1)=gene1[,1]
########选择基因数据信息，准备最后画图需要用到的数据
vlnplot_tmp=x[x$gene%in%gene1[,1],]
vlnplot_tmp$gene=as.character(vlnplot_tmp$gene)
vlnplot_tmp$newgene=factor(gene1[vlnplot_tmp$gene,3],levels=unique(gene1$genename))  ###按照基因顺序设置levels
vlnplot_tmp$vlnplot_cluster=vlnplot_tmp$cluster
#######
tmp_levels=setdiff(levels,c('C1.3','C1.5')) 
# vlnplot_tmp$cluster=factor(vlnplot_tmp$cluster,levels=levels,order=T)
vlnplot_tmp$vlnplot_value=vlnplot_tmp$value  ###赋值vlnplot_value这一列，取的是value列的信息
vlnplot_tmp$vlnplot_value[!vlnplot_tmp$cluster%in%c("C1.2","C1.0","C1.1","C1.4","C12","C14.1","C14.0","C11.0", "C11.1","C8")]=NA ###不画小提琴的亚群的值设置为NA
new_data=vlnplot_tmp
new_data$cluster=factor(new_data$cluster,levels=levels,order=F)
new_data$cluster=factor(new_data$cluster,levels=levels,order=F)
new_data$vlnplot_cluster='new' ### vlnplot_cluster控制画小提琴的亚群，不画的小提琴的亚群都设置为new
new_data$vlnplot_cluster[as.character(new_data$cluster)%in%tmp_levels]=as.character(new_data$cluster[as.character(new_data$cluster)%in%tmp_levels])
new_data$vlnplot_cluster=factor(new_data$vlnplot_cluster,levels=c(tmp_levels,'new'))

a <- ggplot() + #ggnewscale::new_scale()+ #stroke = 0.3
  geom_violin(data=new_data,mapping=aes(x=cluster,y=vlnplot_value,fill=vlnplot_cluster), ## vlnplot_value列是画小提琴的
              scale = "width", trim = TRUE,width=1) +
  geom_jitter(data=new_data, mapping=aes(x=factor(cluster),y=value), ### value列是所有的数据，画jitter点
              width = 0.2,size=0.5,stroke=0.3)+theme_bw()+ #ggnewscale::new_scale("fill")+
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
  scale_fill_manual(values=cluster9_cor_new[tmp_levels]) ###控制颜色，画小提琴的亚群才有颜色，不画小提琴的亚群没有颜色，geom_jitter默认的是黑色的点
###,hjust=1
a1 <- a+theme(axis.text.y = element_text(vjust = 0,size=25,colour='black',face='bold'))
####30个基因长度是27，
height=dim(gene1)[1]*30/27
####3带基因名字的是width=20
ggsave(paste(name,'_lnplot_',gsub('\\-','',Sys.Date()),'.png',sep=""),a1,width=16,height=height,limitsize = F)
ggsave(paste(name,'_lnplot_',gsub('\\-','',Sys.Date()),'.pdf',sep=""),a1,width=16,height=height,limitsize = F)