############################20220422，大群umap颜色修改及点大小更改
setwd('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220421_chongxinhuatu/umap')
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')
Idents(expr)=paste('C',Idents(expr),sep="")
levels_name=c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C11.0','C11.1','C14.0','C14.1')
new_sub_c1_c14=as.character(Idents(expr))
names(new_sub_c1_c14)=colnames(expr)
new_sub_c1_c14[!new_sub_c1_c14%in%levels_name]='other'
Idents(expr)=factor(new_sub_c1_c14[colnames(expr)],levels=c(levels_name,'other'))
labels=c('C1.0 PIs','C1.1 PXy',"C1.2 PPh","C1.3 PXy IN1","C1.4 MC","C1.5 PPh IN1","C11.0 CIs",
         "C11.1 DXy","C14.0 PLIs","C14.1 P",'other')
#"#3A76AF"
cluster9_cor=c("#FCC6BF","#E96E0B", "#F380BC", "#C81618", 
               "#B17BA6","#AF00FF",
               "#7EC97E", "#4BBDC8","#0000FF",'#FFD700','grey')
#####42,178,80  #2AB250
##### 37 189 223  "#25BDDF"
###### 81 219 253 '#51DBFD'  
#######黄色FFD700
pretty_unexpanded <- function(x,n=5){
  #expand_limit
  if(x[1]<=0){
    r2 <- x + c(-x[1],x[1])
  }else{
    r2 <- x + c((x[2]-x[1])*0.04545455,-(x[2]-x[1])*0.04545455)
  }
  pout <-  pretty(r2,n) 
  pout
}
library(scales)
S_sqrt <- function(x){sign(x)*(abs(x)^1.3)}
IS_sqrt <- function(x){(x^2)*sign(x)}
S_sqrt1_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt,breaks=c(-12,-6,0,6))
pumap_8clusters <- DimPlot(expr, reduction = "umap",pt.size=5, repel =TRUE,
                           label = TRUE,label.size=12,shuffle=TRUE)+
  theme(legend.text = element_text(size=28,family = 'Arial',face='bold',colour = 'black'),
        legend.title = element_text(face = "bold", size = 38, colour = "black"),
        axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
        axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
        title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
        strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
  labs(color = "subcluster")+
  guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  scale_color_manual(values=cluster9_cor,labels=labels)+
  scale_x_continuous(trans="S_sqrt1", name="UMAP_1",breaks=c(-12,-6,0,6),labels=c(-12,-6,0,6))
ggsave(paste(gsub('-','',Sys.Date()),'_C1_11_14','UMAP_3.png',sep=""),pumap_8clusters,
       height = 17,width = 25)
pumap_8clusters2 <- DimPlot(expr, reduction = "umap",pt.size=5, repel =TRUE,
                            label = TRUE,label.size=12,shuffle=TRUE)+
  theme(legend.text = element_text(size=28,face='bold',colour = 'black'),
        legend.title = element_text(face = "bold", size = 38, colour = "black"),
        axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
        axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
        title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
        strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
  labs(color = "subcluster")+
  guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  scale_color_manual(values=cluster9_cor,labels=labels)+
  scale_x_continuous(trans="S_sqrt1", name="UMAP_1",breaks=c(-12,-6,0,6),labels=c(-12,-6,0,6))
ggsave(paste(gsub('-','',Sys.Date()),'_C1_11_14','UMAP_3.pdf',sep=""),pumap_8clusters2,
       height = 17,width = 25)
# S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
# IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)  #breaks=c(-12,-6,0,2,4,6,8)

pumap_8clusters <- DimPlot(expr, reduction = "umap",pt.size=5, repel =TRUE,
                           label = TRUE,label.size=12,shuffle=TRUE)+
  theme(legend.text = element_text(size=28,family = 'Arial',face='bold',colour = 'black'),
        legend.title = element_text(face = "bold", size = 38, colour = "black"),
        axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
        axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
        title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
        strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
  labs(color = "subcluster")+
  # scale_x_discrete("type",expand=c(0,20))
  guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  # scale_x_discrete(breaks=c(-12,-6,0,4,8), labels=c(-12,-6,0,4,8))+
  scale_color_manual(values=cluster9_cor,labels=labels)+
  scale_x_continuous(breaks=c(-12,-6,0,6,10),labels=c(-12,-6,0,6,10),
                     limits = c(-12,10))
ggsave(paste(gsub('-','',Sys.Date()),'_C1_11_14','UMAP_1.png',sep=""),pumap_8clusters,
       height = 17,width = 25)
#####
cluster9_cor=c("#99CC00FF","#FFCD00FF","#6EE2FFFF",'#e7de2c',"#7bf1a8","#99991EFF",
               "#38b000", "#e77f0f","#2196f3","#d1b1cb",'grey')
####C1 的6个亚群； C1.0 PIs， C1.1 PXy; C1.2 PPh;  C1.3 PXy IN1; C1.5 PPh IN1; C1.4 MC; C11.0 CIs; C11.1 DXy;  C14.0 PLIs; C14.1 Ph;

pumap_8clusters <- DimPlot(expr, reduction = "umap",pt.size=5, repel =TRUE,
                           label = TRUE,label.size=12,shuffle=TRUE)+
  theme(legend.text = element_text(size=28,face='bold',colour = 'black'),
        legend.title = element_text(face = "bold", size = 38, colour = "black"),
        axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
        axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
        title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
        strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
  labs(color = "subcluster")+
  guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  scale_color_manual(values=cluster9_cor,labels=labels)+
  scale_x_continuous(breaks=c(-12,-6,0,6,10),labels=c(-12,-6,0,6,10),
                     limits = c(-12,10))
ggsave(paste(gsub('-','',Sys.Date()),'_C1_11_14','UMAP_1.pdf',sep=""),pumap_8clusters,
       height = 17,width = 25)
pumap_8clusters <- DimPlot(expr, reduction = "umap",pt.size=5, repel =TRUE,
                           label = TRUE,label.size=12,shuffle=TRUE)+
  theme(legend.text = element_text(size=28,face='bold',colour = 'black'),
        legend.title = element_text(face = "bold", size = 38, colour = "black"),
        axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
        axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
        title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
        strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
  labs(color = "subcluster")+
  guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  scale_color_manual(values=cluster9_cor,labels=labels)+
  scale_x_continuous(trans="S_sqrt", name="UMAP_1",breaks=c(-12,-6,0,6),
                     labels=c(-12,-6,0,6))
ggsave(paste(gsub('-','',Sys.Date()),'_C1_11_14','UMAP_2.pdf',sep=""),pumap_8clusters,height = 17,width = 25)
#######################切片图
total_expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')

expr_sub_1_11_14=subset(total_expr,
                        idents=paste('C',c('1.0','1.1','1.2','1.3','1.4','1.5','11.0','11.1','14.0','14.1'),sep=""))

# cluster9_cor_new=c("#FCC6BF","#E96E0B", "#F380BC", "#C81618", 
#                    "#B17BA6","#AF00FF",
#                    "#7EC97E", "#4BBDC8","#0000FF",'#FFD700')
cluster9_cor_new=c("#99CC00FF", "#FFCD00FF", "#6EE2FFFF", "#e7de2c",   "#7bf1a8",   "#99991EFF",
                   "#38b000",   "#e77f0f",   "#2196f3",   "#d1b1cb")
cor_list=as.list(cluster9_cor_new)
names(cor_list)=paste('C',c('1.0','1.1','1.2','1.3','1.4','1.5','11.0','11.1','14.0','14.1'),sep="")
cor_list[['null']]='#FFFFFF'
total_new_info=as.character(Idents(total_expr))
names(total_new_info)=colnames(total_expr)
# total_new_info[colnames(expr_sub_1_11_14)]=as.character(Idents(expr_sub_1_11_14))
total_new_info[!total_new_info%in%levels(Idents(expr_sub_1_11_14))]='null'
total_expr[['total_new_info']]=factor(total_new_info[colnames(total_expr)],
                                      levels=c(paste('C',c('1.0','1.1','1.2','1.3','1.4','1.5','11.0','11.1','14.0','14.1'),sep=""),'null'))
tmp_expr=total_expr
Idents(tmp_expr)=tmp_expr[['total_new_info']]

p1=SpatialPlot(tmp_expr,crop = TRUE, image.alpha=0,slot='data',combine = FALSE)
p1[[1]]=p1[[1]]+scale_fill_manual(values=c(unlist(cor_list[levels(subset(tmp_expr,orig.ident=='1-N1-1'))])))
p1[[2]]=p1[[2]]+scale_fill_manual(values=c(unlist(cor_list[levels(subset(tmp_expr,orig.ident=='2-N2-1'))])))
p1[[1]]$theme$legend.position='none'
p1[[2]]$theme$legend.position='none'
ggsave(paste(gsub('-','',Sys.Date()),'_4clusters','img.png',sep=""),CombinePlots(p1),width=12,height=6)
ggsave(paste(gsub('-','',Sys.Date()),'_4clusters','img.pdf',sep=""),CombinePlots(p1),width=12,height=6)# c1_data=subset(tmp_expr,idents=c('1.0','1.1','1.2','1.3','1.4','1.5'))
###############################17个cluster的大群umap图
# Idents(tmp_expr)=factor(Idents(tmp_expr),
#                         levels=c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C11.0','C11.1','C14.0','C14.1','null'))
# Idents(tmp_expr)=factor(tmp_expr$total_new_info,levels=c('null',
#                                                          'C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C11.0','C11.1','C14.0','C14.1'))
# labels=c("other","C1.0 PIs",     "C1.1 PXy",     "C1.2 PPh" ,    "C1.3 PXy IN1" ,"C1.4 MC"  ,   
#          "C1.5 PPh IN1", "C11.0 CIs" ,   "C11.1 DXy" ,   "C14.0 PLIs" ,  "C14.1 P")
# cluster9_cor=c("grey","#99CC00FF", "#FFCD00FF", "#6EE2FFFF" ,"#e7de2c" ,  "#7bf1a8" ,  "#99991EFF",
#                "#38b000",   "#e77f0f",   "#2196f3",   "#d1b1cb")
# pumap_8clusters <- DimPlot(tmp_expr, reduction = "umap",pt.size=5, repel =TRUE,
#                            label = TRUE,label.size=12,shuffle=TRUE)+
#   theme(legend.text = element_text(size=28,family = 'Arial',face='bold',colour = 'black'),
#         legend.title = element_text(face = "bold", size = 38, colour = "black"),
#         axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
#         axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
#         title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
#         strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
#   labs(color = "subcluster")+
#   # scale_x_discrete("type",expand=c(0,20))
#   guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
#   guides(color = guide_legend(override.aes = list(size = 10)))+
#   # scale_x_discrete(breaks=c(-12,-6,0,4,8), labels=c(-12,-6,0,4,8))+
#   scale_color_manual(values=cluster9_cor,labels=labels)+
#   scale_x_continuous(breaks=c(-12,-6,0,6,10),labels=c(-12,-6,0,6,10),
#                      limits = c(-12,10))
# ggsave(paste(gsub('-','',Sys.Date()),'_17clusters','UMAP_2.png',sep=""),pumap_8clusters,
#        height = 17,width = 25)
# pumap_8clusters1<- DimPlot(tmp_expr, reduction = "umap",pt.size=5, repel =TRUE,
#                            label = TRUE,label.size=12,shuffle=TRUE)+
#   theme(legend.text = element_text(size=28,face='bold',colour = 'black'),
#         legend.title = element_text(face = "bold", size = 38, colour = "black"),
#         axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
#         axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
#         title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
#         strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
#   labs(color = "subcluster")+
#   # scale_x_discrete("type",expand=c(0,20))
#   guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
#   guides(color = guide_legend(override.aes = list(size = 10)))+
#   # scale_x_discrete(breaks=c(-12,-6,0,4,8), labels=c(-12,-6,0,4,8))+
#   scale_color_manual(values=cluster9_cor,labels=labels)+
#   scale_x_continuous(breaks=c(-12,-6,0,6,10),labels=c(-12,-6,0,6,10),
#                      limits = c(-12,10))
# ggsave(paste(gsub('-','',Sys.Date()),'_17clusters','UMAP_2.pdf',sep=""),pumap_8clusters1,
#        height = 17,width = 25)
# 
# ########只展示1.0，1.4，14.0，11.0 这4个亚群
# newclusters=as.character(Idents(tmp_expr))
# names(newclusters)=colnames(tmp_expr)
# newclusters[!newclusters%in%c('C1.0','C1.4','C14.0','C11.0')]='null'
# tmp_expr[['newclusters']]=factor(newclusters[colnames(tmp_expr)],levels=c('C1.0','C1.4','C11.0','C14.0','null'))
# Idents(tmp_expr)=factor(tmp_expr$newclusters,levels=c('null','C1.0','C1.4','C11.0','C14.0'))
# # labels
# # [1] "C1.0 PIs"     "C1.1 PXy"     "C1.2 PPh"     "C1.3 PXy IN1" "C1.4 MC"     
# # [6] "C1.5 PPh IN1" "C11.0 CIs"    "C11.1 DXy"    "C14.0 PLIs"   "C14.1 P"     [11] "other"
# tmp_labels=c('other','C1.0 PIs','C1.4 MC',"C11.0 CIs","C14.0 PLIs")
# tmp_colors=c("grey","#99CC00FF","#7bf1a8","#38b000","#2196f3")
# pumap_8clusters <- DimPlot(tmp_expr, reduction = "umap",pt.size=5, repel =TRUE,order=TRUE,
#                            label = TRUE,label.size=12,shuffle=TRUE)+
#   theme(legend.text = element_text(size=28,family = 'Arial',face='bold',colour = 'black'),
#         legend.title = element_text(face = "bold", size = 38, colour = "black"),
#         axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
#         axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
#         title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
#         strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
#   labs(color = "subcluster")+
#   # scale_x_discrete("type",expand=c(0,20))
#   guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
#   guides(color = guide_legend(override.aes = list(size = 10)))+
#   # scale_x_discrete(breaks=c(-12,-6,0,4,8), labels=c(-12,-6,0,4,8))+
#   scale_color_manual(values=tmp_colors,labels=tmp_labels)+
#   scale_x_continuous(breaks=c(-12,-6,0,6,10),labels=c(-12,-6,0,6,10),
#                      limits = c(-12,10))
# ggsave(paste(gsub('-','',Sys.Date()),'_4clusters','UMAP_2.png',sep=""),pumap_8clusters,
#        height = 17,width = 25)
# pumap_8clusters1<- DimPlot(tmp_expr, reduction = "umap",pt.size=5, repel =TRUE,order=TRUE,
#                            label = TRUE,label.size=12,shuffle=TRUE)+
#   theme(legend.text = element_text(size=28,face='bold',colour = 'black'),
#         legend.title = element_text(face = "bold", size = 38, colour = "black"),
#         axis.text.x = element_text(face = "bold", size = 28,colour = 'black'), 
#         axis.text.y = element_text(face = "bold",  size = 28,colour = 'black'), 
#         title = element_text(size = 30,face='bold',colour = 'black'), legend.key.height = unit(3,"cm"),
#         strip.text = element_text(face = "bold", size = 26,colour = 'black'))+
#   labs(color = "subcluster")+
#   # scale_x_discrete("type",expand=c(0,20))
#   guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
#   guides(color = guide_legend(override.aes = list(size = 10)))+
#   # scale_x_discrete(breaks=c(-12,-6,0,4,8), labels=c(-12,-6,0,4,8))+
#   scale_color_manual(values=tmp_colors,labels=tmp_labels)+
#   scale_x_continuous(breaks=c(-12,-6,0,6,10),labels=c(-12,-6,0,6,10),
#                      limits = c(-12,10))
# ggsave(paste(gsub('-','',Sys.Date()),'_4clusters','UMAP_2.pdf',sep=""),pumap_8clusters1,
#        height = 17,width = 25)
