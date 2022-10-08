# setwd('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220311_genemonocle_C1_C12_8_11_14')
# load('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220113_C1umap_subcluster_results/monocle_C1/C1_20220118_HSMM_myo.Rdata')
library(monocle)
library('VisiumSpatial')
setwd('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220506_umap_vlnplot')
# setwd('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220316_genemonocle_C1_C12_8_11_14_finalversion')
#HSMM_myo1=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220311_genemonocle_C1_C12_8_11_14/HSMM_myo1.RDS')
load('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_11_8_HSMM_myo.Rdata')
HSMM_myo1=HSMM_myo
# C1_genes=read.delim(file.path('C12_8_11_14','genes.txt'),sep="\t",as.is=T,header=T,check.names=F)
C1_genes=read.delim(file.path('C1_11_8_linecurve.txt'),sep="\t",as.is=T,header=T,check.names=F)
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/1_3_8_9_11_12_13_14.RDS')

x2=pData(HSMM_myo1)
new_expr=as_matrix(expr[['Spatial']]@data)
new_expr=new_expr[,rownames(x2)]

for(i in 1:nrow(C1_genes)){
  new_name=C1_genes[i,2]
  x2=cbind(x2,new_expr[C1_genes[i,1],])
}
colnames(x2)[6:ncol(x2)]=C1_genes[1:nrow(C1_genes),1]
# for(i in 1:nrow(new_tmp_data)){
#   x2=cbind(x2,new_expr[new_tmp_data[i,'gene'],])
# }
# colnames(x2)[70:ncol(x2)]=new_tmp_data[,'gene']


HSMM_myo1@phenoData@data=x2

fi=list()
for(i in 1:nrow(C1_genes)){
  p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = C1_genes[i,1], 
                          cell_size = 1.2) + 
    theme(axis.title = element_text(face = "bold", size = 22, colour = "black"), legend.position = "right", 
          legend.text = element_text(face = "bold", size = 22, 
                                     colour = "black"), legend.title = element_text(face = "bold", 
            size = 25, colour = "black"), axis.text.x = element_text(face = "bold", colour='black',
          size = 22), axis.text.y = element_text(face = "bold", size = 22,colour='black'), 
          title = element_text(size = 20,face='bold',family='Arial'), 
          strip.text = element_text(face = "bold", size = 22),legend.key.height = unit(1.3,"cm"))+
    labs(title=C1_genes[i,2],family = "Arial",face='bold',col='black',size=20)+
    theme(plot.title = element_text(hjust = 0.5,vjust = 0.3)) + #也就加上这一行
    guides(fill=guide_legend(label.hjust=1.5,label.vjust=1.5))+
    guides(shape = guide_legend(override.aes = list(size = 3)))+
    #low="#CDC9C9",mid='#FFFF00', high = "#FF0000"
    scale_colour_gradientn(colours=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ),
                           name=NULL)
  
  # scale_fill_gradientn(colors=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ))
  # scale_color_viridis_c()
  fi[[C1_genes[i,1]]]=p1
}
fi_pdf=list()
for(i in 1:nrow(C1_genes)){
  p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = C1_genes[i,1], 
                          cell_size = 1.2) + 
    theme(axis.title = element_text(face = "bold", size = 22, colour = "black"), legend.position = "right", 
          legend.text = element_text(face = "bold", size = 22, 
                                     colour = "black"), legend.title = element_text(face = "bold", 
     size = 25, colour = "black"), axis.text.x = element_text(face = "bold", colour='black',
      size = 22), axis.text.y = element_text(face = "bold", size = 22,colour='black'), 
          title = element_text(size = 20,face='bold'), 
          strip.text = element_text(face = "bold", size = 22),legend.key.height = unit(1.3,"cm"))+
    labs(title=C1_genes[i,2],face='bold',col='black',size=20)+
    theme(plot.title = element_text(hjust = 0.5,vjust = 0.3)) + #也就加上这一行
    guides(fill=guide_legend(label.hjust=1.5,label.vjust=1.5))+
    guides(shape = guide_legend(override.aes = list(size = 3)))+
    #low="#CDC9C9",mid='#FFFF00', high = "#FF0000"
    scale_colour_gradientn(colours=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ),
                           name=NULL)
  
  # scale_fill_gradientn(colors=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ))
  # scale_color_viridis_c()
  fi_pdf[[C1_genes[i,1]]]=p1
}
library(reshape2)
expr_plotdata=cds_exprs
tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%C1_genes[,1],]

tmp_figure$subcluster=tmp_figure$clusters
tmp_max_value=aggregate(tmp_figure$expectation,by=list(tmp_figure$clusters,tmp_figure$gene_short_name),median)
bar_compute<-function(data){
  upper=median(data)+sd(data)
  lower=median(data)-sd(data)
  return(c(upper,lower))
}
tmp_max=aggregate(tmp_figure$expectation,by=list(tmp_figure$clusters,tmp_figure$gene_short_name),bar_compute)

tmp_max_value=cbind(tmp_max_value,tmp_max[,3])
colnames(tmp_max_value)=c('subcluster','gene','value','tmp_max','tmp_min')
# tmp_max_value$subcluster=factor(as.character(tmp_max_value$subcluster),levels=c("C1.1","C1.3","C1.0","C1.2","C11.0",'C11.1'))
tmp_max_value$subcluster=factor(as.character(tmp_max_value$subcluster),levels=c("C1.0",'C1.4','C11.0','C11.1','C8'))

# tmp_max_value$genecluster=tmp_figure[match(tmp_max_value$gene,tmp_figure$gene_short_name),'newgenecluster']
library(dplyr)
# tmp_max_value=tmp_max_value[order(tmp_max_value$genecluster),]
genes_split=list()
genes_split[['1']]=C1_genes[,1]
rownames(C1_genes)=C1_genes[,1]
for(name in names(genes_split)){
  print(name)
  expr_figure_line=list()
  expr_figure_line_bar=list()
  for(i in genes_split[[name]]){
    tmp_max_value1=data.frame(tmp_max_value)%>%filter(gene==i)
    # print(tmp_max_value)
    p=ggplot(tmp_max_value1,aes(x=subcluster,y=value))
    p=p+geom_line(data=tmp_max_value1,group = 1,mapping=aes(x=subcluster,y=value),size=2.5)+
      theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))+
      labs(title = C1_genes[i,2],x="",y='Log Exp')+geom_point(size = 3)+
      theme(plot.title = element_text(hjust = 0.5,size=25,face='bold'),
            text = element_text(size=22,face='bold',colour="black"),
            axis.text.x = element_text(size=22,angle=90,hjust=1,face='bold',colour="black"),
            axis.text.y = element_text(size=22,face='bold',colour="black"))
      # scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2),labels=c('-2.0','-1.0','0','1.0','2.0'))
    p1=p+geom_errorbar(data=tmp_max_value1,mapping=aes(x=subcluster, ymin=tmp_min, ymax=tmp_max),width=0.3)
    expr_figure_line[[paste(i,sep="")]]=p
    expr_figure_line_bar[[paste(i,sep="")]]=p1
  }
  # require(gridExtra)
  # library(patchwork)
  # expr_line=wrap_plots(expr_figure_line,ncol=5)
  # expr_line=grid.arrange(expr_figure_line[[1]],expr_figure_line[[2]],expr_figure_line[[3]],expr_figure_line[[4]],
  #                        expr_figure_line[[5]],expr_figure_line[[6]],
  #                        ncol = 3)
  # height=(88*round(length(expr_figure_line)/5))/20
  # ggsave(paste('gene_',name,'_',gsub('-','',Sys.Date()),'Vnobar.png',sep=""),expr_line,width=30,height=height)
  # ggsave(paste('gene_',name,'_',gsub('-','',Sys.Date()),'Vnobar.pdf',sep=""),expr_line,width=30,height=height,limitsize = F)
}


## add gene_line 图形放在一起
expr_line=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220310_genecurveline_final/genes_20220316/expr_figure_line.12genes.RDS')
library(patchwork)
new_figure_list=list()
new_figure_list=c(fi[1:5],expr_figure_line[1:5],fi[6:10],expr_figure_line[6:10])
# new_figure_list=c(fi[11:12],expr_line[11:12])
figure_list=wrap_plots(new_figure_list,ncol=5)
ggsave('C1_11_8_10genes_20220506.png',figure_list,width=35,height=20)
new_figure_list=list()
new_figure_list=c(fi_pdf[1:5],expr_figure_line[1:5],fi_pdf[6:10],expr_figure_line[6:10])
# new_figure_list=c(fi[11:12],expr_line[11:12])
figure_list=wrap_plots(new_figure_list,ncol=5)
ggsave('C1_11_8_10genes_20220506.pdf',figure_list,width=35,height=20)

new_figure_list=list()
new_figure_list=c(fi[1:5],expr_figure_line_bar[1:5],fi[6:10],expr_figure_line_bar[6:10])
# new_figure_list=c(fi[11:12],expr_line[11:12])
figure_list=wrap_plots(new_figure_list,ncol=5)
ggsave('C1_11_8_10genes_20220506_bar.png',figure_list,width=35,height=20)
new_figure_list=list()
new_figure_list=c(fi_pdf[1:5],expr_figure_line_bar[1:5],fi_pdf[6:10],expr_figure_line_bar[6:10])
# new_figure_list=c(fi[11:12],expr_line[11:12])
figure_list=wrap_plots(new_figure_list,ncol=5)
ggsave('C1_11_8_10genes_20220506_bar.pdf',figure_list,width=35,height=20)


fi=list()
for(i in 1:nrow(C1_genes)){
  p1=plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = C1_genes[i,1], 
                          cell_size = 1.2) + 
    theme(axis.title = element_text(face = "bold", size = 22, colour = "black"), legend.position = "right", 
          legend.text = element_text(face = "bold", size = 22, 
                                     colour = "black"), 
          legend.title = element_text(face = "bold", 
                                      size = 25, colour = "black"), axis.text.x = element_text(face = "bold", colour='black',
                                                                                               size = 22), axis.text.y = element_text(face = "bold", size = 22,colour='black'), 
          title = element_text(size = 20,face='bold'), 
          strip.text = element_text(face = "bold", size = 22),legend.key.height = unit(1.3,"cm"))+
    labs(title=C1_genes[i,2],face='bold',col='black',size=20)+
    theme(plot.title = element_text(hjust = 0.5,vjust = 0.3)) + #也就加上这一行
    guides(fill=guide_legend(label.hjust=1.5,label.vjust=1.5))+
    guides(shape = guide_legend(override.aes = list(size = 3)))+
    #low="#CDC9C9",mid='#FFFF00', high = "#FF0000"
    scale_colour_gradientn(colours=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ),
                           name=NULL)
  
  # scale_fill_gradientn(colors=c('#EEEEE0','#000080',"#FDAE61", "#F46D43", "#D53E4F" ,"#9E0142" ))
  # scale_color_viridis_c()
  fi[[C1_genes[i,1]]]=p1
}
new_figure_list=list()
# new_figure_list=c(fi[1:5],expr_line[1:5],fi[6:10],expr_line[6:10],fi[11:12],expr_line[11:12])
new_figure_list=c(fi[11:12],expr_line[11:12])
figure_list=wrap_plots(new_figure_list,ncol=2)
ggsave('2genes20220316.pdf',figure_list,width=14,height=10)

# figure_list=wrap_plots(new_