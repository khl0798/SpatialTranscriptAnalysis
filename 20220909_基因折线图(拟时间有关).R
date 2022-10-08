# 5组基因，每组做3个折线图，横坐标如下，用K-means的 方法做
# 折线图的横坐标顺序：C12; C14.1;C14.0, C11.0, C11.1; C8
# 折线图的横坐标顺序： C1.0;  C1.2;  C1.5; C14.0, C14.1;,C12;
# 
# 折线图的横坐标顺序：C1.0;  C1.3; C1.4;  C11.0, C11.1; C8
library(Seurat)
library(monocle)
library(VisiumSpatial)
library(ggplot2)
library(patchwork)

# levels=c("C1.5","C1.2", "C1.0", "C1.3","C1.1", "C1.4")
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
      theme(plot.title = element_text(hjust = 0.5,size=30,face='bold'),
            text = element_text(size=26,face='bold',colour="black"),
            axis.text.x = element_text(size=26,angle=90,hjust=1,face='bold',colour="black"),
            axis.text.y = element_text(size=26,face='bold',colour="black"))+scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2),
                              labels=c('-2.0','-1.0','0','1.0','2.0'))
    expr_figure_line[[gene_express[i,1]]]=p
  }
  return(expr_figure_line)
}
######################20220909,C12_114_11_8拟时序画图 # 拟时序4的横坐标顺序：C1.0;  C1.4;  C11.0, C11.1; C8	
# indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220125_1.0_1.4_11.0_11.1_8/1_11_8_curveline'
# load(file.path(indir,'20220201.rdata'))
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/monocle_20220907/C12_14_11_8'
namelist=list.files(indir,pattern='*txt$')
s=1
levels=c('C12','C14.1','C14.0','C11.0','C11.1','C8')
# gene_express=c()
for(s in seq(1,length(namelist))){
  print(s)
  tmp_gene_express=read.delim(file.path(getwd(),namelist[s]),sep="\t",as.is=T,header=T,check.names=F)
  # if(namelist[s]=='vlnplot20220824_sheet1V1.txt'){
  #   tmp_gene_express$new=paste(tmp_gene_express[,1],tmp_gene_express[,2],sep="|")
  #   tmp_gene_express=tmp_gene_express[,c(1,3)]
  # }
  # if(ncol(tmp_gene_express)>2){
  #   tmp_gene_express=tmp_gene_express[,c(1,2)]
  # }
  tmp_gene_express=tmp_gene_express[!duplicated(tmp_gene_express[,1]),]
  print(paste('gene number is',dim(tmp_gene_express)[1],'and total gene in ',sum(tmp_gene_express[,1]%in%expr_plotdata$gene_short_name)))
  colnames(tmp_gene_express)=c('gene','geneID')
  rownames(tmp_gene_express)=tmp_gene_express[,1]
  gene_express=tmp_gene_express
  tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%gene_express[,1],]
  tmp_figure$clusters=factor(tmp_figure$clusters,levels=levels)
  c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
  figure_list=wrap_plots(c1_linelist,ncol=5)
  id=floor(dim(gene_express)[1]/5)
  if(id*5<dim(gene_express)[1]){
    id=id+1
  }
  cluster_name=paste(paste(levels,collapse="_"),namelist[s],sep="_")
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes1_',5*id,'.png',sep=""),figure_list,
         width=35+3,height=5*id,
         dpi=300,limitsize = F)
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes1_',5*id,'.pdf',sep=""),figure_list,
         width=35,height=5*id,
         dpi=300,limitsize = F)
  # gene_express=rbind(gene_express,tmp_gene_express)
}
levels=c('C1.0','C1.2','C1.5','C14.0','C14.1','C12')
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/monocle_20220907/C1_14_12'
namelist=list.files(indir,pattern='*txt$')
s=1
# gene_express=c()
for(s in seq(1,length(namelist))){
  print(s)
  tmp_gene_express=read.delim(file.path(getwd(),namelist[s]),sep="\t",as.is=T,header=T,check.names=F)
  # if(namelist[s]=='vlnplot20220824_sheet1V1.txt'){
  #   tmp_gene_express$new=paste(tmp_gene_express[,1],tmp_gene_express[,2],sep="|")
  #   tmp_gene_express=tmp_gene_express[,c(1,3)]
  # }
  # if(ncol(tmp_gene_express)>2){
  #   tmp_gene_express=tmp_gene_express[,c(1,2)]
  # }
  tmp_gene_express=tmp_gene_express[!duplicated(tmp_gene_express[,1]),]
  print(paste('gene number is',dim(tmp_gene_express)[1],'and total gene in ',sum(tmp_gene_express[,1]%in%expr_plotdata$gene_short_name)))
  colnames(tmp_gene_express)=c('gene','geneID')
  rownames(tmp_gene_express)=tmp_gene_express[,1]
  gene_express=tmp_gene_express
  tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%gene_express[,1],]
  tmp_figure$clusters=factor(tmp_figure$clusters,levels=levels)
  c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
  figure_list=wrap_plots(c1_linelist,ncol=5)
  id=floor(dim(gene_express)[1]/5)
  if(id*5<dim(gene_express)[1]){
    id=id+1
  }
  cluster_name=paste(paste(levels,collapse="_"),namelist[s],sep="_")
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes1_',5*id,'.png',sep=""),figure_list,
         width=35+3,height=5*id,
         dpi=300,limitsize = F)
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes1_',5*id,'.pdf',sep=""),figure_list,
         width=35,height=5*id,
         dpi=300,limitsize = F)
  # gene_express=rbind(gene_express,tmp_gene_express)
}
# 折线图的横坐标顺序：C1.0;  C1.3; C1.4;  C11.0, C11.1; C8
levels=c('C1.0','C1.3','C1.4','C11.0','C11.1','C8')
indir='/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/monocle_20220907/C1.0_1.3_1.4_11.0_11.1'
namelist=list.files(indir,pattern='*txt$')
s=1
# gene_express=c()
for(s in seq(1,length(namelist))){
  print(s)
  tmp_gene_express=read.delim(file.path(getwd(),namelist[s]),sep="\t",as.is=T,header=T,check.names=F)
  # if(namelist[s]=='vlnplot20220824_sheet1V1.txt'){
  #   tmp_gene_express$new=paste(tmp_gene_express[,1],tmp_gene_express[,2],sep="|")
  #   tmp_gene_express=tmp_gene_express[,c(1,3)]
  # }
  # if(ncol(tmp_gene_express)>2){
  #   tmp_gene_express=tmp_gene_express[,c(1,2)]
  # }
  tmp_gene_express=tmp_gene_express[!duplicated(tmp_gene_express[,1]),]
  print(paste('gene number is',dim(tmp_gene_express)[1],'and total gene in ',sum(tmp_gene_express[,1]%in%expr_plotdata$gene_short_name)))
  colnames(tmp_gene_express)=c('gene','geneID')
  rownames(tmp_gene_express)=tmp_gene_express[,1]
  gene_express=tmp_gene_express
  tmp_figure=expr_plotdata[expr_plotdata$gene_short_name%in%gene_express[,1],]
  tmp_figure$clusters=factor(tmp_figure$clusters,levels=levels)
  c1_linelist=geneline(expr_plotdata,gene_express,cluster_levels = levels)
  figure_list=wrap_plots(c1_linelist,ncol=5)
  id=floor(dim(gene_express)[1]/5)
  if(id*5<dim(gene_express)[1]){
    id=id+1
  }
  cluster_name=paste(paste(levels,collapse="_"),namelist[s],sep="_")
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes1_',5*id,'.png',sep=""),figure_list,
         width=35+3,height=5*id,
         dpi=300,limitsize = F)
  ggsave(paste(cluster_name,'_',gsub('-','',Sys.Date()),'.genes1_',5*id,'.pdf',sep=""),figure_list,
         width=35,height=5*id,
         dpi=300,limitsize = F)
  # gene_express=rbind(gene_express,tmp_gene_express)
}

#############调整图例20220909
#这个拟时序的图例的顺序请调整一下， C1.0/PC;  C1.2/PPh-a;  C1.5/PPh-b; C14.0/PCL, C14.1/DPh,C12/Pf;

p1=plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime", 
                        cell_size = 0.5) + theme(axis.title = element_text(face = "bold", 
                                                                           size = 15, colour = "black"), legend.position = "right", 
                                                 legend.text = element_text(face = "bold", size = 13, 
                                                                            colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                           size = 15, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                    size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                           size = 15), title = element_text(size = 15), strip.text = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                                  size = 15))
pdf(paste(name, "_subset_monocle_Pseudotime.pdf", sep = ""), 
    width = 8, height = 6, onefile = F)
print(p1)
dev.off()
png(paste(name, "_subset_monocle_Pseudotime.png", sep = ""), 
    width = 800, height = 600)
print(p1)
dev.off()
p1=plot_cell_trajectory(HSMM_myo, color_by = "clusters", 
                        show_branch_points = F, cell_size = 0.5) + theme(axis.title = element_text(face = "bold", 
                                                                                                   size = 15, colour = "black"), legend.position = "right", 
                                                                         legend.text = element_text(face = "bold", size = 13, 
                                                                                                    colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                   size = 15, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                            size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                   size = 15), title = element_text(size = 15), strip.text = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                                                          size = 15))
pdf(paste(name, "_subset_monocle_cluster.pdf", sep = ""), 
    width = 8, height = 6, onefile = F)
print(p1)
dev.off()
png(paste(name, "_subset_monocle_cluster.png", sep = ""), 
    width = 800, height = 600)
print(p1)
dev.off()
p1=plot_cell_trajectory(HSMM_myo, color_by = "clusters", 
                        cell_size = 0.5) + facet_wrap(~clusters, nrow = 4) + 
  theme(axis.title = element_text(face = "bold", size = 15, 
                                  colour = "black"), legend.position = "right", legend.text = element_text(face = "bold", 
                                                                                                           size = 13, colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                                     size = 15, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                                              size = 10), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                                     size = 10), title = element_text(size = 15), strip.text = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                                                                            size = 15))

pdf(paste(name, "_subset_monocle_cluster_split.pdf", sep = ""), 
    width = 10, height = 7, onefile = F)
print(p1)
dev.off()
png(paste(name, "_subset_monocle_cluster_split.png", sep = ""), 
    width = 1000, height = 700)
print(p1)
dev.off()
p1=plot_complex_cell_trajectory(HSMM_myo, color_by = "clusters", 
                                show_branch_points = F, cell_size = 0.5) + theme(axis.title = element_text(face = "bold", 
                                                                                                           size = 15, colour = "black"), legend.position = "right", 
                                                                                 legend.text = element_text(face = "bold", size = 13, 
                                                                                                            colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                           size = 15, colour = "black"), strip.text = element_text(face = "bold", 
                                                                                                                                                                                                                   size = 15))

pdf(paste(name, "_subset_monocle_cluster2.pdf", sep = ""), 
    width = 8, height = 6, onefile = F)
print(p1)
dev.off()
png(paste(name, "_subset_monocle_cluster2.png", sep = ""), 
    width = 800, height = 600)
print(p1)
dev.off()