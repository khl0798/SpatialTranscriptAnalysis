 
 args=commandArgs(T)
 file=args[1]
 out=args[2]
 library("gplots")
 library("RColorBrewer") 
 library('pheatmap')
 #path='E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/markergene_20210506/'
 # path="E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20210524/group4_group5_group6_20210528/group6"
 #path='E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20210621/字号变大/字号变大/17 cluster_GO_0506'
 path='E:/项目统计(正在做)/空间转录组/BHT201004/BHT201004_20210412售后分析/20211109_气泡图/heatmap'
 setwd(path)
 sample3 <- read.table(file.path(path,'C11_GO_heatmap.txt'), sep = "\t", row.names= 1,header = T,check.names = F) 
 rownames(sample3)=paste(rownames(sample3),sample3[,'GO_term'],sep=":")
 sample3=sample3[1:nrow(sample3),2:ncol(sample3)]
 sample3[sample3>0.05] <- NA
 sample3[is.na(sample3)]=0
 sample3[sample3<=0.05 & sample3>0.001] <- 2
 sample3[sample3<=0.001 & sample3>0.00001] <- 3
 sample3[sample3<=0.00001 & sample3>0.0000001] <- 4
 sample3[sample3<=0.0000001 & sample3>0.000000001] <- 5
 sample3[sample3<=0.000000001 & sample3>0] <- 6
 sample3[sample3==0] <- 1
 name=rownames(sample3)
 # name[26]='GO:0009834:Plant-type secondary cell wall biogenesis'
 # name[27]='GO:0009833:Plant-type primary cell wall biogenesis'
 # name[29]='GO:0009734:Auxin-activated signaling pathway'
 
#  for (i in 1:nrow(sample3)){
#   if(nchar(name[i])>50){
#     name[i]<-paste(substr(name[i],1,50),' etc...',sep='') 
#   }
# }
rownames(sample3)=name
 # pdf('top40pvalue_go_enrichment.pdf',width=12,height=8,onefile=F)
out='C11go_20211109'  # width=1500,height=800
# colnames(sample3)=paste('C',colnames(sample3),sep="")
png(file.path(path,paste(out, "_enrichment_heatmap_cluster_GOterm_20211109.png", sep = "")), 
    width = 2500, height = 1600)
pheatmap(sample3 ,kmeans_k=NA,
         color=c("white",colorRampPalette(c("#7B68EE","#CD3333"))(5)), 
         scale="none",breaks=NA, border_color="black",
         legend=TRUE, drop_levels = FALSE,cellwidth =45,cellheight=21,
         legend_breaks = 1:6, 
         legend_labels = c("","1e-03<p<=0.05", "1e-05<p<=1e-03", 
                           "1e-07<p<=1e-05", "1e-09<p<=1e-07", "p<=1e-09"),
         show_rownames=T,fontsize_row=18,fontsize_col=18,
         show_colnames=T, main=NA, fontsize=16,
         cluster_rows = T, cluster_cols = T 
)
dev.off()

pdf(file.path(path,paste(out, "_enrichment_heatmap_cluster_GOterm_20211109.pdf", sep = "")), 
    width = 25, height = 16)
pheatmap(sample3,kmeans_k=NA,
         color=c("white",colorRampPalette(c("#7B68EE","#CD3333"))(5)), 
         scale="none",breaks=NA, border_color="black",
         legend=TRUE, drop_levels = FALSE,cellwidth =45,cellheight=21,
         legend_breaks = 1:6, 
         legend_labels = c("","1e-03<p<=0.05", "1e-05<p<=1e-03", 
                           "1e-07<p<=1e-05", "1e-09<p<=1e-07", "p<=1e-09"),
         show_rownames=T,fontsize_row=18,fontsize_col=21,
         show_colnames=T, main=NA, fontsize=15,
         cluster_rows = T, cluster_cols = T 
)
dev.off()  
