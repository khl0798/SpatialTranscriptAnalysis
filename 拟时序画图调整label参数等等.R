library(monocle)
#########主要是拟时序修改了部分画图的参数，最终的字体版本；后面如果文章需要修改拟时序图的时候可以参考字体调整。
HSMM_myo1 <- setOrderingFilter(HSMM_myo, genes355[,1])
HSMM_myo1 <- reduceDimension(HSMM_myo1, max_components = 5, method = "DDRTree")
HSMM_myo1 <- orderCells(HSMM_myo1,reverse=T)
cell_Pseudotime <- data.frame(pData(HSMM_myo1)$Pseudotime)
rownames(cell_Pseudotime) <- rownames(cell_metadata)
write.table(cell_Pseudotime, file = paste(name, "_subset_cell_Pseudotime.xls", sep = ""), sep = "\t", row.names = T, col.names = F)
Dims = t(monocle::reducedDimS(HSMM_myo1))
colnames(Dims) = c("Component 1", "Component 2")
write.table(Dims, file = paste(name, "_subset_cell_reducedDimS.xls", 
                               sep = ""), sep = "\t", row.names = T, col.names = NA)
save(HSMM_myo1, file = paste(name, "_HSMM_myo.Rdata", sep = ""))
Dims = t(monocle::reducedDimS(HSMM_myo1))[,1:2]
colnames(Dims) = c("Component 1", "Component 2")
write.table(Dims, file = paste(name, "_subset_cell_reducedDimS.xls", 
                               sep = ""), sep = "\t", row.names = T, col.names = NA)
##### pseudotime 图
name='test20220322'
pdf(paste(name, "_subset_monocle_Pseudotime.pdf", sep = ""), 
    width = 10, height = 6, onefile = F)
print(plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "Pseudotime", 
                           cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                            size = 15, colour = "black"), legend.position = "right", 
                                                  legend.text = element_text(face = "bold", size = 18, 
                                                                             colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                            size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                     size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                            size = 15), title = element_text(size = 15), 
                                                  strip.text = element_text(face = "bold", size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=2,label.vjust=2))+ ###legend的图例调大，legend label字体调大
        guides(shape = guide_legend(override.aes = list(size = 5))))

dev.off()
png(paste(name, "_subset_monocle_Pseudotime.png", sep = ""), width = 1000, height = 600)
print(plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "Pseudotime", 
                           cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                            size = 15, colour = "black"), legend.position = "right", 
                                                  legend.text = element_text(face = "bold", size = 18, 
                                                                             colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                            size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                     size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                            size = 15), title = element_text(size = 15), 
                                                  strip.text = element_text(face = "bold", size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=2,label.vjust=2))+
        guides(shape = guide_legend(override.aes = list(size = 5))))

dev.off()
###### cluster split 图
pdf(paste(name, "_subset_monocle_cluster2.pdf", sep = ""), width = 11, height = 7, onefile = F)
print(plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "clusters", 
                           cell_size = 0.5) + facet_wrap(~clusters, nrow = 4) + 
        theme(axis.title = element_text(face = "bold", size = 15, 
                                        colour = "black"), legend.position = "right", legend.text = element_text(face = "bold", 
                                                                                                                 size = 18, colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                                           size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                                                    size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                                           size = 15), title = element_text(size = 15), 
              strip.text = element_text(face = "bold", size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
        guides(color = guide_legend(override.aes = list(size = 5))))
# scale_color_manual(values=color_values))
dev.off()
HSMM_myo1 <- reduceDimension(HSMM_myo1, max_components = 6, method = "DDRTree")
HSMM_myo1 <- orderCells(HSMM_myo1,reverse=T)
png(paste(name, "_subset_monocle_cluster2.png", sep = ""), width = 1100, height = 700)
print(plot_cell_trajectory(HSMM_myo1, x=1,y=2,color_by = "clusters", 
                           cell_size = 0.5) + facet_wrap(~clusters, nrow = 4) + 
        theme(axis.title = element_text(face = "bold", size = 15, 
                                        colour = "black"), legend.position = "right", legend.text = element_text(face = "bold", 
                                                                                                                 size = 18, colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                                           size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                                                    size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                                           size = 15), title = element_text(size = 15), 
              strip.text = element_text(face = "bold", size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
        guides(color = guide_legend(override.aes = list(size = 5))))
# scale_color_manual(values=color_values))
dev.off()                                                                                                                                                                                                                       
#### cluster图
png(paste(name, "_subset_monocle_cluster.png", sep = ""), width = 1100, height = 700)
print(plot_cell_trajectory(HSMM_myo1,x=1,y=2, color_by = "clusters", 
                           show_branch_points = F, cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                                                    size = 15, colour = "black"), legend.position = "right", 
                                                                          legend.text = element_text(face = "bold", size = 18, 
                                                                                                     colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                    size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                             size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                    size = 15), title = element_text(size = 15), 
                                                                          strip.text = element_text(face = "bold",size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
        guides(color = guide_legend(override.aes = list(size = 5))))
# scale_color_manual(values=color_values))
dev.off()
pdf(paste(name, "_subset_monocle_cluster.pdf", sep = ""), width = 11, height = 7,onefile = F)
print(plot_cell_trajectory(HSMM_myo1,x=1,y=2, color_by = "clusters", 
                           show_branch_points = F, cell_size = 1) + theme(axis.title = element_text(face = "bold", 
                                                                                                    size = 15, colour = "black"), legend.position = "right", 
                                                                          legend.text = element_text(face = "bold", size = 18, 
                                                                                                     colour = "black"), legend.title = element_text(face = "bold", 
                                                                                                                                                    size = 20, colour = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                                             size = 15), axis.text.y = element_text(face = "bold", 
                                                                                                                                                                                                                                                    size = 15), title = element_text(size = 15), 
                                                                          strip.text = element_text(face = "bold",size = 15),legend.key.height = unit(2,"cm"))+
        guides(fill=guide_legend(label.hjust=3,label.vjust=3))+
        guides(color = guide_legend(override.aes = list(size = 5))))
# scale_color_manual(values=color_values))
dev.off()

diff_test_res2 <- differentialGeneTest(HSMM_myo1[ordering_genes1, ], 
                                       fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 10)
sig_gene_names <- subset(diff_test_res2, qval < 0.05)[, "gene_short_name"]
sig_gene_names = sig_gene_names[!duplicated(sig_gene_names)]
heatmap=plot_genes_branched_heatmap(HSMM_myo[sig_gene_names,],branch_point = 2,cluster_rows = T,cores=12,return_heatmap = T)
anno=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/anno_tf_go_kegg_ko.RDS')
row_cluster = cutree(heatmap$ph_res$tree_row,6)
row_cluster=data.frame(genes=names(row_cluster),genecluster=row_cluster,anno[names(row_cluster),])
write.table(row_cluster, file = paste(name, "_gene_clusters_branch1.xls", sep = ""), sep = "\t", row.names = T)
for(i in 1:6){
  tmp=row_cluster%>%filter(genecluster==i)
  print(head(tmp[,1:2]))
  write.table(tmp[,1:2],paste('genecluster',i,'.branch1.txt',sep=""),sep="\t",quote=F,row.names=F,col.names=F)
}
png(file='branch_point2_heatmap.png',width=800,height=600)
print(heatmap$ph_res)
#################### print(plot_genes_branched_heatmap(HSMM_myo[sig_gene_names,],cluster_rows = T,cores=12)) 此行代码应该暂时用不上
dev.off()
pdf(file='branch_point2_heatmap.pdf',width=8,height=6,onefile = F)
print(heatmap$ph_res)