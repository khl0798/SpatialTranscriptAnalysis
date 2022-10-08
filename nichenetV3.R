library(nichenetr)
library(Seurat)
library(tidyverse)
ligand_target_matrix = readRDS("/data7/konghl/data/SpatialTranscriptome/software/nichenet/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5] #
lr_network = readRDS("/data7/konghl/data/SpatialTranscriptome/software/nichenet/lr_network.rds")
head(lr_network)
## # A tibble: 6 x 4
##   from  to    source         database
##   <chr> <chr> <chr>          <chr>   
## 1 CXCL1 CXCR2 kegg_cytokines kegg    
## 2 CXCL2 CXCR2 kegg_cytokines kegg    
weighted_networks = readRDS("/data7/konghl/data/SpatialTranscriptome/software/nichenet/weighted_networks.rds")
head(weighted_networks$lr_sig)
####seuratObj = readRDS("seuratObj.rds")
seuratObj=readRDS('/data7/konghl/data/SpatialTranscriptome/E20210156/analysis_X3_X4_lognorm/combin.data.RDS')
seuratObj$celltype=paste('C',Idents(seuratObj),sep="")
seuratObj@meta.data %>% head()
seuratObj@meta.data$celltype %>% table() # note that the number of cells of some cell types is very low and should preferably be higher for a real application
# indicated cell types should be cell class identities
# check via: 
# seuratObj %>% Idents() %>% table()
#####包括以下3个函数nichenet_seuratobj_aggregate
######nichenet_seuratobj_cluster_de      nichenet_seuratobj_aggregate_cluster_de
seuratObj=subset(seuratObj,orig.ident=='X4-19S58966-002')
# nichenet_output = nichenet_seuratobj_aggregate(
#   seurat_obj = seuratObj, 
#   receiver = c("C2","C3","C4"), 
#   condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", 
#   sender = c('C5'), 
#   ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
nichenet_output = nichenet_seuratobj_cluster_de(
  seurat_obj = seuratObj, 
  receiver_reference = c('C1'), receiver_affected = c('C2','C3','C4'), 
  sender = c("C5",'C6','C7','C8'), assay_oi='Spatial',
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, 
  weighted_networks = weighted_networks, organism = "human",cutoff_visualization=0)
## [1] "Read in and process NicheNet's networks"
## [
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"

########结果展示
# ligand activity heatmap  配体活性矩阵
library(cowplot)
p1=nichenet_output$ligand_activity_target_heatmap
ggsave('ligand_activity_target_heatmap.png',p1,width=20,height=12)
p2=nichenet_output$ligand_receptor_heatmap
ggsave('ligand_receptor_heatmap.png',p2,width=10,height=8)
p2=nichenet_output$ligand_expression_dotplot
ggsave('ligand_expression_dotplot.png',p2,width=10,height=8)
#######################target基因的的表达值，热图
expression_obj=subset(seuratObj,idents=c('C1','C2','C3','C4','C5','C6','C7','C8'))
expression=t(as.matrix(expression_obj[['Spatial']]@data))
sample_info=data.frame(cell=colnames(expression_obj),tumor=Idents(expression_obj))
geneset_oi=nichenet_output$geneset_oi

expression_df_target = expression[,geneset_oi] %>% data.frame() %>% rownames_to_column("cell") %>% as_tibble() %>% inner_join(sample_info %>% select(cell,tumor), by =  "cell") 

aggregated_expression_target = expression_df_target %>% group_by(tumor) %>% select(-cell) %>% summarise_all(mean)

aggregated_expression_df_target = aggregated_expression_target %>% select(-tumor) %>% t() %>% magrittr::set_colnames(aggregated_expression_target$tumor) %>% data.frame() %>% rownames_to_column("target") %>% as_tibble() 

aggregated_expression_matrix_target = aggregated_expression_df_target %>% select(-target) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_target$target)

vis_target_tumor_expression_scaled = aggregated_expression_matrix_target %>% t() %>% scale_quantile() 
###%>% .[order_tumors,order_targets]
library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
order_tumors=paste('C',1:4,sep="")
order_targets=geneset_oi
vis_target_tumor_expression_scaled=vis_target_tumor_expression_scaled[order_tumors,]

p_target_tumor_scaled_expression = vis_target_tumor_expression_scaled  %>% make_threecolor_heatmap_ggplot("Tumor","Target", low_color = color[1],mid_color = color[50], mid = 0.5, high_color = color[100], legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))
hong ################################配体基因的表达值，热图
top_ligands=nichenet_output$top_ligands
expression_df_target = expression[,top_ligands] %>% data.frame() %>% rownames_to_column("cell") %>% as_tibble() %>% inner_join(sample_info %>% select(cell,tumor), by =  "cell") 

aggregated_expression_target = expression_df_target %>% group_by(tumor) %>% select(-cell) %>% summarise_all(mean)

aggregated_expression_df_target = aggregated_expression_target %>% select(-tumor) %>% t() %>% magrittr::set_colnames(aggregated_expression_target$tumor) %>% data.frame() %>% rownames_to_column("target") %>% as_tibble() 

aggregated_expression_matrix_target = aggregated_expression_df_target %>% select(-target) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_target$target)

vis_target_tumor_expression_scaled = aggregated_expression_matrix_target %>% t() %>% scale_quantile() 
###%>% .[order_tumors,order_targets]
library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
order_tumors=paste('C',5:8,sep="")
order_targets=geneset_oi
vis_target_tumor_expression_scaled=vis_target_tumor_expression_scaled[order_tumors,]

p_target_tumor_scaled_expression = vis_target_tumor_expression_scaled  %>% make_threecolor_heatmap_ggplot("Tumor","ligand", low_color = color[1],mid_color = color[50], mid = 0.5, high_color = color[100], legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))
ggsave('ligand_heatmap.png',p_target_tumor_scaled_expression,width=18,height=7)

##############################配体target基因的信号通路
##为了确定配体和目标靶标之间的信号通路，研究了哪些转录因子最能调节靶基因，并且最靠近配体的下游
##基于集成配体信号传导和基因调控网络中边缘的权重
##确定转录因子和感兴趣的配体之间的最短路径

ligand_tf_matrix =readRDS('/data7/konghl/data/SpatialTranscriptome/software/nichenet/ligand_tf_matrix.rds')
sig_network = readRDS("/data7/konghl/data/SpatialTranscriptome/software/nichenet/signaling_network.rds")
ligands_all = top_ligands[1]   #"TGFB1" # this can be a list of multiple ligands if required
targets_all = nichenet_output$top_targets[1:2]

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, 
                                                     targets_all = targets_all, weighted_networks = weighted_networks)
#####为了更直观的可视化，直接的比较网络的边的权重，对权重直径了标准化
# For better visualization of edge weigths: normalize edge weights to make them comparable 
###between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, 
                          ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

# To render the graph: uncomment following line of code
# library(DiagrammeR)
# library(ggplot2)
# setwd('E:/项目统计(正在做)/空间转录组/E20210156-2/lognormalize_20220423/nichenet/')
# graph_min_max=readRDS('E:/项目统计(正在做)/空间转录组/E20210156-2/lognormalize_20220423/nichenet/graph_min_max.RDS')
# p1=render_graph(graph_min_max, layout = "tree",output='graph')
# 
# 
# 
# png(file='network_ligand_target.png')
# render_graph(graph_min_max, layout = "tree")
# dev.off()
# 
# ggsave('network_ligand_target.png',p1,width=10,height=8)
############把结果导出到cytoscape画图用
output_path = '/data7/konghl/data/SpatialTranscriptome/E20210156/analysis_X3_X4_lognorm/20220426/nichenet/'
write_output = TRUE # change to TRUE for writing output

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>%
              mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}

# networks with information of supporting data sources ('import network' in Cytoscape)
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)

if(write_output){
  data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% 
                    unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
  bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% 
    write_tsv(paste0(output_path,"annotation_table.txt"))
}

####现在，我们将查看哪些收集的数据源支持此网络中的交互。

####data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
#####head(data_source_network) 

######circos图可视化配体-target基因和配体受体调控
#####先根据最强表达的细胞类型预测活性配体，即哪几种细胞类型种，配体活性表达是最高的,
#####我们需要确定每种细胞类型，它们表达的配体比其他细胞类型更强。计算发送细胞中平均配体表达量
avg_expression_ligands = AverageExpression(seuratObj, features = nichenet_output$top_ligands)
#####分配配体给发送细胞
######为了给发送端细胞类型分配配体，我们可以查找哪个发送端细胞类型的表达式高于平均值+ SD。
sender_ligand_assignment = avg_expression_ligands$Spatial[,c('C5','C6','C7','C8'),drop=FALSE] %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})

all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
####如果是general_ligands为空值怎么办
general_ligands = nichenet_output$top_ligands %>% setdiff(unique_ligands)

B_specific_ligands = sender_ligand_assignment$C5 %>% names() %>% setdiff(general_ligands)
NK_specific_ligands = sender_ligand_assignment$C6 %>% names() %>% setdiff(general_ligands)
Mono_specific_ligands = sender_ligand_assignment$C7 %>% names() %>% setdiff(general_ligands)
DC_specific_ligands = sender_ligand_assignment$C8 %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("C5", times = B_specific_ligands %>% length()),
                  rep("C6", times = NK_specific_ligands %>% length()),
                  rep("C7", times = Mono_specific_ligands %>% length()),
                  rep("C8", times = DC_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(B_specific_ligands, NK_specific_ligands, Mono_specific_ligands, DC_specific_ligands, general_ligands))

#ligand_type_indication_df %>% head

#定义感兴趣的配体-目标链接
#为了避免circos图中有太多配体目标链接，我们将只显示权重高于预定义截止值的链接:
##属于最低分数的40%的链接被删除。这并不是说用于这种可视化的边界和其他边界可以根据用户的需要进行更改。
active_ligand_target_links_df = nichenet_output$ligand_target_df %>% mutate(target_type = "test") %>% inner_join(ligand_type_indication_df) 
# if you want ot make circos plots for multiple gene sets,
#combine the different data frames and differentiate which target belongs to which gene set via the target type

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
#################准备circos可视化:给每个片段配体和目标特定的颜色和顺序
grid_col_ligand =c("General" = "lawngreen",
                   "C5" = "royalblue",
                   "C6" = "darkgreen",
                   "C7" = "violet",
                   "C8" = "steelblue2")
grid_col_target =c(
  "test" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

###准备可视化的circos:排序配体和目标
target_order = circos_links$target %>% unique()
ligand_order = c(Mono_specific_ligands, DC_specific_ligands, NK_specific_ligands,B_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

##准备circos可视化:定义不同片段之间的间隙
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "C5") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "C6") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "C7") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "C8") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "test") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

###渲染circos的情节(所有链接相同的透明度)。只有表明每个靶基因的阻滞的宽度与配体-靶的调控电位成正比
####(~支持调控相互作用的先验知识)。
library(circlize)

png('circos.png')
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
dev.off()
circos.clear()

## 绘制circos图(透明度由配体-靶标相互作用的调控潜力决定) 透明度图
png('circos.trans.png')
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,
             transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
dev.off()
circos.clear()


svg("ligand_target_circos.svg", width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()
## 在circos图中可视化优先配体与受体的相互作用
lr_network_top_df = nichenet_output$ligand_receptor_df %>% mutate(receptor_type = "receptor") %>% 
  inner_join(ligand_type_indication_df)
grid_col_ligand =c("General" = "lawngreen",
                   "C5" = "royalblue",
                   "C6" = "darkgreen",
                   "C7" = "violet",
                   "C8" = "steelblue2")
grid_col_receptor =c(
  "receptor" = "darkred")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
### 制备可视化的circos:有序配体和受体

receptor_order = circos_links$receptor %>% unique()
ligand_order = c(Mono_specific_ligands, DC_specific_ligands, NK_specific_ligands,B_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)
## 准备马戏团可视化:定义不同片段之间的间隙

width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "C5") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "C6") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "C7") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "C8") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_receptor,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "receptor") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
)
## 渲染马戏团的情节(所有链接相同的透明度)。只有表明每个受体的阻滞的宽度与配体-受体相互作用的重量成比例(~支持相互作用的先验知识)。

png('circos_receptor.png')
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
dev.off()
circos.clear()
## 渲染circos图(透明程度由配体-受体相互作用的先验相互作用权重决定——正如指示每个受体的块的宽度)

png('circos_receptor.trans.png')
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()
dev.off()














# p=plot_grid(p1,p2)
ligand_activities=nichenet_output$ligand_activities
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
rotated_dotplot = DotPlot(seuratObj %>% subset(celltype %in% sender_celltypes), features = order_ligands_adapted %>% rev(), cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot



p1=DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
library(ggplot2)
ggsave('top_ligand.png',p1,width=8,height=6)
p1=nichenet_output$ligand_target_heatmap
ggsave('ligand_target_pheatmap.png',p1,width=8,height=6)
####或者是修改颜色
p1=nichenet_output$ligand_target_heatmap + scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + xlab("anti-LCMV response genes in CD8 T cells") + ylab("Prioritized immmune cell ligands")
ggsave('ligand_target_pheatmapV1.png',p1,width=8,height=6)
####提取矩阵
p1=nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
#####配体活性和配体及靶基因的热图
p1=nichenet_output$ligand_activity_target_heatmap
ggsave('ligand_activities_ligand_target_heatmap.png',p1,width=15,height=10)
######顶级配体的受体细胞的热图
p1=nichenet_output$ligand_receptor_heatmap
ggsave('ligand_receptor_heatmap.png',p1,width=15,height=10)


