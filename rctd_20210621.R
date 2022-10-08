expr=readRDS('/data7/guohua/Visium/BHT201002_40samples/analysis/combin.data.RDS')
myRCTD11=readRDS('/data7/guohua/10X/BHT201002_sc/SC_VI40_0310/myRCTD.RDS')
cell_type_info <- get_cell_type_info(myRCTD@reference@assays$RNA@counts, 
                                     myRCTD@reference@meta.data$liger_ident_coarse, 
                                     myRCTD@reference@meta.data$nUMI)
# > cell_type_info[[2]]
# [1] "Epithelial cells"    "Endothelial cells"   "Fibroblasts"        
# [4] "Smooth Muscle cells" "NK & T cells"        "B cells"            
# [7] "Myeloid cells"       "Mast cells" 
puck=myRCTD@spatialRNA
# gene_list <- intersect(rownames(cell_type_info[[1]]),rownames(puck@counts))
# gene_list <- gene_list[rowSums(as.matrix(puck@counts[gene_list,])) >= 10]   # 1587
# gene_list <- gene_list[apply(cell_type_info[[1]][gene_list,],1,function(x) max(x)) >= 2e-5]   #1587
# puck <- restrict_puck(puck, names(which(puck@nUMI >= 100)))
# restrict_counts<-function (puck, gene_list, UMI_thresh = 1, UMI_max = 20000){
#   keep_loc = (puck@nUMI >= UMI_thresh) & (puck@nUMI <= UMI_max)
#   puck@counts = puck@counts[gene_list, keep_loc]
#   if (length(puck@cell_labels) > 0) 
#     puck@cell_labels = puck@cell_labels[keep_loc]
#   puck@nUMI = puck@nUMI[keep_loc]
#   return(puck)
# }
# puck <- restrict_counts(puck, gene_list, UMI_max = 200000)
# Epithelial_genes <- names(which(cell_type_info[[1]][gene_list,'Epithelial cells']*2 > 
#             apply(cell_type_info[[1]][gene_list,
#   !(cell_type_info[[2]] %in% c("Epithelial cells"))],1,max)))
# Epithelial_genes <- Epithelial_genes[cell_type_info[[1]][Epithelial_genes,"Epithelial cells"] >= 2e-5]
# 
# # puck_d <- get_decomposed_data(myRCTD@results$results_df, Epithelial_genes, puck, 
#                               myRCTD@results$weights_doublet, cell_type_info)
library(parallel)
library ("foreach") 
library ("doParallel")
library(Matrix)
cl <- makeCluster(16)
registerDoParallel(cl)
# results = foreach(x = c(1:100000),.combine = 'c') %dopar% square(x)

decompose_doublet_fast<-function (bead, weights, gene_list, cell_type_info, type1, type2){
  N_genes = length(gene_list)
  expect_1 = vector(mode = "numeric", length = N_genes)
  expect_2 = vector(mode = "numeric", length = N_genes)
  variance = vector(mode = "numeric", length = N_genes)
  names(expect_1) = gene_list
  names(expect_2) = gene_list
  names(variance) = gene_list
  epsilon = 1e-10
  for (ind in which(bead > 0)) {
    gene = gene_list[ind]
    denom = weights[1] * cell_type_info[[1]][gene, type1] + 
      weights[2] * cell_type_info[[1]][gene, type2] + epsilon
    posterior_1 = (weights[1] * cell_type_info[[1]][gene,type1] + epsilon/2)/denom
    expect_1[[ind]] = posterior_1 * bead[gene]
    expect_2[[ind]] = bead[gene] - posterior_1 * bead[gene]
    variance[[ind]] = posterior_1 * bead[gene] * (1 - posterior_1)
  }
  return(list(expect_1 = expect_1, expect_2 = expect_2, variance = variance))
}
results_df=myRCTD@results$results_df;
gene_list=Epithelial_genes; 
weights_doublet=myRCTD@results$weights_doublet; 

get_decomposed_data<-function (results_df, gene_list, puck, weights_doublet, cell_type_info) 
{
  doublets <- results_df[results_df$spot_class == "doublet_certain", ]
  first_DGE <- Matrix(0, nrow = dim(doublets)[1], ncol = length(gene_list))
  second_DGE <- Matrix(0, nrow = dim(doublets)[1], ncol = length(gene_list))
  rownames(first_DGE) = rownames(doublets)
  rownames(second_DGE) = rownames(doublets)
  colnames(first_DGE) = gene_list
  colnames(second_DGE) = gene_list
  # for (ind in 1:dim(doublets)[1]) {
  decompose_x=function(ind,doublets){
    barcode = rownames(doublets)[ind]
    doub_res <- decompose_doublet_fast(puck@counts[gene_list, 
                                                   barcode], weights_doublet[barcode, ], gene_list, 
                                       cell_type_info, results_df[barcode, "first_type"], 
                                       results_df[barcode, "second_type"])
    new_list=list()
    new_list[[barcode]]=doub_res
    return(new_list)
  }
  x1=foreach(ind=1:10000,.combine='c')%dopar% decompose_x(ind,doublets)
  x2=foreach(ind=10001:20000,.combine='c')%dopar% decompose_x(ind,doublets)
  x3=foreach(ind=20001:30000,.combine='c')%dopar% decompose_x(ind,doublets)
  x4=foreach(ind=30001:dim(doublets)[1],.combine='c')%dopar% decompose_x(ind,doublets)
  stopCluster(cl)
  ### 按照doublets行名对应的矩阵的信息,需要检查一下
  outs=c(x1,x2,x3,x4)
  sum(names(outs)==rownames(doublets)[1:dim(doublets)[1]])
  outs=data.frame(outs)
  first_DGE=t(outs[gene_list,seq(1,ncol(outs),3)])
  second_DGE=t(outs[gene_list,seq(2,ncol(outs),3)])
  rownames(first_DGE)=rownames(doublets)
  rownames(second_DGE)=rownames(doublets)
  
  singlet_id <- results_df$spot_class == "singlet"
  all_DGE <- rbind(first_DGE, second_DGE, t(as.matrix(puck@counts[gene_list, singlet_id])))
  rownames(all_DGE)=c(paste(rownames(doublets),doublets$first_type,sep="_"),
                      paste(rownames(doublets),doublets$second_type,sep="_"),
                      paste(rownames(results_df[singlet_id,]),results_df[singlet_id, "first_type"],sep="_"))
  # all_DGE_1=data.frame(all_DGE,'spot_class'=c(rep('doublet_certain',2*nrow(doublets)),
  #                                             rep('singlet',nrow(results_df[singlet_id,]))))
  
  cell_type_labels <- unlist(list(doublets$first_type, doublets$second_type, 
                                  results_df[singlet_id, "first_type"]))
  coords <- rbind(puck@coords[rownames(doublets), c("x", "y")], 
                  puck@coords[rownames(doublets), c("x", "y")], puck@coords[singlet_id, 
                                                                            c("x", "y")])
  nUMI <- c(weights_doublet[rownames(doublets), "first_type"] * 
              puck@nUMI[rownames(first_DGE)], weights_doublet[rownames(doublets), 
                                                              "second_type"] * puck@nUMI[rownames(second_DGE)], puck@nUMI[singlet_id])
  rownames(coords) = 1:dim(coords)[1]
  names(nUMI) = 1:dim(coords)[1]
  rownames(all_DGE) = 1:dim(coords)[1]
  SpatialRNA<-function (coords = NULL, counts, nUMI = NULL) 
  {
    if (is.null(coords)) {
      coords <- fake_coords(counts)
    }
    if (is.null(nUMI)) {
      nUMI = colSums(counts)
    }
    names(nUMI) <- colnames(counts)
    new("SpatialRNA", coords = coords, counts = counts, nUMI = nUMI)
  }
  ###dataan <- as(as.matrix(datan),"dgCMatrix") # 转换成dgCMatrix

  puck_d1 <- SpatialRNA(coords, Matrix(t(all_DGE),sparse=TRUE), nUMI)
  puck_d1@cell_labels <- cell_type_labels
  names(puck_d1@cell_labels) = 1:dim(coords)[1]
  puck_d1@cell_type_names <- cell_type_info[[2]]
  
  return(puck_d1)
}

my_mod <- function(p) {
  p + xlim(c(3800,5700)) + ylim(c(2100,3500)) + geom_segment(aes(x = 4000, y = 2400, xend = 4154, yend = 2400), color = "black")+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+ theme(legend.position="top")
}

#int_genes <- c("Rgs14","Spink8",'Cpne9','Kcnf1')
int_genes <- c("Rgs14",'Cpne9')
MULT <- 500
p1 <- list()
p2 <- list()
for(ind in 1:length(int_genes)) {
  gene <- int_genes[ind]
  max_val <- 20*mean((puck_d@counts[gene,]/puck_d@nUMI)[cell_barc])
  cur_range <- c(0, MULT*.00087)
  if(ind == 2)
    cur_range <- c(0,MULT*.00025)
  barcodes <- rownames(gene_df)
  raw_df <- puck_d@coords[barcodes,]
  raw_df$val <- as.integer((puck_d@counts[gene,]/puck_d@nUMI)[barcodes] > 0.0002)
  raw_df <- raw_df[raw_df$val > 0,]
  raw_df$val <- factor(raw_df$val)
  gene_df$gene <- puck_d@counts[gene,cell_barc]/puck_d@nUMI[cell_barc]
  fit <- loess(gene ~ x*y, span = 0.8, data = gene_df, degree = 1)
  gene_df$fitted <- fitted(fit)
  plot_val <- pmax(0,gene_df$fitted)
  names(plot_val) <- as.character(which(cell_barc))
  min_val = MULT*min(plot_val);
  max_val = MULT*quantile(plot_val[gene_df$x  >= 3800 & gene_df$x < 5700 & gene_df$y > 2100 & gene_df$y < 3500],0.95)
  p2[[ind]]<- plot_puck_continuous(puck_d,barcodes,MULT*plot_val[names(puck_d@cell_labels)], ylimit = c(min_val,max_val), size = 1.25, alpha = 0.2)
  cur_range <- c(-.00001*MULT, signif(max_val,2))
  p2[[ind]] <- my_mod(p2[[ind]])+ ggplot2::scale_colour_gradientn(paste(gene, "Smoothed"), colors = pals::brewer.orrd(20)[2:20],limits = cur_range, breaks = c(0,signif(max_val,2)), labels = c(0,signif(max_val,2)))
  p2[[ind]] <- p2[[ind]] + geom_point(data=raw_df, aes(fill=val),,size = 0.1) + scale_fill_manual("Raw",values = c("#000000"),labels = c("")) + guides(fill = guide_legend(override.aes = list(size=2)))
}
#ggarrange(p1[[1]], p2[[1]], p1[[2]], p2[[2]],p1[[3]], p2[[3]],p1[[4]], p2[[4]], nrow = 4,ncol=2)
#ggarrange(p1[[1]], p2[[1]],p1[[2]], p2[[2]], nrow = 2,ncol=2)
ggarrange(p2[[1]], p2[[2]],nrow = 1,ncol=2)
my_table=data.frame(sample(1:1000,100),sample(2:2000,100),rnorm(100))
colnames(my_table)=c('x','y','value')
plot <- ggplot2::ggplot(my_table, ggplot2::aes(x = x, y = y))
plot <- plot + ggplot2::geom_point(ggplot2::aes(size = 3, 
                    shape = 16, color = value, stroke = 0), alpha = 0.5)
my_pal = pals::brewer.ylorrd(20)[2:20]

sc <- ggplot2::scale_colour_gradientn(colors = my_pal)
plot <- plot + sc + ggplot2::scale_shape_identity() + 
  ggplot2::theme_classic() + 
  ggplot2::scale_size_identity()
plot <- plot + ggplot2::coord_fixed()
# + ggplot2::xlim(xlim) + 
#   ggplot2::ylim(ylim)
raw_df=data.frame(x=sample(0:3000,100),y=sample(1000:2000,100),val=0)
plot<-plot+geom_point(data=raw_df, aes(fill=val),size = 3)+
  guides(fill = guide_legend(override.aes = list(size=2)))+
    xlim(c(3800,5700)) + ylim(c(2100,3500)) +
  geom_segment(aes(x = 4000, y = 2400, xend = 4154, yend = 2400), color = "black")+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank())+ 
  theme(legend.position="top")
  
plot
