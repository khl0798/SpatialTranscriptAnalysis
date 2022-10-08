library(infercnv)
library(ggplot2)
library(pheatmap)
species='human'
name='X4'
out_dir = getwd()
gene_order_file='/data3/denghj/app/Miniconda3/envs/R3.6/lib/R/library/singlercell2/genes/gene_order_human.txt'
expr=readRDS('/data7/konghl/data/SpatialTranscriptome/E20210156/analysis_X3_X4_lognorm/combin.data.RDS')
expr=subset(expr,orig.ident=='X4-19S58966-002') #X4-19S58966-002 'X3-18S08906-002'
alllist=c(1:6,10)
reflist=c(10,6)
expr <- subset(expr, idents = alllist)
count_matrix <- expr@assays$Spatial@counts
reflist2 <- paste("C", reflist, sep = "")
cell_metadata <- data.frame(clusters = paste("C", Idents(expr), sep = ""))
cell_metadata[, 1] <- as.character(cell_metadata[, 1])
rownames(cell_metadata) <- colnames(expr)
saveRDS(cell_metadata,file.path(out_dir,paste(name,'_cell_metadata.RDS',sep="")))
saveRDS(count_matrix,file.path(out_dir,paste(name,'_count_matrix.RDS',sep="")))
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = count_matrix, 
                                    annotations_file = cell_metadata, delim = "\t", gene_order_file = gene_order_file, 
                                    ref_group_names = reflist2)
infercnv_obj = infercnv::run(infercnv_obj, cutoff = 0.1, output_format = "pdf",
                             out_dir = out_dir, cluster_by_groups = TRUE, plot_steps = FALSE, 
                             denoise = TRUE, HMM = FALSE, no_prelim_plot = TRUE, png_res = 60)