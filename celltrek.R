options(stringsAsFactors = F)
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")
brain_st_cortex <- readRDS("brain_st_cortex.rds")
brain_sc <- readRDS("brain_sc.rds")
## Visualize the ST data
p1=SpatialDimPlot(brain_st_cortex)
ggsave('st_idents.png',p1,width=8,height=6)
brain_traint <- CellTrek::traint(st_data=brain_st_cortex, 
                                 sc_data=brain_sc, sc_assay='RNA',
                                 cell_names='cell_type')
brain_traint=sc_st_int
brain_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=brain_sc, sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=30, ntree=1000, 
                                     dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)
# $celltrek
st_sc_int=brain_traint;int_assay='traint'; sc_data=brain_sc;sc_assay = 'RNA';
reduction='pca'; intp=T;intp_pnt=5000; intp_lin=F;nPCs=30; ntree=1000;
dist_thresh=0.55;top_spot=5;spot_n=5;repel_r=20;repel_iter=20; keep_model=T;
dist_res <- CellTrek:::celltrek_dist(st_sc_int=st_sc_int, int_assay=int_assay, reduction=reduction, intp=intp, intp_pnt=intp_pnt, intp_lin=intp_lin, nPCs=nPCs, ntree=ntree, keep_model=T)

brain_celltrek$cell_type <- factor(brain_celltrek$cell_type, levels=sort(unique(brain_celltrek$cell_type)))

p1=celltrek_vis(brain_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_type:id_new),
                       brain_celltrek@images$anterior1@image, brain_celltrek@images$anterior1@scale.factors$lowres)

brain_sgraph_KL <- CellTrek::scoloc(brain_celltrek_glut, col_cell='cell_type', use_method='KL', eps=1e-50)
brain_celltrek1111 <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=brain_sc, sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                     dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)
## D

brain_celltrek_l5 <- subset(brain_celltrek, subset=cell_type=='L5 IT')
brain_celltrek_l5@assays$RNA@scale.data <- matrix(NA, 1, 1)
brain_celltrek_l5$cluster <- gsub('L5 IT VISp ', '', brain_celltrek_l5$cluster)
brain_celltrek_l5 <- FindVariableFeatures(brain_celltrek_l5)
vst_df <- brain_celltrek_l5@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
nz_test <- apply(as.matrix(brain_celltrek_l5[['RNA']]@data), 1, function(x) mean(x!=0)*100)
hz_gene <- names(nz_test)[nz_test<20]
mt_gene <- grep('^Mt-', rownames(brain_celltrek_l5), value=T)
rp_gene <- grep('^Rpl|^Rps', rownames(brain_celltrek_l5), value=T)
vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(., -vst.variance.standardized)
feature_temp <- vst_df$id[1:2000]
brain_celltrek_l5_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=brain_celltrek_l5, assay='RNA', approach='cc', gene_select = feature_temp, sigm=140, avg_cor_min=.4, zero_cutoff=3, min_gen=40, max_gen=400)
