import os
IN_COLAB = "google.colab" in sys.modules

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

results_folder = './results/lymph_nodes_analysis/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'
# adata_vis = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
vis_filepath='/data7/konghl/data/SpatialTranscriptome/BHT211041/RNA/PC2021019/outs'
def read_visium(vis_path,library_id):
  adata_vis_tmp=sc.read_visium(vis_path,library_id=library_id)
  adata_vis_tmp.obs_names_make_unique()
  # adata_vis_tmp.var_names_make_unique()
  adata_vis_tmp.obs['sample'] = list(adata_vis_tmp.uns['spatial'].keys())[0]
  adata_vis_tmp.var['SYMBOL'] = adata_vis_tmp.var_names
  adata_vis_tmp.var_names = adata_vis_tmp.var['gene_ids']
  adata_vis_tmp.var_names.name = None
  adata_vis_tmp.obs_names=adata_vis_tmp.obs['sample']+'-'+adata_vis_tmp.obs_names
  return(adata_vis_tmp)
sample_name=['PC2021019','PC2021023','PC2021024','PC2021025']
slide_list=[]
for s in sample_name:
  adata_vis_test1=read_visium('/data7/konghl/data/SpatialTranscriptome/BHT211041/RNA/'+s+'/outs',library_id=s)
  slide_list.append(adata_vis_test1)
adata_vis = slide_list[0].concatenate(slide_list[1:], index_unique=None,batch_key='sample',
                uns_merge='unique',batch_categories=sample_name)


# adata_vis=sc.read_visium(vis_filepath)
# adata_vis.var_names_make_unique()

# adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]
# 
# rename genes to ENSEMBL
# adata_vis.var['SYMBOL'] = adata_vis.var_names
# adata_vis.var_names = adata_vis.var['gene_ids']
# adata_vis.var_names.name = None
##########
# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# ref_path='/data7/konghl/data/SpatialTranscriptome/BHT211041/singlecell_analysis/expr_new.h5ad'
ref_path='/data7/konghl/data/SpatialTranscriptome/BHT211041/cell2trek/expr_new_20220403.h5ad'
adata_ref = sc.read_h5ad(ref_path)

# # Use ENSEMBL as gene IDs to make sure IDs are unique and correctly matched
# adata_ref.var['SYMBOL'] = adata_ref.var.index
# adata_ref.var.index = adata_ref.var['GeneID-2'].copy()
# adata_ref.var_names = adata_ref.var['GeneID-2'].copy()
# adata_ref.var.index.name = None
# adata_ref.raw.var['SYMBOL'] = adata_ref.raw.var.index
# adata_ref.raw.var.index = adata_ref.raw.var['GeneID-2'].copy()
# adata_ref.raw.var.index.name = None

from cell2location.utils.filtering import filter_genes
###高可变基因

# selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_vis_plt = adata_vis.copy()

# Log-transform (log(data + 1))
sc.pp.log1p(adata_vis_plt)

hvg_2000=pd.read_csv('/data7/konghl/data/SpatialTranscriptome/BHT211041/cell2trek/hgv2000.csv', index_col=0,header=None)
intersect_ref_hvg2000 = np.intersect1d(adata_ref.var_names, hvg_2000.index)

# Find highly variable genes within each sample
# adata_ref_plt=adata_ref
# adata_ref_plt.var['highly_variable'] = False
# for s in adata_ref_plt.obs['sample'].unique():
# adata_ref_plt_1 = adata_ref_plt[adata_ref_plt.obs['sample'].isin([s]), :]
# sc.pp.highly_variable_genes(adata_ref_plt, min_mean=0.0125, max_mean=5, min_disp=0.5, n_top_genes=2000)
# hvg_list = list(adata_ref_plt_1.var_names[adata_ref_plt_1.var['highly_variable']])
# adata_vis_plt.var.loc[hvg_list, 'highly_variable'] = True
    
# Scale the data ( (data - mean) / sd )
sc.pp.scale(adata_ref_plt, max_value=10)
# PCA, KNN construction, UMAP
# sc.tl.pca(adata_vis_plt, svd_solver='arpack', n_comps=40, use_highly_variable=True)
# sc.pp.neighbors(adata_vis_plt, n_neighbors = 20, n_pcs = 40, metric='cosine')
# sc.tl.umap(adata_vis_plt, min_dist = 0.3, spread = 1)

# with mpl.rc_context({'figure.figsize': [8, 8],
#                      'axes.facecolor': 'white'}):
#     sc.pl.umap(adata_vis_plt, color=['sample'], size=30,
#                color_map = 'RdPu', ncols = 1, #legend_loc='on data',
#                legend_fontsize=10)

# filter the object
# adata_ref = adata_ref[:, selected].copy()
# adata_ref_hvg2000=adata_ref[:,hvg_2000_genes].copy()
# intersect_ref_hvg2000 = np.intersect1d(inf_aver.index, hvg_2000.index)
# intersect_ref_hvg2000_vis=np.intersect1d(intersect_ref_hvg2000,adata_vis.var_names)
# # adata_ref_hvg2000 = adata_ref[:, intersect_ref_hvg2000].copy()
# inf_aver_hvg2000=inf_aver.loc[intersect_ref_hvg2000,:].copy()
# slide_hvg2000=slide[:,intersect_ref_hvg2000].copy()
# adata_ref_hvg2000=adata_ref[:,intersect_ref_hvg2000].copy()

adata_ref_hvg2000=adata_ref[:,intersect_ref_hvg2000].copy()
adata_vis_hvg2000=adata_vis[:,inter]
# prepare anndata for the regression model
scvi.data.setup_anndata(adata=adata_ref_hvg2000,
                        # 10X reaction / sample / batch
                        batch_key='orig.ident',
                        # cell type, covariate used for constructing signatures
                        labels_key='newCelltype',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        # categorical_covariate_keys=[''] #Sample
                       )
scvi.data.view_anndata_setup(adata_ref_hvg2000)

# create and train the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref_hvg2000)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=False)

# plot ELBO loss history during training, removing first 20 epochs from the plot
mod.plot_history(20)
plt.savefig('ref.train.png')

#########
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref_hvg2000 = mod.export_posterior(
    adata_ref_hvg2000, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


mod.plot_QC()
plt.savefig('qc.png')


mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref_hvg2000.varm.keys():
    inf_aver = adata_ref_hvg2000.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref_hvg2000.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref_hvg2000.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref_hvg2000.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref_hvg2000.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# find shared genes and subset both anndata and reference signatures
slide = select_slide(adata_vis, sample_name[0])

intersect = np.intersect1d(slide.var_names, inf_aver.index)
slide = slide[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
scvi.data.setup_anndata(adata=slide, batch_key="sample")
scvi.data.view_anndata_setup(slide)


# create and train the model
mod = cell2location.models.Cell2location(
    slide, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=10,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
)

mod.train(max_epochs=10000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=False)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);
plt.savefig('full_data_training.png')
def select_slide(adata, s, s_col='sample'):
    r""" This function selects the data for one slide from the spatial anndata object.
    
    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param s_col: column in adata.obs listing experiment name for each location
    """
    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]
    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}
    return slide

# adata_vis = mod.export_posterior(#mod.adata.n_obs
#     adata_vis, sample_kwargs={'num_samples': 20, 'batch_size':20 , 'use_gpu': False}
# )
adata_vis_hvg2000 = mod_adata_vis.export_posterior(#mod.adata.n_obs
    adata_vis_hvg2000, sample_kwargs={'num_samples': 1000, 'batch_size':mod_adata_vis.adata.n_obs , 'use_gpu': False}
)
# Save model
mod_adata_vis.save("mod_4samples", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/4samples.h5ad"
adata_vis_hvg2000.write(adata_file)
adata_file
# Examine reconstruction accuracy to assess if there are any issues with mapping
# the plot should be roughly diagonal, strong deviations will signal problems
mod.plot_QC()
plt.savefig('vis_plot_qc.png')
fig = mod.plot_spatial_QC_across_batches()
plt.savefig('vis_fig_qc_batches.png')

#####
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
slide.obs[slide.uns['mod']['factor_names']] = slide.obsm['q05_cell_abundance_w_sf']

# select one slide
from cell2location.utils import select_slide
# slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')
# slide = select_slide(adata_vis, sample_name[0])

# plot in spatial coordinates
type_name=["Ductal_2","Endothelial","T_cell","Endocrine" ,"Plasma_B_cell", "B_cel","Fibroblasts_2","Acinar","Ductal_4" ,"Stellate","Fibroblasts_1", "Ductal_3","Macrophage","Ductal_1"]
# color_index=dict() 
color_index=dict(zip([i for i in range(14)],type_name))

for i,j  in color_index.items():
  slide.obs[j]=slide.obs[i]
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 30})
import matplotlib
SMALL_SIZE = 20
matplotlib.rc('font', size=SMALL_SIZE)
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [9, 10]}):
    sc.pl.spatial(slide, cmap='magma',
                  # show first 8 cell types
                  color=type_name,
                  # color=['Acinar cell', 'Alpha cell', 'Ductal cell', 'Endothelial cell', 'Fibroblast', 'Macrophage', 'Smooth muscle cell'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin='p5', vmax='p95'
                 )
plt.savefig(sample_name[0]+'_vis_predict_cell_abundanceV2.png')
plt.savefig(sample_name[0]+'_vis_predict_cell_abundanceV2.pdf')

# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
# clust_labels = ['Acinar cell', 'Alpha cell', 'Ductal cell', 'Endothelial cell', 'Fibroblast', 'Macrophage', 'Smooth muscle cell']
# clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
clust_labels=['Ductal_2','Endothelial','Fibroblasts_1']
with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_labels, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )
plt.savefig(sample_name[0]+'_3celltype_vis.png')
plt.savefig(sample_name[0]+'_3celltype_vis.pdf')

####下游分析
##### 重新加载vis数据
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
adata_file='data/cell2location_map/sp.h5ad'
adata_vis = sc.read_h5ad(adata_file)
#####重新加载单细胞数据
adata_ref = sc.read_h5ad('data/reference_signatures/sc.h5ad')

# 通过 Leiden 聚类识别离散的组织区域
# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=1.1)

# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")
# compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata_vis, min_dist = 0.3, spread = 1)

# show regions in UMAP coordinates
with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_vis, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20)
    # sc.pl.umap(adata_vis, color=['sample'], size=30,
    #            color_map = 'RdPu', ncols = 2,
    #            legend_fontsize=20)
plt.savefig('cluster_umap_vis.png')
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [8, 5]}):
    sc.pl.spatial(adata_vis, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.6)

plt.savefig('cluster_spatial_vis.png')

######################空间共定位分析
from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(10,12), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{run_name}/CoLocatedComb/'}
)

# Here we plot the NMF weights (Same as saved to `cell_type_fractions_heatmap`)
res_dict['n_fact11']['mod'].plot_cell_type_loadings()
plt.savefig('n_fact11.png')

#########计算spot中基因的表达量
adata_file = f"{run_name}/mod.h5ad"
mod.adata.write(adata_file)
# mod.adata=sc.read_h
# Compute expected expression per cell type
expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata
)

# Add to anndata layers
for i, n in enumerate(mod.factor_names_):
    adata_vis.layers[n] = expected_dict['mu'][i]

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file
