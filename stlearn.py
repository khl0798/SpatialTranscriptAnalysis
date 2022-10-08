import stlearn as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# data_dir='/data7/konghl/data/SpatialTranscriptome/E20210113/RNA/X1012-houai/outs'
data_dir='/data7/konghl/data/SpatialTranscriptome/E20210156/RNA/X4-19S58966-002/outs'
# st.add.image(adata=data,imgpath=data_dir+'/spatial/tissue_hires_image.png', library_id="X4-19S58966-002", visium=True)
# spot_mixtures = pd.read_csv(data_dir+'../../Brad/label_transfer_bc.csv', index_col=0, sep=',')

# Loading raw data #
data = st.Read10X(data_dir)
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath=data_dir+"/spatial/tissue_hires_image.png",
             library_id="X4-19S58966-002", visium=True)

# Basic normalisation #
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data) # NOTE: no log1p

# Adding the label transfer results,  #
###########找出每一行最大的索引，
#   df['max_idx'] = df.idxmax(axis=1) #求一行的最大值对应的索引
#   df['max_val']= df.max(axis=1) #取出该最大值
spot_mixtures = pd.read_csv('/data7/konghl/data/SpatialTranscriptome/E20210156/analysis_X3_X4_lognorm/RCTD_20220412/cell_to_cluster_weights.xls', 
                index_col=0, sep='\t')
new_values=list(spot_mixtures.index)
new_values_id=[i.replace('-1_','-') for i in new_values]
spot_mixtures.index=new_values_id
data.obs_names=[i.replace('-1','-2') for i in data.obs_names]
# spot_mixtures=
# labels = spot_mixtures.loc[:,'predicted.id'].values.astype(str)
# spot_mixtures = spot_m.oixtures.drop(['predicted.id','prediction.score.max'],
                                   # axis=1)
# spot_mixtures.columns = [col.replace('prediction.score.', '')
#                          for col in spot_mixtures.columns]

# Note the format! #
print(labels)
print(spot_mixtures)
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values)) # Check is in correct order

# NOTE: using the same key in data.obs & data.uns
# data.obs['cell_type'] = labels # Adding the dominant cell type labels per spot
# data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures.loc[data.obs_names] # Adding the cell type scores
cluster=pd.read_csv('/data7/konghl/data/SpatialTranscriptome/E20210156/analysis_X3_X4_lognorm/20220426/x4.select5.csv',
sep=",",index_col=0)
data1=data.copy()
data=data[cluster.index,]
data.obs['cluster']=cluster.loc[data.obs_names]
data.obs['cell_type']=labels.loc[data.obs_names]
data.obs['cell_type'] = data.obs['cell_type'].astype('category')

data1.uns['cell_type']=spot_mixtures.loc[data1.obs_names]
data1.obs['cell_type']=data1.uns['cell_type'].idxmax(axis=1)
##### 对obs提取数据，每个spot的注释得细胞类型，目前没有这个细胞类型，一个可能的解决办法是选择scores比例最高的细胞类型
###### 作为该spit的细胞类型
#st.pl.cluster_plot(data, use_label='cell_type')
#plt.savefig('cell_type.png')

#######
###NATMI 使用connectomeDB20220—（2293 manually curated ligand-receptor pairs)
# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(data, lrs, use_label='cell_type',
                  min_spots = 15, #Filter out any LR pairs with no scores for less than min_spots
                  distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=16, # Number of CPUs for parallel. If None, detects & use all available.
                  )
data_copy=data.copy()
data_copy.write('cluster5.h5ad')
#    adata: AnnData          Withthe counts of specified clusters in nearby spots stored as adata.uns['het']
######即统计邻域cluster的细胞类型score值>0.2的数目存储在obsm['cci_het']
st.tools.microenv.cci.het.count(data,use_label='cell_type')
lr_info = data.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
print('\n', lr_info)
#########展示前几个配体受体对得相互作用强度
#### data中的最强得配受体在整个得图种展示
best_lr = data.uns['lr_summary'].index.values[0] # Just choosing one of the top from lr_summary
stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
for i, stat in enumerate(stats):
    fig, axes = plt.subplots(ncols=1, figsize=(6,6))
    st.pl.lr_result_plot(data, use_result=stat, use_lr=best_lr, show_color_bar=True,image_alpha=0.6, ax=axes)
    axes.set_title(f'{best_lr} {stat}')
    plt.savefig('cluster5_top_lr_'+stat+'.png')


st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')

st.pl.lr_diagnostics(data, figsize=(10,2.5))
plt.savefig('diagnostics.png')                   
################# 配受体对得强度，展示形式1
st.pl.lr_plot(data, best_lr, inner_size_prop=0.9,outer_mode='binary', pt_scale=5,
              use_label=None, show_image=True,image_alpha=0.6,figsize=(6.4,8),fig=4,
              sig_spots=False)
plt.tight_layout()
plt.savefig('inner_size_prop.png')
### 展示形式2
st.pl.lr_plot(data, best_lr, outer_size_prop=0.9, outer_mode='binary', pt_scale=20,
              use_label=None, show_image=True,image_alpha=0.6,figsize=(6.4,8),fig=4,
              sig_spots=True)
plt.tight_layout()
plt.savefig('outer_size_prop.png')
##### 展示形式3  报错
# Only significant spots #
st.pl.lr_plot(data, best_lr,
              inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4,
              outer_mode='continuous', pt_scale=60,image_alpha=0.2,
              use_label=None, show_image=True,
              sig_spots=True)
plt.savefig('spot_continuous.png')
######展示形式4
st.pl.lr_plot(data, best_lr,
              inner_size_prop=0.08, middle_size_prop=.3, outer_size_prop=.5,
              outer_mode='binary', pt_scale=40,image_alpha=0.2,
              show_image=True, arrow_width=10, arrow_head_width=10,
              sig_spots=True, show_arrows=True,figsize=(6.4,8))
plt.tight_layout()
plt.savefig('spot_arrow.png')
#######展示形式5
st.pl.lr_plot(data, best_lr,
              inner_size_prop=0.08, middle_size_prop=.3, outer_size_prop=.5,
              outer_mode='binary', pt_scale=150,figsize=(6.4,8),image_alpha=0.7,
              use_label='cell_type', show_image=True,
              sig_spots=True)
plt.tight_layout()
plt.savefig('spot_siglr_celltype.png')

###### 展示形式6
# Showing the rankings of the LR from a global and local perspective.
# Ranking based on number of significant hotspots.
st.pl.lr_summary(data, n_top=50, figsize=(10,7))
plt.savefig('sig_lr_50.png')

# Running the counting of co-occurence of cell types and LR expression hotspots #
st.tl.cci.run_cci(data, 'cell_type', # Spot cell information either in data.obs or data.uns
                  min_spots=3, # Minimum number of spots for LR to be tested.
                  spot_mixtures=True, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0.2, # Spot considered to have cell type if score>0.2
                  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=1000 # Permutations of cell information to get background, recommend ~1000
                 )
###############选择存储
st.tl.cci.run_cci(data, 'cell_type', # Spot cell information either in data.obs or data.uns
                  min_spots=3, # Minimum number of spots for LR to be tested.
                  spot_mixtures=True, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0.2, # Spot considered to have cell type if score>0.2
                  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=100 # Permutations of cell information to get background, recommend ~1000
                 )
#########细胞类型展示形式1
st.pl.cci_check(data, 'cell_type')
plt.savefig('cell_type_diag.png')
# fig, axes = plt.subplots(ncols=2, figsize=(20,6))
# st.pl.lr_result_plot(data, use_result='lr_scores', use_lr=best_lr, image_alpha=0.4,show_color_bar=True, ax=axes[0])
# axes[0].set_title(f'{best_lr} lr_scores')
# st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=best_lr, image_alpha=0.4,show_color_bar=True, ax=axes[1])
# axes[1].set_title(f'{best_lr} lr_sig_scores')
# plt.savefig('top_'+best_lr+"_scores.png")

#####输出lr得表达矩阵文件
labels=pd.read_csv('/data7/konghl/data/SpatialTranscriptome/E20210156/analysis_X3_X4_lognorm/20220426/cluster5.labels.x4.csv',sep=",",index_col=0)
####################
# Visualising the no. of interactions between cell types across all LR pairs #
pos_1 = st.pl.ccinet_plot(data, 'cell_type', return_pos=True)
plt.savefig('sig_celltype_lr.png')

# Just examining the cell type interactions between selected pairs #
lrs = data.uns['lr_summary'].index.values[0:3]
for best_lr in lrs[0:3]:
    st.pl.ccinet_plot(data, 'cell_type', best_lr, min_counts=2,
                         figsize=(10,7.5), pos=pos_1,
                      )
    plt.savefig('sig_celltype_lr_'+best_lr+".png")
######################显著性配受体对玄图
st.pl.lr_chord_plot(data, 'cell_type')
plt.savefig('circos_sig_lr.png')
for lr in lrs:
    st.pl.lr_chord_plot(data, 'cell_type', lr)
    plt.savefig('circos_'+lr+"_sig.png")
########
# This will automatically select the top interacting CCIs and their respective LRs #
st.pl.lr_cci_map(data, 'cell_type', lrs=None, min_total=100, figsize=(15,10))
plt.tight_layout()
plt.savefig('dotplot_sig.png',dpi=600)

##########热图
from matplotlib.pyplot import savefig as fig
st.pl.cci_map(data, 'cell_type',figsize=(6,6))
plt.tight_layout()
plt.savefig('sig_celltype_heatmap.png',dpi=600)

lrs = data.uns['lr_summary'].index.values[0:3]
for lr in lrs[0:3]:
    st.pl.cci_map(data, 'cell_type', lr,figsize=(8,8))
    plt.tight_layout()
    plt.savefig('sig_celltype_heatmap_'+lr+"_png",dpi=600)
##########细胞类型相互作用图
best_lr = lrs[0]

### This will plot with simple black arrows ####
st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode=None,
              pt_scale=40, use_label='cell_type', show_arrows=True,
              show_image=True, sig_spots=False, sig_cci=True,
                 arrow_head_width=40,lr_colors={"ALB_LRP2":'Reds'},
                 arrow_width=10, cell_alpha=.8,image_alpha=0.3
                 )
plt.tight_layout()
plt.savefig('celltype_arrors.png',dpi=600)

### This will colour the spot by the mean LR expression in the spots connected by arrow
st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode=None,
              pt_scale=10, use_label='cell_type', show_arrows=True,
              show_image=True, sig_spots=False, sig_cci=True,
                 arrow_head_width=4, arrow_width=2,image_alpha=0.3,
                 arrow_cmap='Reds', arrow_vmax=1.5,
                 )
plt.tight_layout()
plt.savefig('celltype_arrors.png',dpi=600)
##########################输出表格信息
data.uns['lr_summary'].to_csv('lr_summary.csv')
data.uns['lr_cci_cell_type'].to_csv('lr_cci_cell_type.csv')
outs_info=["lr_scores", "p_vals", "p_adjs", "-log10(p_adjs)", "lr_sig_scores"]
for i in outs_info:
  lr_scores=pd.DataFrame(data.obsm[i])
  lr_scores.index=data.obs_names
  lr_scores.columns=data.uns['lr_summary'].index
  lr_scores.to_csv(i+'.csv')
  

