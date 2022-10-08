import scanpy as sc
import celltypist
from celltypist import models
#adata_2000 = sc.read('celltypist_demo_folder/demo_2000_cells.h5ad',
# backup_url = 'https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/demo_2000_cells.h5ad')
adata_2000=sc.read('demo_2000_cells.h5ad')
# adata_2000.shape
adata_2000.X.expm1().sum(axis = 1)#### 每个细胞的最大值为10000
models.download_models(force_update = True)
# models.models_path
models.models_description()
# Indeed, the `model` argument defaults to `Immune_All_Low.pkl`.
model = models.Model.load(model = 'Immune_All_Low.pkl')
model.cell_types
model.description
# Not run; predict cell identities using this loaded model.
#predictions = celltypist.annotate(adata_2000, model = model, majority_voting = True)
# Alternatively, just specify the model name (recommended as this ensures the model is intact every time it is loaded).
predictions = celltypist.annotate(adata_2000, model = 'Immune_All_Low.pkl', majority_voting = True)
import matplotlib.pyplot as plt
sc.tl.umap(adata)
plt.savefig('test_umap.png')
sc.pl.umap(adata, color = ['cell_type', 'predicted_labels', 'majority_voting'], legend_loc = 'on data')
plt.savefig('test_umap.png')
celltypist.dotplot(predictions, use_as_reference = 'cell_type', use_as_prediction = 'majority_voting')
plt.savefig('test_dotplot.png')
#https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/model_from_immune2000.pkl
#######寻找各细胞类型的驱动基因
model = models.Model.load(model = 'celltypist_demo_folder/model_from_immune2000.pkl')
#######自定义model作为模型迁移数据
adata_400 = sc.read('celltypist_demo_folder/demo_400_cells.h5ad', backup_url = 'https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/demo_400_cells.h5ad')
# The `cell_type` in `adata_2000.obs` will be used as cell type labels for training.
new_model = celltypist.train(adata_2000, labels = 'cell_type', n_jobs = 10, feature_selection = True)
# Save the model.
new_model.write('celltypist_demo_folder/model_from_immune2000.pkl')
new_model = models.Model.load('celltypist_demo_folder/model_from_immune2000.pkl')
# Not run; predict the identity of each input cell with the new model.
#predictions = celltypist.annotate(adata_400, model = new_model, majority_voting = True)
# Alternatively, just specify the model path (recommended as this ensures the model is intact every time it is loaded).
predictions = celltypist.annotate(adata_400,
                                  model = 'celltypist_demo_folder/model_from_immune2000.pkl',
                                  majority_voting = True)
adata = predictions.to_adata()
sc.tl.umap(adata) ##运行umap
sc.pl.umap(adata, color = ['cell_type', 'predicted_labels', 'majority_voting'], legend_loc = 'on data')
plt.savefig('train_umap.png')
# Any model can be inspected.
# Here we load the previously saved model trained from 2,000 immune cells.
model = models.Model.load(model = 'celltypist_demo_folder/model_from_immune2000.pkl')
# model.cell_types
top_3_genes = model.extract_top_markers("Mast cells", 3)
# top_3_genes
# Check expression of the three genes in the training set.
sc.pl.violin(adata_2000, top_3_genes, groupby = 'cell_type', rotation = 90)
plt.savefig('mastcell_top3genes_violin.png')
###################multiple-labels
adata_500 = sc.read('celltypist_demo_folder/demo_500_cells.h5ad', backup_url = 'https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/demo_500_cells.h5ad')
from sklearn.linear_model import SGDClassifier
X = [[0., 0.], [1., 1.]]
y = [0, 1]
clf = SGDClassifier(loss="hinge", penalty="l2", max_iter=5)
clf.fit(X, y)
SGDClassifier(alpha=0.0001, average=False, class_weight=None,
early_stopping=False, epsilon=0.1, eta0=0.0, fit_intercept=True,
l1_ratio=0.15, learning_rate="optimal", loss="hinge", max_iter=5,
n_iter_no_change=5, n_jobs=None, penalty="l2",
power_t=0.5, random_state=None, shuffle=True, tol=None,
validation_fraction=0.1, verbose=0, warm_start=False)
# n_iter=None,

#######多分类标签
from sklearn import datasets
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
iris = datasets.load_iris()
X, y = iris.data, iris.target
OneVsRestClassifier(LinearSVC(random_state=0)).fit(X, y).predict(X)
from tensorflow import keras
from sklearn.svm import SVC
import numpy as np
import pickle
(x_train, y_train), (x_test, y_test) = keras.datasets.mnist.load_data()
#转变数据集形式
x_train_transed=[]
for index in range(len(x_train)):
    x_train_transed.append(x_train[index].reshape(-1))
x_test_transed=[]
for index in range(len(x_test)):
    x_test_transed.append(x_test[index].reshape(-1))
print("转变完了")