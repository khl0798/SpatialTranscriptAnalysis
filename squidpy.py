#!coding: utf-8
import scanpy as sc
import anndata as ad
import squidpy as sq

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

# load the pre-processed dataset
img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()
sq.pl.spatial_scatter(adata, color="cluster")
plt.savefig('cluster.png')

# calculate features for different scales (higher value means more context)
for scale in [1.0, 2.0]:
    feature_name = f"features_summary_scale{scale}"
    sq.im.calculate_image_features(
        adata,
        img.compute(),
        features="summary",
        key_added=feature_name,
        n_jobs=4,
        scale=scale,
    )