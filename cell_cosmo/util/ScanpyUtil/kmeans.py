#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : kmeans.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import warnings
import numpy as np
import pandas as pd
from cell_cosmo.util import runtime, PathUtil
from natsort import natsorted
from sklearn.cluster import KMeans
from .archive_data import archive_data
from cell_cosmo.util.ScanpyUtil.constant import (
    KMEANS_FEATURE_START, KMEANS_FEATURE_NAMES, ClusterSubDir, ClusterFile, assemble_cluster_key
)


def write(adata):
    for axis_type in {'umap', 'tsne'}:
        for feature in KMEANS_FEATURE_NAMES:
            feature_name = assemble_cluster_key(axis_type, feature)
            PathUtil.create_dir_if_not_exists(feature_name)
            df = adata.obs[[feature_name]].copy()
            df.columns = ["Cluster"]
            df['Cluster'] = df['Cluster'].map(lambda x: int(x) + 1)
            # df['Cluster'] = df['Cluster'].map(lambda x: int(x) + 1)
            df.to_csv(f"{feature_name}/{ClusterFile}", index_label="Barcode")


@runtime(__name__)
@archive_data(subdir=ClusterSubDir, write=write)
def kmeans(self):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r"KMeans is known to have a memory leak on Windows with MKL")
        # warnings.simplefilter("ignore")
        # 之前必须运行tsne,umap
        axis_df = self.adata.obsm.to_df()
        axis_columns = set(axis_df.columns.tolist())
        assert 'X_tsne1' in axis_columns and 'X_tsne2' in axis_columns, 'must run tsne first!'
        assert 'X_umap1' in axis_columns and 'X_umap2' in axis_columns, 'must run umap first!'
        df_tsne = axis_df[['X_tsne1', 'X_tsne2']]
        df_umap = axis_df[['X_umap1', 'X_umap2']]
        random_state = 0
        for axis_type, df in {'umap': df_umap, 'tsne': df_tsne}.items():
            for k, feature in enumerate(KMEANS_FEATURE_NAMES, start=KMEANS_FEATURE_START):
                feature_name = assemble_cluster_key(axis_type, feature)
                groups = KMeans(n_clusters=k, random_state=random_state).fit_predict(df)
                self.adata.obs[feature_name] = pd.Categorical(
                    values=groups.astype('U'),
                    categories=natsorted(map(str, np.unique(groups))),
                )
                self.adata.uns[feature_name] = {'param': {'k': k, 'random_state': random_state}}
