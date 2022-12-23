#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : output.py
@Time       : 2022/07/05
@Version    : 1.0
@Desc       : None
"""
import os.path

import pandas as pd
from cell_cosmo.util import runtime, PathUtil
from cell_cosmo.util.ScanpyUtil.constant import (
    GRAPH_BASED_KEY, KMEANS_FEATURE_NAMES, NAME_DICT, assemble_cluster_key
)


@runtime(f"{__name__}.write_cluster")
def _write_cluster(adata, sample):
    axis_df = adata.obsm.to_df()
    df_tsne = axis_df[['X_tsne1', 'X_tsne2']].copy()
    df_umap = axis_df[['X_umap1', 'X_umap2']].copy()

    for axis_type, df in {"tsne": df_tsne, "umap": df_umap}.items():
        # 将聚类信息记录到表格文件
        df[GRAPH_BASED_KEY] = adata.obs[GRAPH_BASED_KEY]
        for feature in KMEANS_FEATURE_NAMES:
            feature_name = assemble_cluster_key(axis_type, feature)
            trans_name = feature_name.replace(axis_type, "")
            df[trans_name] = adata.obs[feature_name]
            # 对聚类团下标都+1
            df[trans_name] = df[trans_name].map(lambda x: int(x) + 1)
        df['Gene_Counts'] = adata.obs.n_genes_by_counts
        df['UMI_Counts'] = adata.obs.total_counts
        df[['UMI_Counts']] = df[['UMI_Counts']].astype(int)
        if axis_type == "tsne":
            name_dict = {'X_tsne1': 'tSNE_1', 'X_tsne2': 'tSNE_2'}
        else:
            name_dict = {'X_umap1': 'UMAP_1', 'X_umap2': 'UMAP_2'}
        df = df.rename(name_dict, axis='columns')
        df.to_csv(f'{sample}_{axis_type}_coord.tsv', sep='\t')


def read_cluster(sample, outdir):
    result = {}
    with PathUtil.chdir(outdir):
        for axis_type in {"tsne", "umap"}:
            df = pd.read_csv(f'{sample}_{axis_type}_coord.tsv', sep='\t')
            # compatible with old version
            if 'Unnamed: 0' in df.columns:
                df.rename(columns={'Unnamed: 0': 'barcode'}, inplace=True)
                df = df.set_index('barcode')
            df.rename(NAME_DICT, axis='columns', inplace=True)
            df['size'] = 5
            df['barcode_index'] = list(range(1, len(df) + 1))
            result[axis_type] = df
    return result


@runtime(f"{__name__}.write_h5ad")
def write_h5ad(adata, h5ad_file):
    adata.write(h5ad_file)


@runtime(__name__)
def output(adata, sample, outdir):
    """将之前所有的步骤运行完毕后再运行这里,输出h5ad,markers,tsne"""
    with PathUtil.chdir(outdir):
        # write_markers(adata, sample)
        _write_cluster(adata, sample)
        # self.df_marker_file_glob = f'{self.out_prefix}_markers'
        # self.df_marker_raw_file = f'{self.out_prefix}_markers_raw.tsv'
    write_h5ad(adata, f"{sample}.h5ad")
