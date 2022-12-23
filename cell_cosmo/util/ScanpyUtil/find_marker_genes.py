#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : find_marker_genes.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import sys
import warnings
import pandas as pd
import scanpy as sc
from cell_cosmo.util import runtime, PathUtil
from .archive_data import archive_data
from .constant import (DiffExpSubDir, DiffExpFile, GRAPH_BASE_DIR,
                       KMEANS_FEATURE_NAMES, FEATURE, FEATURE_COLS,
                       GRAPH_BASED_KEY, L2FC_FN, ADJ_P_VAL_FN, assemble_cluster_key)


def format_df_marker_bak(df_marker):
    PVAL_CUTOFF = 0.05
    avg_logFC_col = "avg_log2FC"  # seurat 4
    if "avg_logFC" in df_marker.columns:  # seurat 2.3.4
        avg_logFC_col = "avg_logFC"
    df_marker = df_marker.loc[:, ["cluster", "gene", avg_logFC_col, "pct.1", "pct.2", "p_val_adj"]]
    df_marker["cluster"] = df_marker["cluster"].apply(lambda x: f"cluster {x}")
    df_marker = df_marker[df_marker["p_val_adj"] < PVAL_CUTOFF]

    return df_marker


def sort_key(t):
    n, d = t
    # 移除df中的cluster列,并将gene设置为index
    d.drop(["group"], axis=1, inplace=True)
    # d = d.round({'L2FC': 3, 'p-value': 3})
    d.set_index(["names"], inplace=True)
    # d['L2FC'] = d['L2FC'].apply(lambda x: format(x, '.2%s'))
    return n


def cmp(col_name: str):
    cols = col_name.split(" ")
    if cols[0] == FEATURE:
        # ID 本來就在 Name 前面
        return 0, cols[1]
    elif cols[0] == "Cluster":
        MAP_DICT3 = {"Mean": 0, "Log2": 1, "Adjusted": 2}
        return 1, int(cols[1]), MAP_DICT3.get(cols[2], sys.maxsize)


def trans_group_by_cluster(df: pd.DataFrame, df0: pd.DataFrame, cdf: pd.DataFrame):
    # 将分组下标从1开始
    df['group'] = df['group'].map(lambda x: f"Cluster {int(x) + 1}")
    # group,names,scores,logfoldchanges,pvals,pvals_adj,pct_nz_group,pct_nz_reference
    df = df[['group', 'names', 'logfoldchanges', 'pvals_adj']]
    df.columns = ['group', 'names', L2FC_FN, ADJ_P_VAL_FN]
    dds = sorted(list(df.groupby("group")), key=sort_key)
    new_df = pd.concat({FEATURE: df0, **{n: sdd for n, sdd in dds}}, axis=1)
    # Cluster 1 Mean Counts,Cluster 1 Log2 fold change,Cluster 1 Adjusted p value,
    new_df.columns = [f"{c[0]} {c[1]}" for c in new_df.columns.tolist()]
    mdf = new_df.join(cdf)
    # 调整标题列的顺序
    order = sorted(mdf.columns.tolist(), key=lambda x: cmp(x))
    mdf = mdf[order]
    return mdf


def get_mean_count(adata, cluster):
    dfs = []
    for g, df in adata.obs[[cluster]].groupby(cluster):
        mean_df = pd.DataFrame(
            data=adata.raw[df.index, adata.var_names].X.A,
            index=df.index,
            columns=adata.var_names,
        ).mean(axis=0)
        # 将分组下标从1开始
        dfs.append(pd.DataFrame(mean_df, columns=[f"Cluster {int(g) + 1} Mean Counts"]))
    ndf = pd.concat(dfs, axis=1)
    return ndf


def write(adata):
    # rank_genes_groups
    PathUtil.create_dir_if_not_exists(GRAPH_BASE_DIR)
    df0 = pd.DataFrame(
        data={FEATURE_COLS[0]: adata.var["gene_ids"],  # ENSG
              FEATURE_COLS[1]: adata.var.index},  # gene name
        index=adata.var.index
    )
    df_markers = sc.get.rank_genes_groups_df(adata, group=None, key=f"{GRAPH_BASED_KEY}_marker")
    count_df = get_mean_count(adata, GRAPH_BASED_KEY)
    trans_df = trans_group_by_cluster(df_markers, df0, count_df)
    trans_df.to_csv(f"{GRAPH_BASE_DIR}/{DiffExpFile}", index=None)
    for dr in ("umap", "tsne"):
        for feature in KMEANS_FEATURE_NAMES:
            feature_name = assemble_cluster_key(dr, feature)
            PathUtil.create_dir_if_not_exists(feature_name)
            df_markers = sc.get.rank_genes_groups_df(adata, group=None, key=f"{feature_name}_marker")
            count_df = get_mean_count(adata, feature_name)
            trans_df = trans_group_by_cluster(df_markers, df0, count_df)
            trans_df.to_csv(f"{feature_name}/{DiffExpFile}", index=None)


def run_rank_genes_groups(adata, group, key, layer):
    sc.tl.rank_genes_groups(
        adata,
        groupby=group,
        key_added=key,
        reference='rest',
        pts=True,
        method="wilcoxon",
        use_raw=False,
        layer=layer
    )


@runtime(__name__)
@archive_data(subdir=DiffExpSubDir, write=write)
def find_marker_genes(self, layer):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r"DataFrame is highly fragmented.")
        # cluster 对应 louvain
        louvain, leiden = "louvain", "leiden"
        run_rank_genes_groups(self.adata, louvain, f"{louvain}_marker", layer)
        run_rank_genes_groups(self.adata, leiden, f"{leiden}_marker", layer)

        for dr in ("umap", "tsne"):
            for feature in KMEANS_FEATURE_NAMES:
                feature_name = assemble_cluster_key(dr, feature)
                run_rank_genes_groups(self.adata, feature_name, f"{feature_name}_marker", layer)
