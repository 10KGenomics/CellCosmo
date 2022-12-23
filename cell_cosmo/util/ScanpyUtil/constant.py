#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : util.py
@Time       : 2022/07/05
@Version    : 1.0
@Desc       : None
"""
# L2FC full Name
L2FC_FN = 'Log2 fold change'
# adj p value full name
ADJ_P_VAL_FN = 'Adjusted p value'
GRAPH_BASED_KEY = "louvain"
TopDirFmt = "%s_analysis"
ClusterSubDir = "clustering"
ClusterFile = "clusters.csv"
DiffExpSubDir = "diffExp"
DiffExpFile = "differential_expression.csv"
CLUSTER_SUFFIX = "_clusters"
GRAPH_BASE_DIR = f"graph{CLUSTER_SUFFIX}"

KMEANS_FEATURE_START = 2
# 2-10
KMEANS_FEATURE_NAMES = [f"kmeans_{k}{CLUSTER_SUFFIX}" for k in range(KMEANS_FEATURE_START, 11)]

NAME_DICT = {
    GRAPH_BASED_KEY: "Graph-based",  # use for df name trans
    GRAPH_BASE_DIR: "Graph-based",  # use for dir name trans
    **{n: f"K-means (K={k})" for k, n in enumerate(
        KMEANS_FEATURE_NAMES, start=KMEANS_FEATURE_START)},
}

FEATURE = "Feature"
FEATURE_COLS = ("ID", "Name")  # ID->gene_ids,Name->gene_name
FEATURE_ID, FEATURE_NAME = (f"{FEATURE} {col}" for col in FEATURE_COLS)


def assemble_cluster_key(dr, feature):
    """dr: umap or tsne, feature: kmeans for 2-9 groups"""
    return f"{dr}{feature}"
