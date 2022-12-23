#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : louvain.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import scanpy as sc
from cell_cosmo.util import runtime, PathUtil
from .archive_data import archive_data
from cell_cosmo.util.ScanpyUtil.constant import (
    GRAPH_BASED_KEY, CLUSTER_SUFFIX, ClusterSubDir, ClusterFile
)


def write(adata):
    write_dir = f"graph{CLUSTER_SUFFIX}"
    PathUtil.create_dir_if_not_exists(write_dir)
    # use .copy() fix SettingWithCopyWarning
    df = adata.obs[[GRAPH_BASED_KEY]].copy()
    df.columns = ['Cluster']
    df['Cluster'] = df['Cluster'].map(lambda x: int(x) + 1)
    df.to_csv(f"{write_dir}/{ClusterFile}", index_label="Barcode")


@runtime(__name__)
@archive_data(subdir=ClusterSubDir, write=write)
def louvain(self):
    sc.tl.louvain(
        self.adata,
        random_state=0,
        key_added=GRAPH_BASED_KEY,
        copy=False
    )
