#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : neighbors.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import scanpy as sc
from cell_cosmo.util import runtime
from .archive_data import archive_data


@runtime(__name__)
@archive_data()
def neighbors(self, n_pcs):
    sc.pp.neighbors(
        self.adata,
        n_neighbors=15,
        n_pcs=n_pcs,
        use_rep=None,
        knn=True,
        random_state=0,
        method='umap',
        metric='euclidean',
        key_added=None,
        copy=False
    )
