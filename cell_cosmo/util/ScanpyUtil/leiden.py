#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : leiden.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import scanpy as sc
from cell_cosmo.util import runtime
from .archive_data import archive_data


@runtime(__name__)
@archive_data()
def leiden(self, resolution):
    sc.tl.leiden(
        self.adata,
        resolution=resolution,
        restrict_to=None,
        random_state=0,
        key_added='leiden',
        adjacency=None,
        directed=True,
        use_weights=True,
        n_iterations=-1,
        partition_type=None,
        neighbors_key=None,
        obsp=None,
        copy=False
    )
