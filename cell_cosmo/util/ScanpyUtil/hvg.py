#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : hvg.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import scanpy as sc
import numpy as np
from cell_cosmo.util import runtime
from .archive_data import archive_data


@runtime(__name__)
@archive_data()
def hvg(self):
    sc.pp.highly_variable_genes(
        self.adata,
        layer=None,
        n_top_genes=None,
        min_disp=0.5,
        max_disp=np.inf,
        min_mean=0.0125,
        max_mean=3,
        span=0.3,
        n_bins=20,
        flavor='seurat',
        subset=False,
        inplace=True,
        batch_key=None,
        check_values=True
    )
