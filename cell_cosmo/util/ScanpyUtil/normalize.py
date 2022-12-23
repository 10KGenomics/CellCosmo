#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : normalize.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import scanpy as sc
from cell_cosmo.util import runtime
from .archive_data import archive_data


@runtime(__name__)
@archive_data()
def normalize(self, layer: str):
    self.adata.raw = self.adata
    sc.pp.normalize_total(
        self.adata,
        target_sum=1e4,
        inplace=True,
    )
    sc.pp.log1p(
        self.adata,
    )
    self.adata.layers[layer] = self.adata.X
