#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : scale.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import scanpy as sc
from cell_cosmo.util import runtime
from .archive_data import archive_data


@runtime(__name__)
@archive_data()
def scale(self):
    sc.pp.scale(
        self.adata,
        zero_center=True,
        max_value=10,
        copy=False,
        layer=None,
        obsm=None
    )
