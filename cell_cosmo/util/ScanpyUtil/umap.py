#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : umap.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import scanpy as sc
import pandas as pd
from cell_cosmo.util import runtime, PathUtil
from .archive_data import archive_data


def write(adata):
    df = adata.obsm.to_df()
    assert isinstance(df, pd.DataFrame)
    cols = df.columns.tolist()
    umap_cols = [c for c in cols if str(c).startswith('X_umap')]
    PathUtil.create_dir_if_not_exists(f"{len(umap_cols)}_components")
    df[umap_cols].to_csv(f"{len(umap_cols)}_components/projection.csv", index_label="Barcode")


@runtime(__name__)
@archive_data(write=write)
def umap(self):
    sc.tl.umap(
        self.adata,
        min_dist=0.5,
        spread=1.0,
        n_components=2,
        maxiter=None,
        alpha=1.0,
        gamma=1.0,
        negative_sample_rate=5,
        init_pos='spectral',
        random_state=0,
        a=None,
        b=None,
        copy=False,
        method='umap',
        neighbors_key=None
    )
