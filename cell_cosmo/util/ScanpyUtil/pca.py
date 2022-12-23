#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : pca.py
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
    pca_cols = [c for c in cols if str(c).startswith('X_pca')]
    PathUtil.create_dir_if_not_exists(f"{len(pca_cols)}_components")
    df[pca_cols].to_csv(f"{len(pca_cols)}_components/projection.csv", index_label="Barcode")
    vr_df = pd.DataFrame(adata.uns['pca']['variance_ratio'], columns=['variance_ratio'])
    vr_df.index = range(1, len(vr_df) + 1)
    vr_df.to_csv(f"{len(pca_cols)}_components/variance.csv", index_label="PC")


@runtime(__name__)
@archive_data(write=write)
def pca(self):
    sc.pp.pca(
        self.adata,
        n_comps=50,
        zero_center=True,
        svd_solver='auto',
        random_state=0,
        return_info=False,
        use_highly_variable=True,
        dtype='float32',
        copy=False,
        chunked=False,
        chunk_size=None
    )
