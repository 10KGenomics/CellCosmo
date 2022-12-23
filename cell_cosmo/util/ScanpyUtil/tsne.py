#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : tsne.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import warnings
import scanpy as sc
import pandas as pd
from cell_cosmo.util import runtime, PathUtil
from .archive_data import archive_data


def write(adata):
    df = adata.obsm.to_df()
    assert isinstance(df, pd.DataFrame)
    cols = df.columns.tolist()
    tsne_cols = [c for c in cols if str(c).startswith('X_tsne')]
    PathUtil.create_dir_if_not_exists(f"{len(tsne_cols)}_components")
    df[tsne_cols].to_csv(f"{len(tsne_cols)}_components/projection.csv", index_label="Barcode")


@runtime(__name__)
@archive_data(write=write)
def tsne(self, n_pcs, thread):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r"In previous versions of scanpy,")
        sc.tl.tsne(
            self.adata,
            n_pcs=n_pcs,
            n_jobs=thread,
            copy=False,
        )
