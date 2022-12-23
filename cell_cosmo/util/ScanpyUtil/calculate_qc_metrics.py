#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : calculate_qc_metrics.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import scanpy as sc
from cell_cosmo.util import runtime, reader
from .archive_data import archive_data


@runtime(__name__)
@archive_data()
def calculate_qc_metrics(self, mt_gene_list, mito_var: str):
    if mt_gene_list:
        mito_genes = []
        for line, *_ in reader(mt_gene_list, ignore_test_env=True):
            if str(line).strip():
                mito_genes.append(str(line).strip())
        self.adata.var[mito_var] = self.adata.var_names.map(lambda x: True if x in mito_genes else False)
        # if not astype(bool), it will be type object and raise an error
        # https://github.com/theislab/anndata/issues/504
        self.adata.var[mito_var] = self.adata.var[mito_var].astype(bool)
    else:
        self.adata.var[mito_var] = self.adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        self.adata,
        qc_vars=[mito_var],
        percent_top=None,
        use_raw=False,
        log1p=False,
        inplace=True
    )
