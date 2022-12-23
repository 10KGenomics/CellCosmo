#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : base_scanpy_runner.py
@Time       : 2022/06/14
@Version    : 1.0
@Desc       : None
"""
import os
import warnings
import scanpy as sc
import pandas as pd
from abc import ABCMeta
from cell_cosmo.output_runner.base_report_runner import BaseReportRunner
from cell_cosmo.util import GenomeUtil, runtime, reader, ScanpyUtil, PathUtil

import glob

# markers adjust p_value
PVAL_CUTOFF = 0.05
# scanpy mitochonrial variable name
MITO_VAR = 'mito'
NORMALIZED_LAYER = 'normalised'
RESOLUTION = 1.2
N_PCS = 25
MITO_GENE_PERCENT_LIST = [5, 10, 15, 20, 50]
# output marker top n in html
MARKER_TOP_N = 50
# marker sort by
MARKER_SORT_BY = 'p_val_adj'


class ScanpyWrapper(BaseReportRunner, metaclass=ABCMeta):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # data
        matrix_file = kwargs.pop("matrix_file")
        genomeDir = kwargs.pop("genomeDir")
        self.adata = sc.read_10x_mtx(
            matrix_file,
            var_names='gene_symbols',
        )
        self.mt_gene_list = GenomeUtil.parse_rna_dir(genomeDir)['mt_gene_list']

        self.kmeans_feature_names = [f"K{k}" for k in range(2, 11)]
        self.cluster_feature_names = ['cluster', 'louvain'] + self.kmeans_feature_names

        # out
        self.df_marker_file = f'{self.out_prefix}_markers.tsv'
        self.df_marker_file_glob = f'{self.out_prefix}_markers'
        self.df_marker_raw_file = f'{self.out_prefix}_markers_raw.tsv'
        self.df_tsne_file = f'{self.out_prefix}_tsne_coord.tsv'
        self.h5ad_file = f'{self.out_prefix}.h5ad'

    @runtime(f"{__name__}.write_mito_stats")
    def write_mito_stats(self):
        mt_pct_var = f'pct_counts_{MITO_VAR}'
        total_cell_number = self.adata.n_obs

        for mito_gene_percent in MITO_GENE_PERCENT_LIST:
            cell_number = sum(self.adata.obs[mt_pct_var] > mito_gene_percent)
            fraction = round(cell_number / total_cell_number * 100, 2)
            self.add_metric(
                name=f'Fraction of cells have mito gene percent>{mito_gene_percent}%',
                value=f'{fraction}%',
            )

    def run(self):
        ScanpyUtil.calculate_qc_metrics(self, self.mt_gene_list, MITO_VAR)
        self.write_mito_stats()
        ScanpyUtil.normalize(self, NORMALIZED_LAYER)  # run sc.pp.log1p at same time
        ScanpyUtil.hvg(self)
        ScanpyUtil.scale(self)
        ScanpyUtil.pca(self)
        ScanpyUtil.neighbors(self, N_PCS)
        ScanpyUtil.tsne(self, N_PCS, self.thread)
        ScanpyUtil.umap(self)
        try:
            ScanpyUtil.leiden(self, RESOLUTION)
            ScanpyUtil.louvain(self)
            ScanpyUtil.kmeans(self)
            ScanpyUtil.find_marker_genes(self, NORMALIZED_LAYER)
        except Exception as e:
            print("============????????")
            print(str(e))
        ScanpyUtil.output(self.adata, self.sample, self.outdir)

    # def get_df(self):
    #     """
    #     return df_tsne, df_marker
    #     """
    #
    #     df_tsne = read_tsne(self.df_tsne_file)
    #     df_marker = pd.read_csv(self.df_marker_file, sep="\t")
    #     df_marker = format_df_marker(df_marker)
    #     return df_tsne, df_marker
    def get_plot_data(self):
        df_dict = ScanpyUtil.read_cluster(self.sample, self.outdir)
        # df_marker = pd.read_csv(self.df_marker_file, sep="\t")
        # df_marker = format_df_marker(df_marker)
        marker_dict = ScanpyUtil.read_marker_genes(self.sample, self.outdir)
        return df_dict, marker_dict
