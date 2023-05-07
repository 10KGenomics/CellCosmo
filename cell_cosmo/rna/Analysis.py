#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : Analysis.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import logging
import os.path

import pandas as pd
from cell_cosmo.util import runtime
# from sklearn.cluster import KMeans
# from cell_cosmo.tools.plotlyplot import TsnePlot
from cell_cosmo.output_runner.base_scanpy_runner import ScanpyWrapper

logger = logging.getLogger(__name__)


class Analysis(ScanpyWrapper):
    _STEP_NAME = "analysis"
    _DISPLAY_TITLE = "Analysis"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def process_df_to_plot_data(self, name, df: pd.DataFrame):
        def _trans_cluster(_feature_name):
            sum_df = df.groupby([_feature_name]).agg("count").iloc[:, 0]
            percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
            res_d = {}
            res_l = []
            for cluster in sorted(df[_feature_name].unique()):
                lbl_name = f"{cluster}({percent_df[cluster]}%)"
                res_d[cluster] = lbl_name
                res_l.append(lbl_name)
            df[_feature_name] = df[_feature_name].map(res_d)
            return res_l

        tsne1, tsn2, umap1, umap2 = 'tSNE_1', 'tSNE_2', 'UMAP_1', 'UMAP_2'
        not_feature_set = {tsne1, tsn2, umap1, umap2, 'size', 'barcode_index'}
        if name == 'tsne':
            ax1, ax2 = tsne1, tsn2
        else:
            ax1, ax2 = umap1, umap2
        features = df.columns.tolist()
        features = [f for f in features if f not in not_feature_set]
        # tSNE_1, tSNE_2, cluster, louvain, K-10, Gene_Counts, UMI_Counts
        # add: size, barcode_index
        X, Y = df[ax1].tolist(), df[ax2].tolist()
        size, bi = df["size"].tolist(), df["barcode_index"].tolist()
        plot_cluster = {"X": X, "Y": Y, "size": size, "barcode_index": bi, "KEY": []}
        plot_gene = {"X": X, "Y": Y, "size": size, "barcode_index": bi}

        for feature in features:
            if feature in {"cluster", "louvain"}:
                # already rename
                continue
            elif feature == "Gene_Counts":
                pass
                # NOTE 按gene数展示聚类关系(n_gene_by_counts)
            elif feature == "UMI_Counts":
                # NOTE 按umi数展示聚类关系(total_counts)
                plot_gene[feature] = df[feature].tolist()
            else:
                # "Graph-based"(leiden:cluster or louvain) and kmeans
                res_list = _trans_cluster(feature)
                plot_cluster[feature] = {
                    "data": df[feature].tolist(),
                    "category_orders": res_list,
                }
                plot_cluster["KEY"].append(feature)

        return {
            "cluster": plot_cluster,
            "gene": plot_gene,
        }

    @runtime(__name__)
    def run(self):
        super(Analysis, self).run()

        self.set_metric_list(metric_list=self.get_metric_list())
        self.add_marker_help()
        df_dict, marker_dict = self.get_plot_data()
        for name in df_dict:
            plot_dict = self.process_df_to_plot_data(name, df_dict[name])
            plot_dict['table_dict'] = {
                'title': 'Marker Genes by Cluster',
                'table': marker_dict[name],
                'id': f"{name}_marker_genes"}
            self.add_data(**{name: plot_dict})

        # todo 流程的最后一步，添加清理大文件的逻辑
        clean_type = os.getenv("CELLCOSMO_LARGE_FILE_CLEAN_STRATEGY", 1)
        if str(clean_type) == '1':
            need_clean_files = [
                f"{self.out_prefix}_Aligned.sortedByCoord.out.bam.featureCounts.bam",
                f"{self.out_prefix}_Aligned.sortedByCoord.out.bam",
                f"{self.out_prefix}_Aligned.sortedByCoord.out.bam.bai",
                f"{self.out_prefix}_Aligned.out.bam",
                f"{self.out_prefix}_no_polyt1.fq.gz",
                f"{self.out_prefix}_no_polyt2.fq.gz",
                f"{self.out_prefix}_1.fq.gz",
                f"{self.out_prefix}_2.fq.gz",
                f"{self.outdir}/.temp.db",
            ]
            for ncf in need_clean_files:
                if os.path.exists(ncf):
                    os.remove(ncf)

    def add_marker_help(self):
        self.add_help_content(
            name='Projection of Cells Colored by UMI counts',
            content=('Shown here are the total Gene counts for each cell-barcode. '
                     'Cells with greater counts likely have higher RNA content than '
                     'cells with fewer counts.The axes correspond to the 2-dimensional '
                     'embedding produced by the t-SNE or UMAP algorithm. In this space, '
                     'pairs of cells that are close to each other have more similar gene '
                     'expression profiles than cells that are distant from each other')
        )
        self.add_help_content(
            name='Projection of Cells by Clustering',
            content=('By an automated clustering algorithm, These valid cells that '
                     'have similar expression profiles will be together. The axes correspond '
                     'to the 2-dimensional embedding produced by the t-SNE or UMAP algorithm. '
                     'In this space, pairs of cells that are close to each other have more '
                     'similar gene expression profiles than cells that are distant from each other. '
                     'The display is limited to a random subset of cells')
        )
        self.add_help_content(
            name='Avg_log2FC',
            content=('Log fold-chage of the average expression between the markers for every '
                     'cluster compared to all remaining cells. A value of 1.0 indicates 2-fold greater '
                     'expression in the cluster of interest')
        )
        self.add_help_content(
            name='P value',
            content='A measure of the statistical significance of the expression difference and is '
                    'based on a negative binomial test. The p-value reported here has been adjusted '
                    'for multiple testing via the Benjamini-Hochberg procedure'
        )

    def collect_matrix(self):
        # merge Demultiplexing&Trimming to Sequencing,step is barcode&cutadapt
        self.merge_step("barcode", "cutadapt", "sequence")
        self.set_summary_step()
