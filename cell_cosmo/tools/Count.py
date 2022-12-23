#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : Count.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import os
import sys
import shutil
import random
import logging
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
from scipy.sparse import coo_matrix
from cell_cosmo.util import (
    runtime, bam2count_table, reader,
    ImageUtil, GenomeUtil, GTFDictUtil, PathUtil
)
# from cell_cosmo.util.CountDictUtil import CountDictUtil
from cell_cosmo.util.BarcodeCorrectUtil import correct_barcode_and_umi
from cell_cosmo.tools.aa.cell_calling_3 import cell_calling_3
from cell_cosmo.tools.aa import get_plot_elements
from cell_cosmo.tools.plotlyplot import Line_plot, Violin_plot  # TODO
from cell_cosmo.output_runner.base_report_runner import BaseReportRunner

logger = logging.getLogger(__name__)
RAW_MATRIX_DIR_SUFFIX = ['raw_feature_bc_matrix', 'all_matrix']
FILTERED_MATRIX_DIR_SUFFIX = ['filtered_feature_bc_matrix', 'matrix_10X']
MATRIX_FILE_NAME = ['matrix.mtx']
FEATURE_FILE_NAME = ['genes.tsv']
BARCODE_FILE_NAME = ['barcodes.tsv']
TOOLS_DIR = os.path.dirname(__file__)
random.seed(0)
np.random.seed(0)

# downsample.csv
READ_FRACTION = 'read_fraction'
MEDIAN_GENE_NUMBER = 'median_gene_number'
READ_SATURATION = 'read_saturation'
UMI_SATURATION = 'umi_saturation'

# Plot axis title in HTML
X_TITLE = 'Read Fraction'
SATURATION_Y_TITLE = 'Sequencing Saturation(%)'
MEDIAN_GENE_Y_TITLE = 'Median Genes per Cell'


class Count(BaseReportRunner):
    _DISPLAY_TITLE = "Cells"
    _STEP_NAME = "count"

    # barcode_diff, barcode_correct_limit,
    # umi_diff, umi_correct_limit,
    def __init__(self,
                 bam, genomeDir, n_umi_filter,
                 barcode_correct_limit, umi_correct_limit,
                 force_cell_num, expected_cell_num,
                 cell_calling_method, **kwargs):
        super(Count, self).__init__(**kwargs)
        self.n_umi_filter = n_umi_filter
        self.barcode_correct_limit = barcode_correct_limit
        self.umi_correct_limit = umi_correct_limit
        self.force_cell_num = force_cell_num
        self.cell_calling_method = cell_calling_method
        self.expected_cell_num = int(expected_cell_num)
        self.bam = bam

        # set
        self.gtf_file = GenomeUtil.parse_rna_dir(genomeDir)['gtf']
        self.mt_gene_list = GenomeUtil.parse_rna_dir(genomeDir)['mt_gene_list']
        self.gtf_dict = GTFDictUtil.GTFDict(self.gtf_file)
        self.down_sample_dict = {}
        self.qc_plot = {}
        self.qc_plot_name_map = {
            'n_genes_by_counts': 'nFeature_RNA',
            'total_counts': 'nCount_RNA',
            'pct_counts_mito': 'percent.mt'
        }

        # output files
        self.umi_correct_log = f'{self.out_prefix}_umi_correct.log'
        self.barcode_correct_log = f'{self.out_prefix}_barcode_correct.log'
        self.marked_count_file = f'{self.out_prefix}_counts.tsv'
        self.raw_matrix_dir = f'{self.out_prefix}_{RAW_MATRIX_DIR_SUFFIX[0]}'
        self.cell_matrix_dir = f'{self.out_prefix}_{FILTERED_MATRIX_DIR_SUFFIX[0]}'
        self.down_sample_file = f'{self.out_prefix}_downsample.tsv'
        # matrix
        self.CB_describe = None
        self.CB_total_Genes = None
        self.CB_reads_count = None
        self.reads_mapped_to_transcriptome = None

    def get_df_line(self):
        df_line = pd.read_csv(self.down_sample_file, sep="\t", header=0)
        df_line.rename(columns={
            READ_FRACTION: X_TITLE,
            MEDIAN_GENE_NUMBER: MEDIAN_GENE_Y_TITLE,
            UMI_SATURATION: SATURATION_Y_TITLE,
        }, inplace=True)

        return df_line

    @runtime(__name__)
    def run(self):
        df = bam2count_table(self.bam, self.out_prefix)

        df = correct_barcode_and_umi(
            df, self.outdir, self.sample,
            n_umi_filter=self.n_umi_filter,
            filter_limit=self.barcode_correct_limit,
            percent=self.umi_correct_limit,
            thread=self.thread)
        df.rename(columns={
            "barcode": "Barcode",
            "umi": "UMI",
            "gene_id": "geneID"
        }, inplace=True)

        # df_sum
        df_sum = Count.get_df_sum(df)

        # export all matrix
        self.write_matrix_10X(df, self.raw_matrix_dir)

        # call cells
        cell_bc, _threshold = self.cell_calling(df_sum)

        # get cell stats
        self.CB_describe = self.get_cell_stats(df_sum, cell_bc)

        # export cell matrix
        df_cell = df.loc[df['Barcode'].isin(cell_bc), :]
        self.write_matrix_10X(df_cell, self.cell_matrix_dir)
        # must set after cell_matrix_dir data write in
        self.set_qc_plot()
        (self.CB_total_Genes, self.CB_reads_count, self.reads_mapped_to_transcriptome) = self.cell_summary(
            df, cell_bc)

        # downsampling
        cell_bc = set(cell_bc)
        self.downsample(df_cell)

    def set_qc_plot(self):
        adata = sc.read_10x_mtx(
            self.cell_matrix_dir,
            var_names='gene_symbols',
        )
        mito_var = 'mito'
        if self.mt_gene_list:
            mito_genes = []
            for line, *_ in reader(self.mt_gene_list, ignore_test_env=True):
                if str(line).strip():
                    mito_genes.append(str(line).strip())
            adata.var[mito_var] = adata.var_names.map(lambda x: True if x in mito_genes else False)
            # if not astype(bool), it will be type object and raise an error
            # https://github.com/theislab/anndata/issues/504
            adata.var[mito_var] = adata.var[mito_var].astype(bool)
        else:
            adata.var[mito_var] = adata.var_names.str.upper().str.startswith('MT-')

        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=['mito'],
            percent_top=None,
            use_raw=False,
            log1p=False,
            inplace=True
        )
        self.qc_plot = {
            "option": [],
            "plot_data": [],
        }

        for key, rename_key in self.qc_plot_name_map.items():
            df_plot = adata.obs[[key]].copy()
            df_plot.columns = [rename_key]
            data = Violin_plot(df_plot, rename_key).get_plotly_div()
            class_name = f"qc_plot_{rename_key.replace('.', '_')}"
            self.qc_plot['option'].append({
                'title': rename_key,
                'class_name': class_name
            })
            self.qc_plot['plot_data'].append({
                'data': data,
                'class_name': class_name
            })

    @runtime(f"{__name__}.cell_calling")
    def cell_calling(self, df_sum):
        cell_calling_method = self.cell_calling_method

        if (self.force_cell_num is not None) and (self.force_cell_num != 'None'):
            cell_bc, UMI_threshold = self.force_cell(df_sum)
        elif cell_calling_method == 'auto':
            cell_bc, UMI_threshold = self.auto_cell(df_sum)
        elif cell_calling_method == 'EmptyDrops_CR':
            cell_bc, UMI_threshold = self.emptydrop_cr_cell(df_sum)
        else:
            raise Exception
        return cell_bc, UMI_threshold

    @runtime(f"{__name__}.force_cell")
    def force_cell(self, df_sum):

        force_cell_num = int(self.force_cell_num)

        df_barcode_count = df_sum.groupby(
            ['UMI']).size().reset_index(
            name='barcode_counts')
        sorted_df = df_barcode_count.sort_values("UMI", ascending=False)
        sorted_df["barcode_cumsum"] = sorted_df["barcode_counts"].cumsum()
        num_points = sorted_df.shape[0]
        index_low, index_high = None, None
        for i in range(num_points):
            if sorted_df.iloc[i, :]["barcode_cumsum"] >= force_cell_num:
                index_low = i - 1
                index_high = i
                break
        df_sub = sorted_df.iloc[index_low: index_high + 1, :]
        distance = abs(df_sub["barcode_cumsum"] - force_cell_num)
        actual_index = np.argmin(distance)
        threshold = df_sub.iloc[actual_index, :]['UMI']
        cell_bc = Count.get_cell_bc(df_sum, threshold, col='UMI')
        return cell_bc, threshold

    @staticmethod
    def find_threshold(df_sum, idx):
        return int(df_sum.iloc[idx - 1, df_sum.columns == 'UMI'])

    @staticmethod
    def get_cell_bc(df_sum, threshold, col='UMI'):
        return list(df_sum[df_sum[col] >= threshold].index)

    @runtime(f"{__name__}.auto_cell")
    def auto_cell(self, df_sum):
        idx = int(self.expected_cell_num * 0.01)
        barcode_number = df_sum.shape[0]
        idx = int(min(barcode_number, idx))
        if idx == 0:
            sys.exit("cell number equals zero!")
        # calculate read counts threshold
        threshold = int(Count.find_threshold(df_sum, idx) * 0.1)
        threshold = max(1, threshold)
        cell_bc = Count.get_cell_bc(df_sum, threshold)

        return cell_bc, threshold

    @runtime(f"{__name__}.emptydrop_cr_cell")
    def emptydrop_cr_cell(self, df_sum):
        cell_bc, initial_cell_num = cell_calling_3(self.raw_matrix_dir, self.expected_cell_num)
        threshold = Count.find_threshold(df_sum, initial_cell_num)
        return cell_bc, threshold

    @staticmethod
    @runtime(f"{__name__}.get_df_sum")
    def get_df_sum(df, col='UMI'):
        def num_gt2(x):
            return pd.Series.sum(x[x > 1])

        df_sum = df.groupby('Barcode').agg({
            'count': ['sum', num_gt2],
            'UMI': 'count',
            'geneID': 'nunique'
        })
        df_sum.columns = ['readcount', 'UMI2', 'UMI', 'geneID']
        df_sum = df_sum.sort_values(col, ascending=False)
        return df_sum

    def get_cell_stats(self, df_sum, cell_bc):
        df_sum.loc[:, 'mark'] = 'UB'
        df_sum.loc[df_sum.index.isin(cell_bc), 'mark'] = 'CB'
        df_sum.to_csv(self.marked_count_file, sep='\t')
        CB_describe = df_sum.loc[df_sum['mark'] == 'CB', :].describe()
        return CB_describe

    @runtime(f"{__name__}.write_matrix_10X")
    def write_matrix_10X(self, df, matrix_dir):
        if not os.path.exists(matrix_dir):
            os.mkdir(matrix_dir)

        df_UMI = df.groupby(['geneID', 'Barcode']).agg({'UMI': 'count'})
        mtx = coo_matrix((df_UMI.UMI, (df_UMI.index.codes[0], df_UMI.index.codes[1])))
        gene_id = df_UMI.index.levels[0].to_series()
        # add gene symbol
        gene_name = gene_id.apply(lambda x: self.gtf_dict[x])
        genes = pd.concat([gene_id, gene_name], axis=1)
        genes.columns = ['gene_id', 'gene_name']

        barcodes = df_UMI.index.levels[1].to_series()
        genes.to_csv(f'{matrix_dir}/{FEATURE_FILE_NAME[0]}', index=False, sep='\t', header=False)
        barcodes.to_csv(f'{matrix_dir}/{BARCODE_FILE_NAME[0]}', index=False, sep='\t', header=False)
        mmwrite(f'{matrix_dir}/{MATRIX_FILE_NAME[0]}', mtx)

    @runtime(f"{__name__}.cell_summary")
    def cell_summary(self, df, cell_bc):

        df.loc[:, 'mark'] = 'UB'
        df.loc[df['Barcode'].isin(cell_bc), 'mark'] = 'CB'
        CB_total_Genes = df.loc[df['mark'] == 'CB', 'geneID'].nunique()
        CB_reads_count = df.loc[df['mark'] == 'CB', 'count'].sum()
        reads_mapped_to_transcriptome = df['count'].sum()
        return CB_total_Genes, CB_reads_count, reads_mapped_to_transcriptome

    @runtime(f"{__name__}.collect_matrix")
    def collect_matrix(self):
        self.get_summary()

        df_line = self.get_df_line()

        line_saturation = Line_plot(df_line=df_line, title="Sequencing Saturation", x_title=X_TITLE,
                                    y_title=SATURATION_Y_TITLE, y_range=[0, 100], section=False).get_plotly_div()

        self.add_data(line_saturation=line_saturation)
        self.add_help_content(
            name='Sequencing Saturation Plot',
            content=('The plot shows the Sequencing Saturation metric as a function of downsampled '
                     'sequencing depth in mean reads per cell, up to the observed sequencing depth. '
                     'Sequencing Saturation is a measure of the observed library complexity, '
                     'and approaches 1.0 (100%) when all converted mRNA transcripts have been sequenced. '
                     'The slope of the curve near the endpoint can be interpreted as an upper bound to '
                     'the benefit to be gained from increasing the sequencing depth beyond this point')
        )

        line_median = Line_plot(df_line=df_line, title="Median Genes per Cell", x_title=X_TITLE,
                                y_title=MEDIAN_GENE_Y_TITLE).get_plotly_div()
        self.add_data(line_median=line_median)
        self.add_help_content(
            name='Median Genes per Cell Plot',
            content=('The plot shows the Median Genes per Cell as a function of downsampled '
                     'sequencing depth in mean reads per cell, up to the observed sequencing depth. '
                     'The slope of the curve near the endpoint can be interpreted as an upper bound '
                     'to the benefit to be gained from increasing the sequencing depth beyond this point.')
        )
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.marked_count_file))
        self.add_help_content(
            name='Barcode Rank Plot',
            content=('This graph shows the distribution of UMI counts in each barcode. '
                     'Barcodes can be determined to be cell-associated based on their UMI '
                     'counts or by their expression profiles. Therefore, the graph contains '
                     'both cell-associated (Colored Regions) and background-associated '
                     'barcodes (Gray Regions)')
        )

        self.add_data(qc_plot=self.qc_plot)
        self.add_help_content(
            name="nFeature_RNA",
            content='The violin plot shows gene numbers of all cell-associated barcodes.'
        )
        self.add_help_content(
            name="nCount_RNA",
            content='The violin plot shows UMI counts of all cell-associated barcodes.'
        )
        self.add_help_content(
            name="nFeature_RNA",
            content='The violin plot shows mitochondrial content of all cell-associated barcodes.'
        )

    def get_summary(self):
        CB_describe = self.CB_describe
        CB_reads_count = self.CB_reads_count
        reads_mapped_to_transcriptome = self.reads_mapped_to_transcriptome
        CB_total_Genes = self.CB_total_Genes

        estimated_cells = int(CB_describe.loc['count', 'readcount'])
        self.add_metric(
            name='Estimated Number of Cells',
            value=estimated_cells,
            help_info='The number of barcodes considered as cell-associated'
        )

        fraction_reads_in_cells = round(float(CB_reads_count) / reads_mapped_to_transcriptome * 100, 2)
        self.add_metric(
            name='Fraction Reads in Cells',
            value=fraction_reads_in_cells,
            display=f'{fraction_reads_in_cells}%',
            help_info='The fraction of valid-barcode, uniquely-mapped-to-transcriptome '
                      'reads with cell-associated barcodes'
        )

        try:
            # BarcodeSplitterRunner : _STEP_NAME = "barcode"
            valid_read_number = self.get_slot_key(
                slot='metrics',
                step_name='barcode',
                key='Valid Library Reads',
            )
        except KeyError:
            logger.warning('Will not output `Mean Reads per Cell`')
        else:
            mean_reads_per_cell = int(valid_read_number / estimated_cells)
            self.add_metric(
                name='Mean Reads per Cell',
                value=mean_reads_per_cell,
                help_info='The number of Valid reads divided by the estimated number of cells'
            )

        median_umi_per_cell = int(CB_describe.loc['50%', 'UMI'])
        self.add_metric(
            name='Median UMI counts per Cell',
            value=median_umi_per_cell,
            help_info='The median number of UMI counts per cell-associated barcode'
        )

        median_genes_per_cell = int(CB_describe.loc['50%', 'geneID'])
        self.add_metric(
            name='Median Genes per Cell',
            value=median_genes_per_cell,
            help_info='The median number of genes detected per cell-associated barcode.'
        )

        total_genes = int(CB_total_Genes)
        self.add_metric(
            name='Total Genes Detected',
            value=total_genes,
            help_info='The number of genes with at least one UMI count in any cell'
        )

        umi_saturation = round(self.down_sample_dict['umi_saturation'][-1], 2)
        read_saturation = round(self.down_sample_dict['read_saturation'][-1], 2)
        self.add_metric(
            name='Sequencing Saturation',
            value=umi_saturation,
            display=f'{umi_saturation}%',
            help_info=(
                'The fraction of reads originating from an already-observed UMI. '
                'This is a function of library complexity and sequencing depth.'
            )
        )

    @staticmethod
    def sub_sample(fraction, df_cell, cell_read_index):
        """
        umi_saturation = 1 - n_deduped_reads / n_umis
        read_saturation = 1 - n_deduped_reads / n_reads
        Currently the html report shows umi_saturation.

        n_deduped_reads = Number of unique (valid cell-barcode, valid UMI, gene) combinations among confidently mapped reads.
        n_umis = Total number of (confidently mapped, valid cell-barcode, valid UMI) UMIs.
        n_reads = Total number of (confidently mapped, valid cell-barcode, valid UMI) reads.
        Args:
            fration: subsmaple fration
            df_cell: in cell df with (Barcode geneID UMI count)
            cell_read_index: df_cell repeat index
        """
        cell_read = df_cell['count'].sum()
        frac_n_read = int(cell_read * fraction)
        subsample_read_index = cell_read_index[:frac_n_read]
        index_dedup, counts = np.unique(subsample_read_index, return_counts=True)
        n_count_once = np.sum(counts == 1)
        # total = UMI
        umi_total = len(index_dedup)
        umi_saturation = round((1 - n_count_once / umi_total) * 100, 2)
        read_total = frac_n_read
        read_saturation = round((1 - n_count_once / read_total) * 100, 2)

        # gene median
        df_cell_subsample = df_cell.loc[index_dedup,]
        geneNum_median = float(df_cell_subsample.groupby(
            'Barcode').agg({'geneID': 'nunique'}).median())

        return umi_saturation, read_saturation, geneNum_median

    @runtime(f"{__name__}.downsample")
    def downsample(self, df_cell):
        """saturation and median gene
        """
        cell_read_index = np.array(df_cell.index.repeat(df_cell['count']), dtype='int32')
        np.random.shuffle(cell_read_index)

        downsample_dict = {
            READ_FRACTION: [0],
            UMI_SATURATION: [0],
            READ_SATURATION: [0],
            MEDIAN_GENE_NUMBER: [0],
        }

        for fraction in np.arange(0.1, 1.1, 0.1):
            umi_saturation, read_saturation, geneNum_median = Count.sub_sample(
                fraction, df_cell, cell_read_index)
            fraction = round(fraction, 1)
            umi_saturation = round(umi_saturation, 2)
            read_saturation = round(read_saturation, 2)
            downsample_dict[READ_FRACTION].append(fraction)
            downsample_dict[UMI_SATURATION].append(umi_saturation)
            downsample_dict[READ_SATURATION].append(read_saturation)
            downsample_dict[MEDIAN_GENE_NUMBER].append(geneNum_median)

            self.add_metric(
                name=f'Read Fraction {fraction} read_saturation',
                value=read_saturation,
                show=False,
            )

            self.add_metric(
                name=f'Read Fraction {fraction} umi_saturation',
                value=umi_saturation,
                show=False,
            )

        df_downsample = pd.DataFrame(downsample_dict,
                                     columns=[READ_FRACTION, MEDIAN_GENE_NUMBER, UMI_SATURATION, READ_SATURATION])
        df_downsample.to_csv(self.down_sample_file, index=False, sep='\t')
        self.down_sample_dict = downsample_dict
