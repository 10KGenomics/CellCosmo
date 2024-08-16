#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : star_solo.py
@Time       : 2024/07/09
@Version    : 1.0
@Desc       : None
"""
import os
import re
import sys
import gzip
import pysam
import logging
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
from abc import ABCMeta
from collections import Counter

from anndata import AnnData

from cell_cosmo.tools import utils
from cell_cosmo.output_runner.base_runner import BaseRunner
from cell_cosmo.output_runner import BaseReportRunner
from cell_cosmo.util import runtime, GenomeUtil
from cell_cosmo.tools.Count import Count
from cell_cosmo.tools.plotlyplot import Pie_plot
from cell_cosmo.tools.aa import get_plot_elements
from cell_cosmo.tools.chemistry import LibraryInfo, get_sequence_by_pattern
from cell_cosmo.tools import FILTERED_MATRIX_DIR_SUFFIX, COUNTS_FILE_NAME
from cell_cosmo.tools.matrix import CountMatrix
from cell_cosmo.tools.plotlyplot import Line_plot, Violin_plot
from cell_cosmo.util import reader

logger = logging.getLogger(__name__)
# downsample.csv
READ_FRACTION = 'read_fraction'
MEDIAN_GENE_NUMBER = 'median_gene_number'
READ_SATURATION = 'read_saturation'
UMI_SATURATION = 'umi_saturation'

# Plot axis title in HTML
X_TITLE = 'Read Fraction'
SATURATION_Y_TITLE = 'Sequencing Saturation(%)'
MEDIAN_GENE_Y_TITLE = 'Median Genes per Cell'


def get_solo_pattern(library_info: LibraryInfo) -> (str, str, str):
    """
    Returns:
        solo_type
        cb_str
        umi_str
    """
    if len(library_info.pattern_U) != 1:
        raise Exception(f"Error: Wrong pattern:{library_info.pattern_str}. \n "
                        f"Solution: fix pattern so that UMI only have 1 position.\n")

    ul, ur = library_info.pattern_U[0]
    umi_len = ur - ul

    if len(library_info.pattern_C) == 1:
        solo_type = 'CB_UMI_Simple'
        l, r = library_info.pattern_C[0]
        cb_start = l + 1
        cb_len = r - l
        umi_start = ul + 1
        cb_str = f'--soloCBstart {cb_start} --soloCBlen {cb_len} '
        umi_str = f'--soloUMIstart {umi_start} --soloUMIlen {umi_len} '
    else:
        solo_type = 'CB_UMI_Complex'
        cb_pos = ' '.join([f'0_{s}_0_{e - 1}' for s, e in library_info.pattern_C])
        umi_pos = f'0_{ul}_0_{ur - 1}'
        cb_str = f'--soloCBposition {cb_pos} '
        umi_str = f'--soloUMIposition {umi_pos} --soloUMIlen {umi_len} '
    return solo_type, cb_str, umi_str


class _O(BaseReportRunner, metaclass=ABCMeta):
    def __init__(self, **kwargs):
        super(_O, self).__init__(**kwargs)
        # move file to outs
        self.outs = []
        self.outs_dir = f'{self.outdir}/../outs'
        utils.check_mkdir(self.outdir)
        utils.check_mkdir(self.outs_dir)

    def _move_files(self):
        for f in self.outs:
            cmd = ''
            if not os.path.exists(f):
                sys.stderr.write(f'WARNING: output file {f} not found! The pipeline may have failed.\n')
                continue
            elif os.path.isfile(f):
                cmd = f'mv -f {f} {self.outs_dir}'
            elif os.path.isdir(f):
                cmd = f'set -e; cp -r {f} {self.outs_dir}; rm -r {f}'
            subprocess.check_call(cmd, shell=True)

    def _clean_up(self):
        # 重写父类,调整 新增执行步骤
        self._move_files()
        self._add_content_data()
        self._add_content_metric()
        self.process_starsolo_demultiplexing_step()
        self.set_summary_step()
        self._write_stat()
        self._dump_content()
        self._render_html()
        self._mtx2tsv()


class Starsolo(_O):
    _STEP_NAME = "star"
    _DISPLAY_TITLE = "Mapping"

    def __init__(self, fq1, fq2, genomeDir, out_unmapped,
                 outFilterMatchNmin, consensus_fq, STAR_param,
                 pattern=None, chemistry_config=None, chemistry_name=None, **kwargs):
        # 删除了 outFilterMultimapNmax 参数
        super(Starsolo, self).__init__(**kwargs)

        self.fq1 = fq1
        self.fq2 = fq2
        self.library_info = LibraryInfo(
            chemistry_name, chemistry_config, pattern,
            auto_init_db=False
        )
        self.solo_type, self.cb_str, self.umi_str = get_solo_pattern(self.library_info)
        self.whitelist_str = " ".join(self.library_info.barcode_files)

        self.genomeDir = genomeDir
        self.out_unmapped = out_unmapped

        self.soloFeatures = kwargs.get("soloFeatures")
        self.SAM_attributes = kwargs.get("SAM_attributes")
        self.soloCellFilter = kwargs.get("soloCellFilter")
        self.outSAMtype = kwargs.get("outSAMtype")

        self.outFilterMatchNmin = int(outFilterMatchNmin)
        # self.multi_max = int(outFilterMultimapNmax)
        self.STAR_param = STAR_param
        self.consensus_fq = consensus_fq

        # parse
        self.genome = GenomeUtil.parse_dir(self.genomeDir)
        self.stat_prefix = 'Reads'
        if self.consensus_fq:
            self.stat_prefix = 'UMIs'

        # out
        self.outPrefix = f'{self.outdir}/{self.sample}_'
        self.STAR_map_log = f'{self.outPrefix}Log.final.out'
        self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
        self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'

        # output files
        self.solo_out_dir = f'{self.outdir}/{self.sample}_Solo.out/'
        solo_dir = f'{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS'
        self.raw_matrix = f'{solo_dir}/raw'
        self.filtered_matrix = f'{solo_dir}/filtered'
        self.summary_file = f'{solo_dir}/Summary.csv'
        self.outs = [self.raw_matrix, self.filtered_matrix, self.STAR_bam]

    def STAR(self):
        cmd = [
            'STAR',
            '--genomeDir', self.genomeDir,
            '--readFilesIn', self.fq2 + ' ' + self.fq1,
            f'--soloCBwhitelist {self.whitelist_str}',
            f'--soloCellFilter {self.soloCellFilter}',
            '--outFileNamePrefix', self.outPrefix,
            '--runThreadN', str(self.thread),
            '--outFilterMatchNmin', str(self.outFilterMatchNmin),
            f'--soloFeatures {self.soloFeatures}',
            f'--soloType {self.solo_type}',
            self.cb_str,  # https://github.com/alexdobin/STAR/issues/1607
            self.umi_str,
            # length of the barcode read 1 equal to sum of soloCBlen+soloUMIlen 0 not defined, do not check
            # only one match in whitelist with 1 mismatched base allowed
            '--soloCBmatchWLtype 1MM',
            # cell reads stats https://github.com/alexdobin/STAR/issues/1501
            '--soloCellReadStats Standard',
            '--soloBarcodeReadLength 0',
            # '--outFilterMultimapNmax', str(self.multi_max),
            # '--soloBarcodeMate', '0',
            # '--clip5pNbases','74 0',
        ]
        if self.outSAMtype is not None:
            cmd.append(f"--outSAMtype {self.outSAMtype}")
            if self.SAM_attributes is not None:
                cmd.append(f"--outSAMattributes {self.SAM_attributes}")
        else:
            cmd.append(f"--outSAMtype None")

        if self.fq1[-3:] == ".gz":
            # 苹果电脑 zcat 异常
            cmd += ['--readFilesCommand', 'zcat']
        # if self.out_unmapped:
        #     cmd += ['--outReadsUnmapped', 'Fastx']
        cmd = ' '.join(cmd)
        if self.STAR_param:
            # STAR_param 异常添加双引号或单引号的处理
            self.STAR_param = str(self.STAR_param).lstrip(
                "'").lstrip('"').rstrip("'").rstrip('"')
            cmd += (" " + self.STAR_param)
        logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def gzip_matrix(self):
        cmd = f'gzip {self.raw_matrix}/*; gzip {self.filtered_matrix}/*'
        subprocess.check_call(cmd, shell=True)

    @runtime(__name__)
    def run(self):
        self.STAR()
        self.gzip_matrix()
        pass

    def collect_matrix(self):
        super().collect_matrix()

    def get_Q30_cb_UMI(self):
        fq1_list = self.fq1.split(",")

        cb_10k, umi_10k, cb_10k_100k, umi_10k_100k = Counter(), Counter(), Counter(), Counter()
        n = 0
        with pysam.FastxFile(fq1_list[0], persist=False) as fq1:
            for entry in fq1:
                n += 1
                if n > 10 ** 5:
                    break
                qual = entry.quality
                cb_qual = get_sequence_by_pattern(qual, self.library_info.pattern_C)
                umi_qual = get_sequence_by_pattern(qual, self.library_info.pattern_U)
                if n <= 10 ** 4:
                    cb_10k.update(cb_qual)
                    umi_10k.update(umi_qual)
                else:
                    cb_10k_100k.update(cb_qual)
                    umi_10k_100k.update(umi_qual)

        cb_qual_counter = cb_10k
        umi_qual_counter = umi_10k
        if cb_10k_100k:
            cb_qual_counter = cb_10k_100k
            umi_qual_counter = umi_10k_100k

        q30_cb = sum([cb_qual_counter[k] for k in cb_qual_counter
                      if k >= chr(30 + 33)]) / float(sum(cb_qual_counter.values()))
        q30_umi = sum([umi_qual_counter[k] for k in umi_qual_counter
                       if k >= chr(30 + 33)]) / float(sum(umi_qual_counter.values()))
        return q30_cb, q30_umi


class Mapping(_O):
    _STEP_NAME = "star"
    _DISPLAY_TITLE = "Mapping"

    def __init__(self, genomeDir, **kwargs):
        super(Mapping, self).__init__(**kwargs)
        solo_dir = f'{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS'

        self.cellReadsStats = f'{solo_dir}/CellReads.stats'

        self.filtered_matrix = f'{self.outs_dir}/{FILTERED_MATRIX_DIR_SUFFIX}'
        self.counts_file = f'{solo_dir}/{COUNTS_FILE_NAME}'
        self.outs = [self.counts_file]
        self.genome = GenomeUtil.parse_rna_dir(genomeDir)['genome_name']
        # wait to set
        self.valid = None
        self.corrected = None

    def run(self):
        df = pd.read_csv(self.cellReadsStats, sep='\t', header=0, index_col=0)
        df = df.iloc[1:, ]  # skip first line cb not pass whitelist
        df_count = df.loc[:, ['nUMIunique', 'countedU']]  # keep dataframe format
        df_count.rename(columns={'nUMIunique': 'UMI'}, inplace=True)

        cols = [
            'cbMatch', 'cbPerfect', 'genomeU', 'genomeM', 'exonic',
            'intronic', 'exonicAS', 'intronicAS', 'countedU'
        ]
        df = df.loc[:, cols]
        s = df.sum()

        # json does not recognize NumPy data types.
        # TypeError: Object of type int64 is not JSON serializable
        valid = int(s['cbMatch'])
        perfect = int(s['cbPerfect'])
        corrected = valid - perfect
        self.valid = valid
        self.corrected = corrected
        genomeU = int(s['genomeU'])
        genomeM = int(s['genomeM'])
        mapped = genomeU + genomeM
        exonic = int(s['exonic'])
        intronic = int(s['intronic'])
        antisense = int(s['exonicAS'] + s['intronicAS'])
        intergenic = mapped - exonic - intronic - antisense
        countedU = int(s['countedU'])
        del df
        self.add_metric(
            name='genome',
            value=self.genome,
        )

        self.add_metric(
            name='Uniquely Mapped Reads',
            value=genomeU,
            total=valid,
            help_info='Fraction of reads that mapped uniquely to the genome.'
        )

        self.add_metric(
            name='Multi-Mapped Reads',
            value=genomeM,
            total=valid,
            help_info='Fraction of Reads that mapped to multiple locations in the genome.'
        )
        unique_transcriptome = countedU / valid
        # self.add_metric(
        #     name='Reads mapped uniquely to Transcriptome',
        #     value=countedU,
        #     total=valid,
        #     help_info='Reads that mapped to a unique gene in the transcriptome. '
        #               'These reads are used for UMI counting.'
        # )
        self.add_metric(
            name='Mapped Reads Assigned To Exonic Regions',
            value=exonic,
            total=mapped,
            help_info='Bases that mapped uniquely to a coding base or a UTR base of the genome.',
        )
        self.add_metric(
            name='Mapped Reads Assigned To Intronic Regions',
            value=intronic,
            total=mapped,
            help_info='Bases that mapped uniquely to an intronic base of the genome, and not a coding or UTR base.',
        )
        self.add_metric(
            name='Mapped Reads Assigned To Intergenic Regions',
            value=intergenic,
            total=mapped,
            help_info='Bases that mapped uniquely to an intergenic base of the genome, not align to any gene.',
        )
        self.add_metric(
            name='Mapped Reads Assigned Antisense To Gene',
            value=antisense,
            total=mapped,
            help_info='Reads that assigned to the opposite strand of genes',
        )
        df_count.sort_values(by='UMI', ascending=False, inplace=True)
        cbs = CountMatrix.read_barcodes(self.filtered_matrix)
        df_count['mark'] = 'UB'
        for cb in cbs:
            df_count.loc[cb, 'mark'] = 'CB'
        df_count.to_csv(self.counts_file, sep='\t', index=True)

        region_plot = {'regions': ['Exonic Regions', 'Intronic Regions', 'Intergenic Regions'],
                       'values': [exonic, intronic, intergenic]}
        df_region = pd.DataFrame(region_plot)
        region_pie = Pie_plot(df_region=df_region).get_plotly_div()
        self.add_data(region_pie=region_pie)

    def get_vc(self):
        return self.valid, self.corrected

    def collect_matrix(self):
        pass


class Cells(_O):
    _STEP_NAME = "count"
    _DISPLAY_TITLE = "Cells"

    def __init__(self, valid_reads, genomeDir, **kwargs):
        super(Cells, self).__init__(**kwargs)
        solo_dir = f'{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS'
        self.summary_file = f'{solo_dir}/Summary.csv'

        self.counts_file = f'{self.outs_dir}/{COUNTS_FILE_NAME}'
        self.cell_matrix_dir = f'{self.outs_dir}/filtered'
        self.cell_matrix_raw_dir = f'{self.outs_dir}/raw'
        self.mt_gene_list = GenomeUtil.parse_rna_dir(genomeDir)['mt_gene_list']
        self.valid_reads = valid_reads
        self._parse_summary_add_metrics()

    def _parse_summary_add_metrics(self):
        df = pd.read_csv(self.summary_file, index_col=0, header=None)
        s = df.iloc[:, 0]
        self.n_cells = int(s["Estimated Number of Cells"])
        self.mean_used_reads_per_cell = int(s["Mean Reads per Cell"])
        self.median_umi_per_cell = int(s["Median UMI per Cell"])
        self.median_genes_per_cell = int(s["Median GeneFull_Ex50pAS per Cell"])
        self.total_genes = int(s["Total GeneFull_Ex50pAS Detected"])
        self.saturation = float(s["Sequencing Saturation"])
        # 3 args
        self.fraction_reads_in_cells = float(s["Fraction of Unique Reads in Cells"])
        self.n_reads = int(s["Number of Reads"])
        self.q30_RNA = float(s["Q30 Bases in RNA read"])

    def run(self):
        pass

    def _get_plots(self):
        adata = sc.read_10x_mtx(
            self.cell_matrix_dir,
            var_names='gene_symbols',
        )
        qc_plot = self._get_qc_plot(adata)
        # df_line = self._get_df_line(adata)
        df_line = None
        return qc_plot, df_line

    def _get_qc_plot(self, adata: AnnData):
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
        qc_plot = {
            "option": [],
            "plot_data": [],
        }
        qc_plot_name_map = {
            'n_genes_by_counts': 'nFeature_RNA',
            'total_counts': 'nCount_RNA',
            'pct_counts_mito': 'percent.mt'
        }
        for key, rename_key in qc_plot_name_map.items():
            df_plot = adata.obs[[key]].copy()
            df_plot.columns = [rename_key]
            data = Violin_plot(df_plot, rename_key).get_plotly_div()
            class_name = f"qc_plot_{rename_key.replace('.', '_')}"
            qc_plot['option'].append({
                'title': rename_key,
                'class_name': class_name
            })
            qc_plot['plot_data'].append({
                'data': data,
                'class_name': class_name
            })
        return qc_plot

    def _down_sample(self, df_cell):
        """saturation and median gene"""
        cell_read_index = np.array(df_cell.index.repeat(df_cell['count']), dtype='int32')
        np.random.shuffle(cell_read_index)

        down_sample_dict = {
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
            down_sample_dict[READ_FRACTION].append(fraction)
            down_sample_dict[UMI_SATURATION].append(umi_saturation)
            down_sample_dict[READ_SATURATION].append(read_saturation)
            down_sample_dict[MEDIAN_GENE_NUMBER].append(geneNum_median)

        return pd.DataFrame(down_sample_dict, columns=[
            READ_FRACTION, MEDIAN_GENE_NUMBER,
            UMI_SATURATION, READ_SATURATION])

    def _get_df_line(self, adata):
        # df = adata.to_df()
        # df = df.melt(var_name='geneID', value_name='count', ignore_index=False)
        # df = df.loc[df['count'] != 0, :]
        # df = df.reset_index(names=['Barcode'])

        def read_mtx(filepath):
            with gzip.open(filepath, 'rt') as fh:
                # 过滤三行
                fh.readline()
                fh.readline()
                fh.readline()
                df_ = pd.DataFrame(
                    data=[str(line).strip().split(' ') for line in fh.readlines()],
                    columns=['geneID', 'Barcode', 'count']
                )
                df_["count"] = df_["count"].astype(int)
                return df_

        df = read_mtx(os.path.join(self.cell_matrix_raw_dir, 'matrix.mtx.gz'))

        # down sample
        df_line = self._down_sample(df)
        df_line.rename(columns={
            READ_FRACTION: X_TITLE,
            MEDIAN_GENE_NUMBER: MEDIAN_GENE_Y_TITLE,
            UMI_SATURATION: SATURATION_Y_TITLE,
        }, inplace=True)

        return df_line

    def collect_matrix(self):
        qc_plot, df_line = self._get_plots()

        if df_line:
            self.saturation = float(df_line.loc[df_line[X_TITLE] == 1.0, SATURATION_Y_TITLE]) / 100
            self.median_genes_per_cell = int(df_line.loc[df_line[X_TITLE] == 1.0, MEDIAN_GENE_Y_TITLE])

        self.add_metric(
            name='Estimated Number of Cells',
            value=self.n_cells,
            help_info='The number of barcodes considered as cell-associated'
        )

        self.add_metric(
            name='Fraction Reads in Cells',
            value=self.fraction_reads_in_cells,
            display='{:.2%}'.format(self.fraction_reads_in_cells),
            help_info='The fraction of valid-barcode, uniquely-mapped-to-transcriptome '
                      'reads with cell-associated barcodes'
        )

        mean_reads_per_cell = int(self.valid_reads / self.n_cells)
        self.add_metric(
            name='Mean Reads per Cell',
            value=mean_reads_per_cell,
            help_info='The number of Valid reads divided by the estimated number of cells'
        )

        self.add_metric(
            name='Median UMI counts per Cell',
            value=self.median_umi_per_cell,
            help_info='The median number of UMI counts per cell-associated barcode'
        )

        self.add_metric(
            name='Median Genes per Cell',
            value=self.median_genes_per_cell,
            help_info='The median number of genes detected per cell-associated barcode.'
        )

        self.add_metric(
            name='Total Genes Detected',
            value=self.total_genes,
            help_info='The number of genes with at least one UMI count in any cell'
        )

        self.add_metric(
            name='Sequencing Saturation',
            value=self.saturation,
            display='{:.2%}'.format(self.saturation),
            help_info=(
                'The fraction of reads originating from an already-observed UMI. '
                'This is a function of library complexity and sequencing depth.'
            )
        )

        if df_line:
            line_saturation = Line_plot(
                df_line=df_line, title="Sequencing Saturation", x_title=X_TITLE,
                y_title=SATURATION_Y_TITLE, y_range=[0, 100], section=False
            ).get_plotly_div()

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

        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.counts_file))
        self.add_help_content(
            name='Barcode Rank Plot',
            content=('This graph shows the distribution of UMI counts in each barcode. '
                     'Barcodes can be determined to be cell-associated based on their UMI '
                     'counts or by their expression profiles. Therefore, the graph contains '
                     'both cell-associated (Colored Regions) and background-associated '
                     'barcodes (Gray Regions)')
        )
        self.add_data(qc_plot=qc_plot)
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


class Demultiplexing(_O):
    _STEP_NAME = "barcode"
    _DISPLAY_TITLE = "Demultiplexing"

    def __init__(self, valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_rna, **kwargs):
        super(Demultiplexing, self).__init__(**kwargs)
        self.valid_reads = valid_reads
        self.n_reads = n_reads
        self.corrected = corrected
        self.q30_cb = q30_cb
        self.q30_umi = q30_umi
        self.q30_RNA = q30_rna

    def run(self):
        pass

    def collect_matrix(self):
        self.add_metric(
            name='Number of Raw Reads',
            value=self.n_reads,
            help_info='Total number of read pairs from Fastq file'
        )
        self.add_metric(
            name='Valid Library Reads',
            value=self.valid_reads,
            total=self.n_reads,
            help_info='Fraction of reads that conform to the expected library structure'
        )
        self.add_metric(
            name='Corrected Barcodes',
            value=self.corrected,
            total=self.valid_reads,
            help_info='fraction of valid reads with corrected barcodes. '
                      'Barcodes are corrected to the whitelist sequence that is 1 Hamming-distance away.'
        )
        self.add_metric(
            name='Q30 of Barcode',
            value=self.q30_cb,
            display='{:.2%}'.format(self.q30_cb),
            help_info='Fraction of barcode bases with quality scores over Q30',
        )

        self.add_metric(
            name='Q30 of UMI',
            value=self.q30_umi,
            display='{:.2%}'.format(self.q30_umi),
            help_info='Fraction of UMI bases with quality scores over Q30',
        )

        self.add_metric(
            name='Q30 of RNA Read',
            value=self.q30_RNA,
            display='{:.2%}'.format(self.q30_RNA),
            help_info='Fraction of RNA read bases with quality scores over Q30',
        )
