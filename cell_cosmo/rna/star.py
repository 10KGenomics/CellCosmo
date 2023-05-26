#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : star.py
@Time       : 2022/06/08
@Version    : 1.0
@Desc       : None
"""
import re
import logging
import pandas as pd
from cell_cosmo.output_runner.base_star_runner import BaseSTARRunner
from cell_cosmo.util import runtime, GenomeUtil
from cell_cosmo.tools.plotlyplot import Pie_plot

logger = logging.getLogger(__name__)


class Star(BaseSTARRunner):
    _STEP_NAME = "star"
    _DISPLAY_TITLE = "Mapping"

    def __init__(self, picard_mem, **kwargs):
        super(Star, self).__init__(**kwargs)
        self.picard_mem = picard_mem
        # reset parent attr
        self.genome = GenomeUtil.parse_rna_dir(self.genomeDir)
        self.refflat = self.genome['refflat']

        self.ribo_log = f'{self.out_prefix}_ribo_log.txt'
        self.ribo_run_log = f'{self.out_prefix}_ribo_run.log'
        self.picard_region_log = f'{self.out_prefix}_region.log'
        self.plot = None
        # self.stats = pd.Series()
        self.df_region = pd.DataFrame()

    @runtime(__name__)
    def run(self):
        super().run()
        self.picard()

    # def ribo(self):
    #     # TODO remove bbduk.sh and use picard ribo bases
    #     human_ribo_fa = f'{ROOT_PATH}/data/rRNA/human_ribo.fasta'
    #     cmd = (
    #         f'bbduk.sh '
    #         f'in1={self.fq} '
    #         f'ref={human_ribo_fa} '
    #         f'stats={self.ribo_log} '
    #         f'overwrite=t '
    #         f'> {self.ribo_run_log} 2>&1 '
    #     )
    #     self.debug_subprocess_call(cmd)

    @runtime(f"{__name__}.picard")
    def picard(self):
        cmd = [
            'picard',
            f'-Xmx{self.picard_mem}G',
            f'-XX:ParallelGCThreads={str(self.thread)}',
            'CollectRnaSeqMetrics',
            'I=%s' % self.STAR_bam,
            'O=%s' % self.picard_region_log,
            'REF_FLAT=%s' % self.refflat,
            'STRAND=NONE',
            'VALIDATION_STRINGENCY=SILENT']
        cmd = ' '.join(cmd)
        logger.info(cmd)
        self.debug_subprocess_call(cmd)

    def collect_matrix(self):
        """
        add picard region bases
        add region plot
        if debug, add ribosomal RNA reads percent
        """
        super().collect_matrix()

        with open(self.picard_region_log, 'r') as picard_log:
            region_dict = {}
            for line in picard_log:
                if not line:
                    break
                if line.startswith('## METRICS CLASS'):
                    header = picard_log.readline().strip().split('\t')
                    data = picard_log.readline().strip().split('\t')
                    region_dict = dict(zip(header, data))
                    break

        total = float(region_dict['PF_ALIGNED_BASES'])
        exonic_regions = int(region_dict['UTR_BASES']) + \
                         int(region_dict['CODING_BASES'])
        intronic_regions = int(region_dict['INTRONIC_BASES'])
        intergenic_regions = int(region_dict['INTERGENIC_BASES'])

        self.add_metric(
            name='Read Bases Mapped to Exonic',
            value=exonic_regions,
            total=total,
            help_info='Bases that mapped uniquely to a coding base or a UTR base of the genome'
        )
        self.add_metric(
            name='Read Bases Mapped to Intronic',
            value=intronic_regions,
            total=total,
            help_info='Bases that mapped uniquely to an intronic base of the genome, and not a coding or UTR base'
        )
        self.add_metric(
            name='Read Bases Mapped to Intergenic',
            value=intergenic_regions,
            total=total,
            help_info='Bases that mapped uniquely to an intergenic base of the genome, not align to any gene'
        )

        # ribo
        """
        if self.debug:
            with open(self.ribo_log, 'r') as ribo_log:
                for line in ribo_log:
                    if line.find('#Matched') != -1:
                        items = line.split()
                        Reads_Mapped_to_rRNA = int(items[1])
                    if line.find('#Total') != -1:
                        items = line.split()
                        Reads_Total = int(items[1])
                self.add_metric(
                    name=f'{self.stat_prefix} Mapped to rRNA',
                    value=Reads_Mapped_to_rRNA,
                    total=Reads_Total,
                    help_info='Number of reads or umis that mapped to rRNA'
                )
        """

        region_plot = {'regions': ['Exonic Regions', 'Intronic Regions', 'Intergenic Regions'],
                       'values': [exonic_regions, intronic_regions, intergenic_regions]}
        self.df_region = pd.DataFrame(region_plot)
        region_pie = Pie_plot(df_region=self.df_region).get_plotly_div()
        self.add_data(region_pie=region_pie)
