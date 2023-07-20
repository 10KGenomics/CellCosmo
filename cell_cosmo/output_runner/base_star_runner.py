#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : base_star_runner.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.util import runtime
from cell_cosmo.tools import utils
import re
from cell_cosmo.util import GenomeUtil
from cell_cosmo.output_runner import BaseReportRunner
import subprocess
import logging

logger = logging.getLogger(__name__)


class BaseSTARRunner(BaseReportRunner):
    """
    base class for STAR
    """
    _STEP_NAME = None
    _DISPLAY_TITLE = None

    def __init__(self, genomeDir, fq, out_unmapped,
                 outFilterMultimapNmax, outFilterMatchNmin, consensus_fq,
                 STAR_param, add_prefix=None, **kwargs):
        super(BaseSTARRunner, self).__init__(**kwargs)

        self.fq = fq
        self.genomeDir = genomeDir
        self.out_unmapped = out_unmapped

        self.outFilterMatchNmin = int(outFilterMatchNmin)
        self.multi_max = int(outFilterMultimapNmax)
        self.STAR_param = STAR_param
        self.consensus_fq = consensus_fq

        # parse
        self.genome = GenomeUtil.parse_dir(self.genomeDir)
        self.stat_prefix = 'Reads'
        if self.consensus_fq:
            self.stat_prefix = 'UMIs'

        # out
        self.outPrefix = f'{self.outdir}/{self.sample}_'
        if add_prefix:
            self.outPrefix += add_prefix + '_'
        self.STAR_map_log = f'{self.outPrefix}Log.final.out'
        self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
        self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'

    @runtime(f'{__name__}.STAR')
    def STAR(self):
        cmd = [
            'STAR',
            '--runThreadN', str(self.thread),
            '--genomeDir', self.genomeDir,
            '--readFilesIn', self.fq,
            '--outFilterMultimapNmax', str(self.multi_max),
            '--outFileNamePrefix', self.outPrefix,
            '--outSAMtype', 'BAM', 'Unsorted',  # controls sort by Coordinate or not
            '--outFilterMatchNmin', str(self.outFilterMatchNmin)
        ]
        if self.out_unmapped:
            cmd += ['--outReadsUnmapped', 'Fastx']
        if self.fq[-3:] == ".gz":
            # 苹果电脑 zcat 异常
            cmd += ['--readFilesCommand', 'zcat']
        cmd = ' '.join(cmd)
        if self.STAR_param:
            # STAR_param 异常添加双引号或单引号的处理
            self.STAR_param = str(self.STAR_param).lstrip(
                "'").lstrip('"').rstrip("'").rstrip('"')
            cmd += (" " + self.STAR_param)
        logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.STAR()
        # TODO 数据收集
        # self.get_star_metrics()
        self.sort_bam()
        self.index_bam()

    @runtime(f'{__name__}.sort_bam')
    def sort_bam(self):
        utils.sort_bam(
            self.unsort_STAR_bam,
            self.STAR_bam,
            threads=self.thread,
        )

    @runtime(f'{__name__}.index_bam')
    def index_bam(self):
        utils.index_bam(self.STAR_bam)

    def collect_matrix(self):
        """collect matrix"""

        with open(self.STAR_map_log, 'r') as map_log:
            # number amd percent
            unique_reads_list = []
            multi_reads_list = []
            total_reads = 0
            for line in map_log:
                if line.strip() == '':
                    continue
                if re.search(r'Uniquely mapped reads', line):
                    unique_reads_list.append(line.strip().split()[-1])
                if re.search(r'of reads mapped to too many loci', line):
                    multi_reads_list.append(line.strip().split()[-1])
                if re.search(r'Number of input reads', line):
                    total_reads = int(line.strip().split()[-1])

        unique_reads = int(unique_reads_list[0])
        multi_reads = int(multi_reads_list[0])

        self.add_metric(
            name='Genome',
            value=self.genome['genome_name'],
        )
        self.add_metric(
            name=f'Uniquely Mapped {self.stat_prefix}',
            value=unique_reads,
            total=total_reads,
            # help_info='reads that mapped uniquely to the genome'
            help_info='Fraction of reads that mapped uniquely to the genome'
        )
        self.add_metric(
            name=f'Multi-Mapped {self.stat_prefix}',
            value=multi_reads,
            total=total_reads,
            # help_info='reads that mapped to multiple locations in the genome'
            help_info='Fraction of Reads that mapped to multiple locations in the genome'
        )
