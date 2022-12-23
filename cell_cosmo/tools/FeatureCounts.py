#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : FeatureCounts.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import os
import re
from cell_cosmo.output_runner import BaseReportRunner
from cell_cosmo.util import GenomeUtil, runtime
from cell_cosmo.tools import utils, Samtools


class FeatureCounts(BaseReportRunner):
    _STEP_NAME = "featureCounts"
    # _DISPLAY_TITLE = "featureCounts"
    _DISPLAY_TITLE = "Quantification"

    def __init__(self, genomeDir, input_file, gtf_type, featureCounts_param, **kwargs):
        super(FeatureCounts, self).__init__(**kwargs)
        self.gtf_type = gtf_type
        self.input = input_file
        # set
        self.gtf = GenomeUtil.parse_rna_dir(genomeDir)['gtf']
        self.featureCounts_param = featureCounts_param

        # out files
        input_basename = os.path.basename(self.input)
        self.featureCounts_bam = f'{self.outdir}/{input_basename}.featureCounts.bam'
        self.name_sorted_bam = f'{self.out_prefix}_name_sorted.bam'
        self.featureCount_log_file = f'{self.out_prefix}.summary'

    def run_featureCounts(self):
        cmd = (
            'featureCounts '
            '-s 1 '
            f'-a {self.gtf} '
            f'-o {self.out_prefix} '  # not bam
            '-R BAM '
            f'-T {self.thread} '
            f'-t {self.gtf_type} '
            f'{self.input} '
        )
        if self.featureCounts_param:
            cmd += (" " + self.featureCounts_param)
        self.debug_subprocess_call(cmd)

    @runtime(__name__)
    def run(self):
        self.run_featureCounts()
        samtools_runner = Samtools(
            in_bam=self.featureCounts_bam,
            out_bam=self.featureCounts_bam,
            threads=self.thread,
            debug=self.debug
        )
        samtools_runner.add_tag(self.gtf)
        samtools_runner.temp_sam2bam(by='coord')
        samtools_runner.samtools_sort(
            in_file=self.featureCounts_bam,
            out_file=self.name_sorted_bam,
            by='name',
        )

    def collect_matrix(self):
        metrics_strs = ['Assigned', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity']
        metrics_numbers = {}
        metrics_compiled = {}

        for metrics_str in metrics_strs:
            raw_str = re.escape(metrics_str) + r'.*?(\d+)'
            compiled = re.compile(raw_str, flags=re.S)
            metrics_compiled[metrics_str] = compiled

        with open(self.featureCount_log_file, 'r') as fh:

            for line in fh:
                line = line.strip()
                if not line:
                    continue

                for metrics_str in metrics_compiled:
                    compiled = metrics_compiled[metrics_str]
                    match = compiled.search(line)
                    if match:
                        metrics_numbers[metrics_str] = int(match.group(1))
                        break

            total = sum(metrics_numbers.values())

            self.add_metric(
                name='Assigned',
                value=metrics_numbers['Assigned'],
                total=total,
                help_info='Reads that can be successfully assigned without ambiguity'
            )
            self.add_metric(
                name='Unassigned_NoFeatures',
                value=metrics_numbers['Unassigned_NoFeatures'],
                total=total,
                help_info='Alignments that do not overlap any feature'
            )
            self.add_metric(
                name='Unassigned_Ambiguity',
                value=metrics_numbers['Unassigned_Ambiguity'],
                total=total,
                help_info='Alignments that overlap two or more features'
            )
