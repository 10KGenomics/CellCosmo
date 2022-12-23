#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : Sample.py
@Time       : 2022/06/10
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.output_runner.base_report_runner import BaseReportRunner
from cell_cosmo.util import runtime
from cell_cosmo import __VERSION__


class Sample(BaseReportRunner):
    _STEP_NAME = "sample"
    _DISPLAY_TITLE = "Sample"

    def __init__(self, chemistry, **kwargs):
        super(Sample, self).__init__(**kwargs)
        self.assay_description = self.assay_text
        self.version = __VERSION__
        self.chemistry = chemistry

    @runtime(__name__)
    def run(self):
        # TODO,如果为默认的auto,则通过fq1文件自动识别chemistry
        if self.chemistry == 'auto':
            chemistry = self.chemistry
        else:
            chemistry = self.chemistry
        self.chemistry = chemistry

    def collect_matrix(self):
        self.add_metric(
            name='Sample ID',
            value=self.sample,
            help_info="Sample name."
        )
        self.add_metric(
            name='Omics',
            value=self.assay_description,
            help_info="Types of omics technology."
        )
        self.add_metric(
            name='Chemistry',
            value=self.chemistry,
            display=f"{self.chemistry}",
            help_info='Chemistry version.',
        )
        self.add_metric(
            name='Software Version',
            value=self.version,
        )
