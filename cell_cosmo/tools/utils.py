#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : utils.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import subprocess
import logging

logger = logging.getLogger(__name__)


def sort_bam(input_bam, output_bam, threads=1):
    cmd = (
        f'samtools sort {input_bam} '
        f'-o {output_bam} '
        f'--threads {threads} '
    )
    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


def index_bam(input_bam, samtools_index_param=None):
    if samtools_index_param is None:
        samtools_index_param = ""  # -c -m 4
    cmd = f"samtools index {samtools_index_param} {input_bam}"
    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
