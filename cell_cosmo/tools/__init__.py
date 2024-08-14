#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : __init__.py.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : 第三方软件调用,软件开发功能调用
"""
from cell_cosmo.tools.Samtools import Samtools

# count
OUTS_DIR = 'outs'
RAW_MATRIX_DIR_SUFFIX = 'raw'
FILTERED_MATRIX_DIR_SUFFIX = 'filtered'
MATRIX_FILE_NAME = 'matrix.mtx.gz'
FEATURE_FILE_NAME = 'features.tsv.gz'
BARCODE_FILE_NAME = 'barcodes.tsv.gz'
STAR_BAM_SUFFIX = 'Aligned.out.bam'
TAG_BAM_SUFFIX = 'aligned_posSorted_addTag.bam'
STARSOLO_BAM_SUFFIX = 'Aligned.sortedByCoord.out.bam'
COUNTS_FILE_NAME = 'counts.tsv'

__all__ = (
    'Samtools'
)
