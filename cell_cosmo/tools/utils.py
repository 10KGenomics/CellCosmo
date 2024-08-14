#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : utils.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import os
import logging
import subprocess
import pandas as pd

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


def get_matrix_file_path(matrix_dir, file_name):
    """
    compatible with non-gzip file
    """
    non_gzip = file_name.strip('.gz')
    file_path_list = [f'{matrix_dir}/{file_name}', f'{matrix_dir}/{non_gzip}']
    for file_path in file_path_list:
        if os.path.exists(file_path):
            return file_path


def read_one_col(file):
    """
    Read file with one column. Strip each line.
    Returns col_list, line number
    """
    df = pd.read_csv(file, header=None)
    col1 = list(df.iloc[:, 0])
    col1 = [item.strip() for item in col1]
    num = len(col1)
    return col1, num


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")
