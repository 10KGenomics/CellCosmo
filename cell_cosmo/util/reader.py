#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : reader.py
@Time       : 2022/05/30
@Version    : 1.0
@Desc       : None
"""
import os
import gzip
from xopen import xopen
# from cell_cosmo.constant import FileType
import itertools


def reader(filepath: str, ignore_test_env=False):
    """
    1. 该方法会自动通过后缀识别是否为压缩格式,目前仅支持识别gz压缩格式
    2. 该方法返回tuple形式的数据,
        - 对 .fa|.fasta 文件 2行为 1 个 tuple
        - 对 .fq|.fastq 文件 4行为 1 个 tuple
        - 对 .txt[default] 文件 1行为 1 个 tuple
    3. 该方法识别 CELL_COSMO_TEST_ENV_READ_N_ROWS 环境变量，用于测试时读取指定行的数据

    :param filepath: the file path need to read
    :param ignore_test_env: 忽略`CELL_COSMO_TEST_ENV_READ_N_ROWS`变量
    :return: yield tuple data one by one
    """
    # 如果设置环境变量,则仅读指定行数量的内容,用来提高测试时运行的速度
    envkey = "CELL_COSMO_TEST_ENV_READ_N_ROWS"
    read_n_rows = os.environ.get(envkey, -1)
    if ignore_test_env:
        read_n_rows = -1
    if read_n_rows != -1:
        # 必须为正整数,且最小为1
        if not str(read_n_rows).isdigit():
            print(f"{envkey} in env must specify an integer,replace {read_n_rows} by default 100")
            read_n_rows = 100
        read_n_rows = int(read_n_rows)
        if read_n_rows < 1:
            print(f"{envkey} in env must >= 1,replace {read_n_rows} by default 100")
            read_n_rows = 100
    # 根据文件后缀自动识别open方式
    sub_suffix = filepath
    mode = "r"
    if filepath.endswith(".gz"):
        open_func = gzip.open
        mode = "rt"
        sub_suffix = sub_suffix[:-len(".gz")]
    else:
        open_func = open
    # 根据文件后缀自动识别zip的行数
    sub_suffix = sub_suffix.split(".")[-1]
    if sub_suffix in {"fa", "fasta"}:
        n = 2
    elif sub_suffix in {"fq", "fastq"}:
        n = 4
    else:
        n = 1
    with open_func(filepath, mode) as fh:
        for row in itertools.zip_longest(*[fh] * n):
            if read_n_rows != -1 and read_n_rows == 0:
                return
            elif read_n_rows != -1:
                read_n_rows -= 1
            if row[-1] is None:
                # 兼容文件异常换行的情况
                continue
            yield tuple(r.strip() for r in row)


