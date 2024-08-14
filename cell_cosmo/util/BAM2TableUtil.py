#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : BAM2TableUtil.py
@Time       : 2022/06/13
@Version    : 1.0
@Desc       : None
"""
import os
import random
import string
import time
import pysam
import pandas as pd
from itertools import groupby
from contextlib import contextmanager
from cell_cosmo.util.BarcodeCorrectUtil import ConstNS
from cell_cosmo.util import runtime
from tqdm import tqdm


def key_func(raw):
    # 从bam文件中提取的序列名称
    # 期望的字符串模式 {barcode}_{umi}_{raw_line_no}
    # 使用barcode作为key
    return raw.query_name.split('_', 1)[0]


@contextmanager
def get_temp_file():
    while True:
        letter = string.ascii_letters + string.digits
        temp = "".join(random.choices(letter, k=8))
        temp_file = f".{temp}"
        if not os.path.exists(temp_file):
            break
    try:
        yield temp_file
    except Exception as e:
        print(str(e))
    if os.path.exists(temp_file):
        os.remove(temp_file)


@runtime(f"{__name__}.bam2count_table")
def bam2count_table(bam, out_prefix) -> pd.DataFrame:
    count_detail_raw_file = f'{out_prefix}_count_detail.raw.tsv'

    # bam 文件经过name排序,这里迭代出来barcode相同的一批数据
    before_barcode = None

    save = pysam.set_verbosity(0)
    sam_file = pysam.AlignmentFile(bam, "rb")
    pysam.set_verbosity(save)
    data = []
    # bam 文件经过name排序,这里迭代出来barcode相同的一批数据
    for barcode, g in tqdm(groupby(iter(sam_file), key_func), desc="load bam:"):
        if barcode == before_barcode:
            raise Exception(f"the bam file {bam} must sort by name,"
                            f"and seq name start with `@barcode_umi...`")
        before_barcode = barcode
        for seg in g:
            (_, umi) = seg.query_name.split('_')[:2]
            if not seg.has_tag('XT'):
                continue
            # 只统计包含gene id的数据
            gene_id = seg.get_tag('XT')
            data.append([barcode, gene_id, umi])

    cols = ConstNS.get_0_cols()
    df = pd.DataFrame(data=data, columns=cols)
    # count 和 size 在这里结果相同,因为没有NA值
    df = df.groupby(cols)[ConstNS.umi].count().reset_index(name=ConstNS.count)
    # todo async output
    # df.to_csv(count_detail_raw_file, sep="\t", index=False)
    _out(count_detail_raw_file, df)

    assert df.shape[0] != 0, f"No data left!!"
    return df


def _out(output, df: pd.DataFrame):
    compression_type = os.getenv("CELLCOSMO_COMPRESSION_STRATEGY", 1)
    if str(compression_type) == "1" and not output.endswith(".gz"):
        df.to_csv(f"{output}.gz", sep="\t", compression="gzip", index=False)
    else:
        df.to_csv(f"{output}.gz", sep="\t", index=False)
