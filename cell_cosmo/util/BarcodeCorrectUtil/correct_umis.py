#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : correct_umis.py
@Time       : 2022/08/04
@Version    : 1.0
@Desc       : None
"""
import os
import logging
import itertools
import pandas as pd
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
import time
from cell_cosmo.util.BarcodeCorrectUtil.db_util import DBUtil
from cell_cosmo.util import distance
from collections import defaultdict
from threading import Thread

logger = logging.getLogger(__name__)


def multiprocess_callback(pbar: tqdm, res):
    pbar.update()


def correct_umi(umi_dict, percent=0.1):
    """
    Correct umi_dict in place.
    Args:
        umi_dict: {umi_seq: umi_count}
        percent: if hamming_distance(low_seq, high_seq) == 1 and
            low_count / high_count < percent, merge low to high.
    Returns:
        n_corrected_umi: int
        n_corrected_read: int
    """
    n_corrected_umi = 0
    n_corrected_read = 0

    # sort by value(UMI count) first, then key(UMI sequence)
    umi_arr = sorted(
        umi_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    while True:
        # break when only highest in umi_arr
        if len(umi_arr) == 1:
            break
        umi_low = umi_arr.pop()
        low_seq = umi_low[0]
        low_count = umi_low[1]

        for umi_kv in umi_arr:
            high_seq = umi_kv[0]
            high_count = umi_kv[1]
            if float(low_count / high_count) > percent:
                break
            if distance.hamming_distance(low_seq, high_seq) == 1:
                n_low = umi_dict[low_seq]
                n_corrected_umi += 1
                n_corrected_read += n_low
                # merge
                umi_dict[high_seq] += n_low
                del (umi_dict[low_seq])
                break
    return n_corrected_umi, n_corrected_read


class CorrectUMIs:
    def __init__(self, db_util: DBUtil, df: pd.DataFrame, output: str, percent=0.1, thread=4):
        """df: count_detail type data"""
        self.tasks = df.groupby(DBUtil.barcode)
        self.db_util = db_util
        self.percent = percent
        self.thread = thread
        self.output = output
        self.diff_num = 1

    def out(self, df: pd.DataFrame):
        compression_type = os.getenv("CELLCOSMO_COMPRESSION_STRATEGY", 1)
        output = self.output
        if str(compression_type) == "1" and not output.endswith(".gz"):
            output = f"{self.output}.gz"
        df.to_csv(output, sep="\t", compression="gzip", index=False)

    def correct(self):
        data = []
        for barcode, df in tqdm(self.tasks, desc='step3: correct umi in barcode'):
            gene_umi_dict = defaultdict(lambda: defaultdict(int))
            for ind, seg in df.iterrows():
                gene_id = seg.gene_id
                umi = seg.umi
                count = seg['count']  # seg.count 是method，不是count列的值
                gene_umi_dict[gene_id][umi] += count
            for gene_id in gene_umi_dict:
                correct_umi(gene_umi_dict[gene_id], percent=self.percent)

            for gene_id in gene_umi_dict:
                for umi in gene_umi_dict[gene_id]:
                    count = gene_umi_dict[gene_id][umi]
                    data.append([barcode, gene_id, umi, count])
        df = pd.DataFrame(data=data, columns=DBUtil.COUNT_DETAIL_COLS)
        # 使用异步线程避免输出文件阻塞
        thr = Thread(target=self.out, args=(df,))
        thr.start()
        return df
