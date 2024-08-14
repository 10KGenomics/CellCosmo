#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : yield_batch_data.py
@Time       : 2024/07/05
@Version    : 1.0
@Desc       : None
"""
import dnaio
import logging
from xopen import xopen
from itertools import zip_longest

logger = logging.getLogger(__name__)


def yield_batch_data(fq1: str, fq2: str, batch_size=1024 * 1024 * 100):
    fq1_list = fq1.split(",")
    fq2_list = fq2.split(",")
    fq_num = len(fq1_list)
    if fq_num != len(fq2_list):
        raise Exception('fq1 and fq2 must be the same file number!')

    start_index = 0

    for i in range(fq_num):
        with (xopen(fq1_list[i], "rb") as fh1,
              xopen(fq2_list[i], "rb") as fh2):
            for chunk_id, (chunk1, chunk2) in enumerate(
                    dnaio.read_paired_chunks(fh1, fh2, batch_size)
            ):
                chunk1 = str(chunk1.tobytes(), encoding='utf-8').strip().split("\n")
                chunk2 = str(chunk2.tobytes(), encoding='utf-8').strip().split("\n")
                batch_data = list(zip_longest(*[iter(chunk1)] * 4, *[iter(chunk2)] * 4))
                yield chunk_id, start_index, batch_data
                start_index += len(batch_data)
                # if start_index > 200 * 10000:
                #     break
