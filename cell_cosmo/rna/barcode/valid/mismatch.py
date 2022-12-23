#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : mismatch.py
@Time       : 2022/06/01
@Version    : 1.0
@Desc       : None
"""
from itertools import combinations, product


def yield_all_mismatch_seq(seq, n_mismatch=1, bases='ACGTN'):
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for pos in combinations(range(seq_len), n_mismatch):
        seq_pos = [[base] for base in seq]
        for loc in pos:
            seq_pos[loc] = list(bases)
        for poss in product(*seq_pos):
            yield ''.join(poss)


def yield_concat_seq(seqs_list):
    # 将多个barcode库按顺序进行组合
    for concat_seq in product(*seqs_list):
        concat_seq = ''.join(concat_seq)
        yield concat_seq
