#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : mixin_mismatch_dict.py
@Time       : 2022/06/08
@Version    : 1.0
@Desc       : None
"""
import abc
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


class MixinMismatchDict(mataclass=abc.ABCMeta):
    def _get_mismatch_dict(self, file,n_allow):
        a_set = set()
        a_dct = {}
        self.s, dir_path = params
        if self._files == "" and self.is_checked:
            raise Exception(f"use barcode valid reads,but not config barcode files")

        if self.is_checked:
            # 这种情形下初始化barcode文件，触发对barcode文库信息的校验
            # self.xxbarcode_files = []
            filepaths = [self.join_path(dir_path, f) for f in self._files.split(",")]
            pattern_list = self.pattern.get_pattern("C")
            assert len(filepaths) == len(pattern_list), \
                f"pattern string has {len(pattern_list)} C," \
                f"but config {len(filepaths)} barcode files"
            seqs = []
            for f, (s, e) in zip(filepaths, pattern_list):
                seqs.append(self.read_list(f, e - s))
            for seq in yield_concat_seq(seqs):
                a_set.add(seq)
                for cb in yield_all_mismatch_seq(seq, n_allow):
                    a_dct[cb] = seq

    @abc.abstractmethod
    def get_mismatch_dict(self):
        pass
