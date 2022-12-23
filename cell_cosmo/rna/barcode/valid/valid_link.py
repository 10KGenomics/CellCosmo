#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : valid_link.py
@Time       : 2022/06/01
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.barcode.valid.base_valid import BaseValid, pp
from cell_cosmo.rna.barcode.valid.mixin_cfg_files import CfgFilesMixin
from cell_cosmo.rna.barcode.valid.mismatch import yield_all_mismatch_seq, yield_concat_seq
from cell_cosmo.util.runtime import runtime


class ValidLink(BaseValid, CfgFilesMixin):
    __filename__ = "no_linker"

    def __init__(self, pattern, **kwargs):
        super(ValidLink, self).__init__(
            pattern=pattern,
            need_valid=kwargs.pop("use_link_valid_reads", False),
            need_output=kwargs.pop("output_no_linker", True),  # TODO false mark
            n_allow=kwargs.pop("link_mismatch_num_allow", 2),
            **kwargs
        )
        if self.need_valid:
            assert pp.L in pattern.pattern_dict, f"check linker but code {pp.L} not in pattern"
        self._files = None
        self.link_set = set()
        self.link_mismatch_dit = {}

    def _valid(self, r1, r2) -> bool:
        extract_reads = self.pattern.get_seq_str(r1[1], pp.L)
        passed = extract_reads in self.link_set or extract_reads in self.link_mismatch_dit
        self._write(r1, r2, not passed)
        return passed

    @property
    def files(self):
        return self._files

    @files.setter
    @runtime('cell_cosmo.rna.barcode.get_link_mismatch_dict')
    def files(self, params):
        self._files, dir_path = params
        if self._files == "" and self.need_valid:
            raise Exception(f"use link valid reads,but not config link files")

        if self.need_valid:
            filepaths = [self.join_path(dir_path, f) for f in self._files.split(",")]
            pattern_list = self.pattern.get_pattern("L")
            assert len(filepaths) == len(pattern_list), \
                f"pattern string has {len(pattern_list)} L," \
                f"but config {len(filepaths)} link files"
            seqs = []
            for f, (s, e) in zip(filepaths, pattern_list):
                seqs.append(self.read_list(f, e - s))
            for seq in yield_concat_seq(seqs):
                self.link_set.add(seq)
                for cb in yield_all_mismatch_seq(seq, self.n_allow):
                    self.link_mismatch_dit[cb] = seq
