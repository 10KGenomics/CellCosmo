#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : valid_barcode.py
@Time       : 2022/06/01
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.barcode.valid.base_valid import BaseValid, pp
from cell_cosmo.rna.barcode.valid.mixin_cfg_files import CfgFilesMixin
from cell_cosmo.rna.barcode.valid.mismatch import yield_all_mismatch_seq, yield_concat_seq
from cell_cosmo.util.runtime import runtime


class ValidBarcode(BaseValid, CfgFilesMixin):
    __filename__ = "no_barcode"

    def __init__(self, pattern, **kwargs):
        super(ValidBarcode, self).__init__(
            pattern=pattern,
            need_valid=kwargs.pop("use_barcode_valid_reads", False),
            need_output=kwargs.pop("output_error_barcode", True),  # TODO false mark
            n_allow=kwargs.pop("allow_barcode_diff_num", 1),
            **kwargs
        )
        self._files = None
        self.barcode_set = set()
        self.barcode_mismatch_dit = {}
        if self.need_valid:
            assert pp.C in pattern.pattern_dict, f"check barcode but code {pp.C} not in pattern"

    def _valid(self, r1, r2) -> bool:
        extract_reads = self.pattern.get_seq_str(r1[1], pp.C)
        passed = extract_reads in self.barcode_set or extract_reads in self.barcode_mismatch_dit
        self._write(r1, r2, not passed)
        return passed

    @property
    def files(self):
        return self._files

    @files.setter
    @runtime('cell_cosmo.rna.barcode.get_barcode_mismatch_dict')
    def files(self, params):
        self._files, dir_path = params
        if self._files == "" and self.need_valid:
            raise Exception(f"use barcode valid reads,but not config barcode files")

        if self.need_valid:
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
                self.barcode_set.add(seq)
                for cb in yield_all_mismatch_seq(seq, self.n_allow):
                    self.barcode_mismatch_dit[cb] = seq
