#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : valid_polyt.py
@Time       : 2022/06/01
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.barcode.valid.base_valid import BaseValid, pp


class ValidPolyt(BaseValid):
    __filename__ = "no_polyt"

    def __init__(self, pattern, **kwargs):
        super(ValidPolyt, self).__init__(
            pattern=pattern,
            need_valid=kwargs.pop("use_polyt_valid_reads", False),
            need_output=kwargs.pop("output_no_polyt", True),  # TODO false mark
            n_allow=-1,
            **kwargs
        )
        if self.need_valid:
            assert pp.T in pattern.pattern_dict, f"check polyt but code {pp.T} not in pattern"
            self.rate = kwargs.pop("polyt_rate", 0.7)
            self.n_allow = int(pattern.get_pattern_len(pp.T) * self.rate)

    def _valid(self, r1, r2) -> bool:
        extract_reads = self.pattern.get_seq_str(r1[1], pp.T)
        passed = extract_reads.count("T") >= self.n_allow
        self._write(r1, r2, not passed)
        return passed
