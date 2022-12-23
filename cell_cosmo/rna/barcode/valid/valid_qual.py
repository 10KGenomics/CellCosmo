#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : valid_qual.py
@Time       : 2022/06/01
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.barcode.valid.base_valid import BaseValid, pp


class ValidQual(BaseValid):
    __filename__ = "low_qual"

    def __init__(self, pattern, offset=33, **kwargs):
        self.offset = offset
        self.low_qual = kwargs.pop("low_qual", 0)
        self.low_num = kwargs.pop("low_num", 2)
        super(ValidQual, self).__init__(
            pattern=pattern,
            need_valid=self.low_qual != 0,
            need_output=kwargs.pop("output_low_qual", True),  # TODO false mark
            n_allow=self.low_num,
            **kwargs
        )
        if self.need_valid:
            assert all([c in pattern.pattern_dict for c in {pp.T, pp.U}]), \
                f"check qual but code {pp.T} or {pp.U} not in pattern"

    def _valid(self, r1, r2) -> bool:
        qual_str = self.pattern.get_seq_str(r1[3], [pp.C, pp.U])
        low_qual_num = len([
            q for q in qual_str
            if self.qual_int(q, self.offset) < self.low_qual])
        passed = low_qual_num <= self.n_allow
        self._write(r1, r2, not passed)
        return passed

    @staticmethod
    def ord2chr(q, offset=33):
        return chr(int(q) + offset)

    @staticmethod
    def qual_int(char, offset=33):
        return ord(char) - offset
