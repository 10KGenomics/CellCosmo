#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : base_valid.py
@Time       : 2022/06/01
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.barcode import parse_pattern as pp
from cell_cosmo.rna.barcode.valid.base_output import BaseOutput


class BaseValid(BaseOutput):

    def __init__(self, need_valid: bool, pattern: pp.PatternParser,
                 n_allow: int, need_output: bool = False, **kwargs):
        super(BaseValid, self).__init__(need_valid, need_output, **kwargs)
        self.pattern = pattern
        self.n_allow = n_allow

    def _valid(self, r1, r2) -> bool:
        raise NotImplementedError

    def valid(self, r1, r2) -> bool:
        if not self.need_valid:
            return True
        return self._valid(r1, r2)
