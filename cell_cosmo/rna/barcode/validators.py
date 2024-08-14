#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : validators.py
@Time       : 2024/07/06
@Version    : 1.0
@Desc       : None
"""
import logging
from .models import (
    ValidParam, ValidPolytParam, ValidBarcodeParam,
    ValidLinkParam, ValidQualParam, CleanedReadsParam
)

logger = logging.getLogger(__name__)


class Validators:
    def __init__(self, n_polyt, batch_size, **kwargs):
        self.outdir = kwargs.pop("outdir")  # must in
        self.sample = kwargs.pop("sample")  # must in
        self.batch_size = batch_size
        self.is_gz = kwargs.pop("gzip", False)
        self.prefix = f"{self.outdir}/{self.sample}"
        self.suffix = ".gz" if self.is_gz else ""
        # self.temp_done_state_format = ".temp.%s.done_state"
        # self.temp_done_state_regexp = ".temp.*.done_state"
        self.temp_state_file = f"{self.outdir}/.run_chunk.state"

        with open(self.temp_state_file, "w") as fh:
            fh.write("0")

        self.params_polyt = ValidPolytParam(
            self.prefix, self.suffix,
            polyt_length=n_polyt, **kwargs
        )
        self.params_barcode = ValidBarcodeParam(
            self.prefix, self.suffix, **kwargs
        )
        self.params_link = ValidLinkParam(
            self.prefix, self.suffix, **kwargs
        )
        self.params_qual = ValidQualParam(
            self.prefix, self.suffix, **kwargs
        )
        self.params_cleaned = CleanedReadsParam(
            self.prefix, self.suffix
        )

        self.list = [
            self.params_cleaned,
            self.params_barcode,
            self.params_polyt,
            self.params_link,
            self.params_qual,
        ]
