#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MetaBarcode.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta as C, CMDBase


class MetaBarcode(CMDBase):
    def __init__(self, fq1: str, fq2: str, outdir: str, sample: str):
        super(MetaBarcode, self).__init__(
            section="barcode",
            meta=[
                C("--chemistry-name"),
                C("--chemistry-config"),
                C("--pattern", required=True),
                C("--use-link-valid-reads", default="False", is_flag=True),
                C("--use-barcode-valid-reads", default="False", is_flag=True),
                C("--use-polyt-valid-reads", default="False", is_flag=True),
                C("--allow-link-diff-num", default="2"),
                C("--allow-barcode-diff-num", default="1"),
                C("--low-qual", default="0"),
                C("--low-num", default="2"),
                C("--polyt-rate", default="0.7"),
                C("--output-r1", default="False", is_flag=True),
                C("--gzip", default="False", is_flag=True),
                C("--thread", default="4"),
            ],
            cmd=f"CellCosmo rna barcode --fq1 {fq1} --fq2 {fq2} -o {outdir} -s {sample}".split()
        )
