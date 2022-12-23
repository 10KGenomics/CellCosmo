#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MetaSample.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta as C, CMDBase


class MetaSample(CMDBase):
    def __init__(self, fq1: str, outdir: str, sample: str):
        super(MetaSample, self).__init__(
            section="sample",
            meta=[
                C("--chemistry", default="auto"),
            ],
            cmd=f"CellCosmo rna sample "
                f"-o {outdir} "
                f"-s {sample}".split()
        )
