#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MetaCount.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta as C, CMDBase


class MetaCount(CMDBase):
    def __init__(self, bam, genomeDir, outdir, sample):
        super(MetaCount, self).__init__(
            section="count",
            meta=[
                C("--expected-cell-num", "3000"),
                C("--cell-calling-method", "EmptyDrops_CR", choice=['auto', 'EmptyDrops_CR']),
                C("--n-umi-filter", "0"),
                C("--barcode-correct-limit", "0.01"),
                C("--umi-correct-limit", "0.1"),
                C("--force-cell-num"),
                C("--thread", default="4")
            ],
            cmd=(f"CellCosmo rna count "
                 f"-b {bam} "
                 f"--genomeDir {genomeDir} "
                 f"-o {outdir} "
                 f"-s {sample}".split())
        )
