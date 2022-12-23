#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MetaAnalysis.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta as C, CMDBase


class MetaAnalysis(CMDBase):
    def __init__(self, matrix_file, genomeDir, outdir, sample):
        super(MetaAnalysis, self).__init__(
            section="analysis",
            meta=[
                C("--thread", default="4")
            ],
            cmd=(f"CellCosmo rna analysis  "
                 f"--matrix-file {matrix_file} "
                 f"--genomeDir {genomeDir} "
                 f"-o {outdir} "
                 f"-s {sample}".split())
        )
