#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MetaStar.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta as C, CMDBase


class MetaStar(CMDBase):
    def __init__(self, fq, genomeDir, outdir, sample):
        super(MetaStar, self).__init__(
            section="star",
            meta=[
                C("--consensus-fq", default="False", is_flag=True),
                C("--out-unmapped", default="False", is_flag=True),
                C("--out-filter-match-n-min", default="0"),
                C("--out-filter-multimap-n-max", default="1"),
                C("--picard-mem", default="30"),
                C("--star-mem", default="30"),
                C("--star-param"),
                C("--thread", default="4"),
            ],
            cmd=f"CellCosmo rna star "
                f"--fq {fq} "
                f"--genomeDir {genomeDir} "
                f"-o {outdir} "
                f"-s {sample}".split()
        )
