#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MetaFeatureCounts.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta as C, CMDBase


class MetaFeatureCounts(CMDBase):
    def __init__(self, inp, genomeDir, outdir, sample):
        super(MetaFeatureCounts, self).__init__(
            section="featureCounts",
            meta=[
                C("--gtf-type", default="exon"),  # TODO range？
                C("--feature-counts-param", default=""),  # -ignoreGroupsWithoutExons
                C("--thread", default="4"),
            ],
            cmd=f"CellCosmo rna featureCounts "
                f"--input {inp} "
                f"--genomeDir {genomeDir} "
                f"-o {outdir} "
                f"-s {sample}".split()
        )
