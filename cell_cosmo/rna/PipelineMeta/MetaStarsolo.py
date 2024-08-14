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


class MetaStarsolo(CMDBase):
    def __init__(self, fq1, fq2, genomeDir, outdir, sample, **kwargs):
        super(MetaStarsolo, self).__init__(
            section="starsolo",
            meta=[
                C("--chemistry-name"),
                C("--chemistry-config"),
                C("--pattern"),
                # C("--consensus-fq", default="False", is_flag=True),
                C("--out-unmapped", default="False", is_flag=True),
                C("--out-filter-match-n-min", default="50"),
                C("--out-filter-multimap-n-max", default="1"),
                C("--star-mem", default="30"),
                C("--star-param", default=""),
                C("--solo_cell_filter_method", default="EmptyDrops_CR"),
                C("--solo_cell_filter_n_expect", default="3000"),
                C("--solo_cell_filter_args", default="0.99 10 45000 90000 500 0.01 20000 0.001 10000"),
                C("--out_sam_type", default="BAM SortedByCoordinate"),
                C("--solo-features", default="Gene GeneFull_Ex50pAS"),
                C("--sam_attributes", default="NH HI nM AS CR UR CB UB GX GN"),
                C("--thread", default="4"),
            ],
            cmd=f"CellCosmo rna starsolo "
                f"--fq1 {fq1} "
                f"--fq2 {fq2} "
                f"--genomeDir {genomeDir} "
                f"-o {outdir} "
                f"-s {sample}".split()
        )
