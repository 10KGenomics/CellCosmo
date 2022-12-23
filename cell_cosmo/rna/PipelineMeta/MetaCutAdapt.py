#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MetaCutAdapt.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta as C, CMDBase


class MetaCutAdapt(CMDBase):
    def __init__(self, fq, outdir, sample):
        super(MetaCutAdapt, self).__init__(
            section="cutadapt",
            meta=[
                C("--adapter-fasta"),
                C("--minimum-length", "20"),
                C("--nextseq-trim", "20"),
                C("--overlap", "10"),
                C("--insert", "150"),
                C("--cutadapt-param", ""),
                C("--gzip", default="False", is_flag=True),
                C("--thread", default="4")
            ],
            cmd=f"CellCosmo rna cutadapt "
                f"--fq {fq} "
                f"-o {outdir} "
                f"-s {sample}".split()
        )
