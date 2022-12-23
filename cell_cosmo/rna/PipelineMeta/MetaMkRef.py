#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MetaMkRef.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta as C, CMDBase


class MetaMkRef(CMDBase):
    def __init__(self):
        super(MetaMkRef, self).__init__(
            section="mkref",
            meta=[
                C("--genome-name", required=True),
                C("--fasta", required=True),
                C("--gtf", required=True),
                C("--mt-gene-list"),
                C("--thread", default="6"),
                C("--genomeSAindexNbases", default="14"),
                C("--gene-name-as-name2", default="False", is_flag=True),
            ],
            cmd=f"CellCosmo rna mkref".split()
        )
