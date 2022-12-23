#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : Analysis.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.commands import common
from cell_cosmo.rna.Analysis import Analysis


@click.command()
@common.genomeDir
@click.option('--matrix-file', required=True,
              help="Matrix_10X directory from step count.")
@common.sample
@common.outdir
@common.thread4
# @common.debug
def analysis(**kwargs):
    kwargs.setdefault("subparser_assay", "rna")
    with Analysis(**kwargs) as runner:
        runner.run()
