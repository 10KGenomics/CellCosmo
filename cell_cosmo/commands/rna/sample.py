#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : sample.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.commands import common
from cell_cosmo.tools.Sample import Sample


@click.command()
@click.option('-f', "--fq1", help="read1 fq file")
# TODO choices for chemistry
@click.option('-c', "--chemistry", default='auto', help="chemistry version,#TEST,current use pattern instead")
@common.outdir
@common.sample
# @common.thread4
# @common.debug
def sample(**kwargs):
    kwargs.setdefault("subparser_assay", "rna")
    with Sample(**kwargs) as runner:
        runner.run()
