#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : star.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.commands import common
from cell_cosmo.rna.star import Star


@click.command()
@click.option('--fq', required=True, help="R2 fastq file.")
@click.option('--consensus-fq', is_flag=True, default=False,
              help="A indicator that the input fastq has been consensused.")
@common.genomeDir
@click.option('--out-unmapped', is_flag=True, default=False, help="Output unmapped reads.")
@click.option('--out-filter-match-n-min', 'outFilterMatchNmin', default=0,
              help="Alignment will be output only if the number of matched bases"
                   "is higher than or equal to this value.")
@click.option('--out-filter-multimap-n-max', 'outFilterMultimapNmax', default=1,
              help="How many places are allowed to match a read at most.")
@click.option('--star-mem', default=30,
              help="Maximum memory that STAR can use.")
@click.option('--picard-mem', default=30,
              help="Maximum memory that picard can use.")
@click.option('--star-param', "STAR_param", default="",
              help='Additional parameters for the called software. Need to '
                   'be enclosed in quotation marks. For example: '
                   '`--{software}_param "--param1 value1 --param2 value2"`.')
@common.outdir
@common.sample
@common.thread4
# @common.debug
def star(**kwargs):
    kwargs.setdefault("subparser_assay", "rna")
    with Star(**kwargs) as runner:
        runner.run()
