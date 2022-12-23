#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : count.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.commands import common
from cell_cosmo.tools.Count import Count


# # barcode_diff_limit' and 'umi_diff_limit
@click.command()
@common.genomeDir
@click.option('-b', '--bam', required=True, help="BAM file from feature_counts.")
@click.option('--expected-cell-num', default=3000, help="Expected cell number.")
# choices=['auto', 'EmptyDrops_CR'],
@click.option('--cell-calling-method', default="EmptyDrops_CR", help="Choose from [`auto`, `EmptyDrops_CR`]")
# @click.option('--barcode-diff', default=1,
#               help="diff base low than this value,the barcode will be merged."
#                    "if set 0,while not correct barcode")
# @click.option('--barcode-correct-limit', default=0.1,
#               help="if correct barcode,low count/high count need less than this value."
#                    "if set 1,merge low to high for all match case.")
# @click.option('--umi-diff', default=1,
#               help="diff base low than this value,the umi will be merged."
#                    "if set 0,while not correct umi")
# @click.option('--umi-correct-limit', default=0.1,
#               help="if correct umi,low count/high count need less than this value."
#                    "if set 1,merge low to high for all match case.")
@click.option('--n-umi-filter', default=0,
              help="when correct barcode and umi, the umis count less than this value will "
                   "be discard.set this params will accelerated running speed but "
                   "some data will be discarded.")
@click.option('--barcode-correct-limit', default=0.01,
              help="when barcode correct,low_umis_count/high_umis_count need less than this value."
                   "if set 1,merge low to high for all match case.")
@click.option('--umi-correct-limit', default=0.1,
              help="when umi correct,low count/high count need less than this value."
                   "if set 1,merge low to high for all match case.")
@click.option('--force-cell-num', default=None, help="Force the cell number to be this number.")
@common.sample
@common.outdir
@common.thread4
# @common.debug
def count(**kwargs):
    # n_umi_filter=20, filter_limit=0.01, percent=0.1,

    kwargs.setdefault("subparser_assay", "rna")
    with Count(**kwargs) as runner:
        runner.run()
