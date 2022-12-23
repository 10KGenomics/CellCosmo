#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : featureCounts.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.commands import common
from cell_cosmo.tools.FeatureCounts import FeatureCounts


@click.command()
@common.sample
@click.option("--input", "input_file", required=True, help="BAM file path.")
@common.outdir
@click.option("--gtf-type", default="exon", help="Specify feature type in GTF annotation")
@common.genomeDir
@common.thread4
# @common.debug
@click.option("--feature-counts-param", "featureCounts_param", default="",
              help="""Additional parameters for the 
called software. Need to be enclosed in quotation marks. For example: 
`--{software}_param "--param1 value1 --param2 value2"`.""")
def featureCounts(**kwargs):
    kwargs.setdefault("subparser_assay", "rna")
    with FeatureCounts(**kwargs) as runner:
        runner.run()
