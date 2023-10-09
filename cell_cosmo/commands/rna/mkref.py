#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : _MKRefBase.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.rna.MKRef4RNA import MKRef4RNA
from cell_cosmo.commands import common


@click.command()
@common.genome_name
@common.fasta
@click.option('--gtf', required=True, help=" Genome gtf file. Use absolute path or relative path to `genomeDir`.")
@click.option('--mt-gene-list', default=None, help="""Mitochondria gene list file. 
Use absolute path or relative path to `genomeDir`.
It is a plain text file with one gene per line.
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.""")
@click.option('--genomeSAindexNbases', "genomeSAindexNbases",  # genomeSAindexNbases 参数名会转换成小写,这里添加别名
              default=14, help="STAR param: genomeSAindexNbases.")
@click.option('--gene-name-as-name2', is_flag=True, default=False,
              help="control the param of `-geneNameAsName2` in gtfToGenePred")
@click.option('--star-param', "STAR_param", default="",
              help='Additional parameters for the called software. Need to '
                   'be enclosed in quotation marks. For example: '
                   '`--{software}_param "--param1 value1 --param2 value2"`.')
@common.thread6
# @common.dry_run
def mkref(**kwargs):
    kwargs.setdefault("subparser_assay", "rna")
    with MKRef4RNA(**kwargs) as runner:
        runner.run()
