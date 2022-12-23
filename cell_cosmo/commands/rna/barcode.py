#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : barcode.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.rna.barcode.runner import BarcodeSplitterRunner
from cell_cosmo.commands import common


@click.command()
@click.option('--fq1', required=True,
              help="Fastq file containing library information, "
                   "the file suffix must be .fastq.gz or .fq.gz")
@click.option('--fq2', required=True,
              help="Another fastq file that paired with fq1")
@common.sample
@click.option('-n', '--chemistry-name',
              help="Chemistry name,If this parameter is specified, "
                   "the program will load the chemistry information corresponding "
                   "to the name from the preset database")
@click.option('-c', '--chemistry-config',
              help="Chemistry configuration file path, "
                   "the same level directory contains chemistry information "
                   "such as barcode library and link string")
@click.option('-p', '--pattern', help="""\b
The pattern of sequences,eg:`C5U3L15U3C6U3L6C5T30`,
The number after the letter represents the number of bases
- `C`: barcode
- `L`: linker
- `U`: UMI
- `T`: poly T""")
@click.option('--use-link-valid-reads', is_flag=True, default=False,
              help="Validate reads using link sequences")
@click.option('--use-barcode-valid-reads', is_flag=True, default=False,
              help="Validate reads using barcode sequences")
@click.option('--use-polyt-valid-reads', is_flag=True, default=False,
              help='Validate reads using polyt sequences')
@click.option('--allow-link-diff-num', type=click.INT, default=2,
              help="mismatch number with link sequence,"
                   "this parameter only takes effect when `--use-link-valid-reads` is specified")
@click.option('--allow-barcode-diff-num', type=click.INT, default=1,
              help='mismatch number with barcode library,'
                   'this parameter is invalid when `--use-barcode-valid-reads` is specified')
@click.option('--low-qual', type=click.INT, default=0,
              help="The barcode and UMI whose phred value are lower than "
                   "--low-qual will be regarded as low-quality bases.")
@click.option('--low-num', type=click.INT, default=2,
              help="The maximum number of low qual bases allowed in barcode and UMI.")
@click.option('--polyt-rate', type=click.FLOAT, default=0.7,
              help="The proportion of T bases in polyT that need to be satisfied,range is [0,1]."
                   "When specified as 0, the limit of polyT is ignored;"
                   "When specified as 1, it means all T bases in polyT;"
                   "When the specified pattern contains T configuration, "
                   "the default value of this parameter is 2/3, otherwise the default value is 0")
@click.option('--gzip', is_flag=True, default=False,
              help="Output fastq file in compressed format")
@click.option('--output-r1', is_flag=True, default=False,
              help="Output the R1 sequence corresponding to the valid sequence")
@common.outdir
# @common.debug
@common.thread4
def barcode(**kwargs):
    """
    split barcode

    """
    kwargs.setdefault("subparser_assay", "rna")
    # TODO 有参和无参对于参数的校验提示,
    # 1. 对于无参，(--use-xx-valid-reads)不能为true,-c参数不能存在,(仍然可以使用polyt和qual验证reads?)
    # 2. 对于有参, -c参数必须存在
    with BarcodeSplitterRunner(**kwargs) as runner:
        runner.run()
