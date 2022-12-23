#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : common.py
@Time       : 2022/06/08
@Version    : 1.0
@Desc       : None
"""
import click

outdir = click.option('-o', "--outdir", required=True, help="Output directory.")
sample = click.option('-s', "--sample", required=True, help="Sample name.")
thread4 = click.option('-t', "--thread", default=4, help="Thread to use.")
debug = click.option('-d', "--debug", is_flag=True,
                     help="If this argument is used, may output additional file for debugging.")

thread6 = click.option('-t', '--thread', default=6, help="Threads to use.")
genome_name = click.option('--genome-name', required=True, default=None, help="genome name.")
dry_run = click.option('--dry-run', is_flag=True, default=False, help="Only write config file and exit.")
fasta = click.option('--fasta', required=True,
                     help="Genome fasta file. Use absolute path or relative path to `genomeDir`.")

# 必须添加别名,否则会变为小写
genomeDir = click.option('--genomeDir', "genomeDir", required=True, help="Genome directory.")
