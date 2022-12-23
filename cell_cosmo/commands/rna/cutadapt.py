#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : cutadapt_test.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.commands import common
from cell_cosmo.tools.cutadapt import Cutadapt


@click.command()
@click.option('--fq', required=True, help="Required. R2 reads from step Barcode.")
@click.option('--gzip', is_flag=True, default=False, help="Output gzipped fastq files.")
@click.option('--adapter-fasta', help="Additional adapter fasta file.")
@click.option('--minimum-length', default=20, help="Discard processed reads that are shorter than LENGTH.")
@click.option('--nextseq-trim', default=20, help="""Quality trimming of reads using two-color chemistry (NextSeq).
Some Illumina instruments use a two-color chemistry to encode the four bases.
This includes the NextSeq and the NovaSeq.
In those instruments, a `dark cycle` (with no detected color) encodes a G.
However, dark cycles also occur when sequencing `falls off` the end of the fragment.
The read then contains a run of high-quality, but incorrect `G` calls at its 3' end.
""")
@click.option("--overlap", default=10, help="""Since Cutadapt allows partial matches between the read and the adapter sequence,
short matches can occur by chance, leading to erroneously trimmed bases.
For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter.
To reduce the number of falsely trimmed bases, the alignment algorithm requires that
at least {overlap} bases match between adapter and read.
""")
@click.option("--insert", default=150, help="Read2 insert length.")
@click.option("--cutadapt-param", default="", help="Other cutadapt parameters. For example, --cutadapt_param '-g AAA' ")
@common.outdir
@common.sample
@common.thread4
# @common.debug
def cutadapt(**kwargs):
    kwargs.setdefault("subparser_assay", "rna")
    with Cutadapt(**kwargs) as runner:
        runner.run()
