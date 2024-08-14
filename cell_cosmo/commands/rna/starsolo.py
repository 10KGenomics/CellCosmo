#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : starsolo.py
@Time       : 2024/07/09
@Version    : 1.0
@Desc       : None
"""
import click
import logging
from cell_cosmo.commands import common
from cell_cosmo.rna.starsolo import Cells
from cell_cosmo.rna.starsolo import Mapping
from cell_cosmo.rna.starsolo import Starsolo
from cell_cosmo.rna.starsolo import Demultiplexing

logger = logging.getLogger(__name__)
SAM_attributes = 'NH HI nM AS CR UR CB UB GX GN '
SOLO_CELL_FILTER_ARGS = "0.99 10 45000 90000 500 0.01 20000 0.001 10000"


@click.command()
@click.option('--fq1', required=True, help="R1 fastq file. Multiple files are separated by comma.")
@click.option('--fq2', required=True, help="R2 fastq file. Multiple files are separated by comma.")
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
@click.option('--consensus-fq', is_flag=True, default=False,
              help="A indicator that the input fastq has been consensused.")
@common.genomeDir
@click.option('--star-mem', default=30,
              help="Maximum memory that STAR can use.")
@click.option('--star-param', "STAR_param", default="",
              help='Additional parameters for the called software. Need to '
                   'be enclosed in quotation marks. For example: '
                   '`--star_param "--param1 value1 --param2 value2"`.')
@click.option('--out-filter-match-n-min', 'outFilterMatchNmin', default=50,  # 0 or 50?
              help="Alignment will be output only if the number of matched bases"
                   "is higher than or equal to this value.")
@click.option('--sam_attributes', 'SAM_attributes',
              default=SAM_attributes,
              help=f"Additional attributes(other than {SAM_attributes}) to be added to SAM file")
@click.option('--solo-features', 'soloFeatures', type=str,
              default="Gene GeneFull_Ex50pAS",
              help=f"The same as the soloFeatures argument in STARsolo")
@click.option('--out_sam_type', 'outSAMtype',
              default="BAM SortedByCoordinate",
              help=f"type of SAM/BAM output")
@click.option('--solo_cell_filter_method', default="EmptyDrops_CR", help=f"Cellcalling Method")
@click.option('--solo_cell_filter_n_expect', default=3000, help=f"Expect_num")
@click.option('--solo_cell_filter_args', default=SOLO_CELL_FILTER_ARGS, help='solo cell filter args')
@click.option('--adapter_3p',
              default="AAAAAAAAAAAA",
              help=f"Adapter sequence to clip from 3 prime. Multiple sequences are seperated by space")
@click.option('--out-unmapped', is_flag=True, default=False, help="Output unmapped reads.")
@common.outdir
@common.sample
@common.thread4
# @common.debug
def starsolo(**kwargs):
    kwargs.setdefault("subparser_assay", "rna")

    sa = kwargs.get('SAM_attributes')
    out_sam_type = kwargs.get('outSAMtype', None)
    if str(out_sam_type).strip().lower() in {'none', 'null', 'false', ''}:
        if sa != SAM_attributes:  # 用户修改了默认值，此时该参数不生效
            logger.warning(f"--out_sam_type设置为空时，将忽略 --sam_attributes 参数")
        kwargs['outSAMtype'] = None
        kwargs['SAM_attributes'] = None
    else:
        kwargs['SAM_attributes'] = sa

    solo_cell_filter_method = kwargs.pop('solo_cell_filter_method')
    solo_cell_filter_n_expect = kwargs.pop('solo_cell_filter_n_expect')
    solo_cell_filter_args = kwargs.pop('solo_cell_filter_args')
    kwargs['soloCellFilter'] = (f'{solo_cell_filter_method} '
                                f'{solo_cell_filter_n_expect} '
                                f'{solo_cell_filter_args}')

    # for k, v in kwargs.items():
    #     print("==>", k, v)
    with Starsolo(**kwargs) as runner:
        runner.run()
        q30_cb, q30_umi = runner.get_Q30_cb_UMI()
        # 0.9211 0.89505

    with Mapping(**kwargs) as runner:
        runner.run()
        valid_reads, corrected = runner.get_vc()
    #
    with Cells(valid_reads, **kwargs) as runner:
        runner.run()
        n_reads, q30_RNA = runner.n_reads, runner.q30_RNA
    with Demultiplexing(valid_reads, n_reads, corrected,
                        q30_cb, q30_umi, q30_RNA, **kwargs) as runner:
        runner.run()
