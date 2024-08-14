#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : starsolo_fake_run.py
@Time       : 2024/06/22
@Version    : 1.0
@Desc       : None
"""
import os
import shutil
import time
from cell_cosmo.commands import common
from cell_cosmo.rna.starsolo import Starsolo
from cell_cosmo.rna.starsolo import Mapping
from cell_cosmo.rna.starsolo import Cells
from cell_cosmo.rna.starsolo import Demultiplexing
import logging

fmt_str = '%(asctime)s,%(msecs)03d %(levelname)s %(name)s >> %(message)s'
logging.basicConfig(level=logging.INFO, format=fmt_str)


def starsolo_fake_run():
    target = "/mnt/data/万乘"
    outdir = "starsolo_fake_run/result"
    os.chdir(target)
    # if os.path.exists(outdir):
    #     # 删除 历史运行输出
    #     shutil.rmtree(outdir)
    # os.makedirs(outdir)

    kwargs = {
        "fq1": "/mnt/data/万乘/mb-1w_L04_R1.fq.gz",
        "fq2": "/mnt/data/万乘/mb-1w_L04_R2.fq.gz",
        "outdir": outdir,
        "sample": "MB",
        "chemistry_name": "DemoV1",
        "genomeDir": "/mnt/data/cellcosmo/20230706_lockdb/Mouse/STAR_CellCosmo",
        "picard_mem": 5,  # todo 删除
        "star_mem": 5,
        "thread": 6,
        "STAR_param": None,
        "consensus_fq": False,
        "outFilterMatchNmin": 50,
        "soloCellFilter": "EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000",
        "SAM_attributes": "NH HI nM AS CR UR CB UB GX GN ",
        "soloFeatures": "Gene GeneFull_Ex50pAS",

        "out_unmapped": False,
        # "samtools_index_param": None,  # todo 删除
        "allow_barcode_diff_num": 1,
        "subparser_assay": "rna",
    }
    start = time.time()
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
        # runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)

    end = time.time()
    print(end - start)


if __name__ == '__main__':
    starsolo_fake_run()
