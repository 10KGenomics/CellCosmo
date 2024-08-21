#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : barcode.py
@Time       : 2024/06/22
@Version    : 1.0
@Desc       : None
"""
import os
import shutil
import time
from cell_cosmo.rna.barcode.runner import BarcodeSplitterRunner
import logging

fmt_str = '%(asctime)s,%(msecs)03d %(levelname)s %(name)s >> %(message)s'
logging.basicConfig(level=logging.INFO, format=fmt_str)


def test_1():
    target = "/mnt/data/万乘"
    outdir = "T1w/result"
    os.chdir(target)
    if os.path.exists(outdir):
        # 删除 历史运行输出
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    kwargs = {
        "fq1": "/mnt/data/万乘/mb-1w_L04_R1.fq.gz",
        "fq2": "/mnt/data/万乘/mb-1w_L04_R2.fq.gz",
        "outdir": outdir,
        "sample": "MB",
        "chemistry_name": "DemoV1",
        "use_link_valid_reads": True,
        "use_barcode_valid_reads": True,
        "use_polyt_valid_reads": True,
        "gzip": True,
        "thread": 32,
        "chemistry_config": None,
        "allow_link_diff_num": 2,
        "allow_barcode_diff_num": 1,
        "low_qual": 0,
        "low_num": 2,
        "polyt_rate": 0.7,
        "output_r1": False,
        "pattern": None,
        "subparser_assay": "rna",
    }
    start = time.time()
    print("start to run it ..")
    with BarcodeSplitterRunner(**kwargs) as runner:
        runner.run()
    end = time.time()
    print(end - start)


if __name__ == '__main__':
    test_1()