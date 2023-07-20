#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : BaseOut.py
@Time       : 2023/07/20
@Version    : 1.0
@Desc       : None
"""
import os
import pandas as pd


class BaseOut:
    def __init__(self, outdir, sample):
        self.outdir = outdir
        self.sample = sample

    def to_csv(self, filename: str, df: pd.DataFrame, sep="\t"):
        if sep == "\t":
            if not filename.endswith(".tsv"):
                filename = f"{filename}.tsv"
        elif sep == ",":
            if not filename.endswith(".csv"):
                filename = f"{filename}.csv"
        else:
            raise Exception(f"不支持该文件分割符输出：`{sep}`")

        compression_type = os.getenv("CELLCOSMO_COMPRESSION_STRATEGY", 1)
        if str(compression_type) == "1":
            # 输出压缩格式
            df.to_csv(f"{self.outdir}/{self.sample}_{filename}.gz", sep=sep, compression="gzip", index=False)
        else:
            df.to_csv(f"{self.outdir}/{self.sample}_{filename}", sep=sep, index=False)
