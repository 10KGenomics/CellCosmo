#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : stat.py
@Time       : 2024/06/27
@Version    : 1.0
@Desc       : None
"""
import os
import dnaio


def _stat_file_reads(f1, f2):
    s1 = os.path.basename(f1)
    s2 = os.path.basename(f2)
    with dnaio.open(f1, f2, mode="r") as fh:
        num = 0
        for _ in fh:
            num += 1
        print(f"{s1}/{s2} count is: {num}")


def stat(dirpath):
    print(f"===========>>>>> stat filepath: {dirpath}")
    filenames = set()
    for f in os.listdir(dirpath):
        if not f.endswith(".gz"):
            continue
        if f.endswith("1.fq.gz"):
            filenames.add(f[:-len("1.fq.gz")])
    filenames = list(filenames)
    filenames.sort()
    for fn in filenames:
        f1 = f"{dirpath}/{fn}1.fq.gz"
        f2 = f"{dirpath}/{fn}2.fq.gz"
        _stat_file_reads(f1, f2)


def main():
    stat("/mnt/data/万乘/T1w/result")
    stat("/mnt/data/万乘/O1w/result")


if __name__ == '__main__':
    main()
