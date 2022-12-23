#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : BarcodeStatAfterCorrect.py
@Time       : 2022/07/08
@Version    : 1.0
@Desc       : None
"""
import os
import json
import pandas as pd


def get_stat_data(input_file):
    # the file columns is pos,sorted_group
    _cols = ["intended_barcode", "neighbor_barcode", "intended_size",
             "neighbor_size", "position", "intended_base", "neighbor_base"]
    df = pd.read_table(input_file, header=0)
    assert df.columns.tolist() == _cols, \
        f"barcode correct raw file first line cols is `{'`,`'.join(_cols)}`"
    pos, sorted_groups = _cols
    print(df.head(100))
    base_sort = "ACGT"
    xAxis_data = [f"{i}{b}" for i in range(1, 17) for b in base_sort]
    print(xAxis_data)
    data = {k: [0 for _ in range(16 * 4)] for k in base_sort}
    for current_pos, ddf in df.groupby(pos):
        stat = {}
        for bs in ddf[sorted_groups].tolist():
            bc_all = [bc for bc in bs.split(";")]
            k, v = bc_all[0].split(":")
            intended = k[current_pos - 1]
            i_num = int(v)
            if intended not in stat:
                stat[intended] = {}
            for bc in bc_all[1:]:
                k, v = bc.split(":")
                base = k[current_pos - 1]
                num = int(v)
                stat[intended].setdefault(base, 0)
                stat[intended][base] += num
        print(current_pos, stat)
        for base in base_sort:
            for i, bb in enumerate(base_sort):
                data[base][(current_pos - 1) * 4 + i] = stat[base].get(bb, 0)
    #         # 1G
    #         for i, bb in enumerate(base_sort):
    #             data[key][i] = stat[base].get(bb, 0)
    # print(json.dumps(xAxis_data))
    for na, dd in data.items():
        # print(na, dd)
        # , "label": {"show": "true"}
        x = {"name": na, "type": "bar", "stack": "total", "data": dd}
        print(json.dumps(x), ",")
