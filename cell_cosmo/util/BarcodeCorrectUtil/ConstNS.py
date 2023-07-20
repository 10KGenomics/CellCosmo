#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : const_ns.py
@Time       : 2023/07/19
@Version    : 1.0
@Desc       : None
"""

# 字符串常量
barcode = "barcode"
gene_id = "gene_id"
umi = "umi"
count = "count"

intended_barcode = "intended_barcode"
intended_size = "intended_size"
intended_base = "intended_base"

neighbor_barcode = "neighbor_barcode"
neighbor_size = "neighbor_size"
neighbor_base = "neighbor_base"

position = "position"
filter_reason = "filter_reason"

percent = "percent"
neighbors = "neighbors"
n_neighbors = "n_neighbors"
neighbor_size_total = "neighbor_size_total"


def get_0_cols():
    return [barcode, gene_id, umi]


def get_1_cols():
    # COUNT_DETAIL_COLS = [barcode, gene_id, umi, count]
    return [barcode, gene_id, umi, count]
