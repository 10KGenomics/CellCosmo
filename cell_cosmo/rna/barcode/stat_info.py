#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : stat_info.py
@Time       : 2024/07/03
@Version    : 1.0
@Desc       : None
"""
from collections import Counter


class StatInfo:
    def __init__(self):
        self.total_num = 0
        self.clean_num = 0
        self.num_for_no_link = 0
        self.num_for_no_barcode = 0
        self.num_for_no_polyt = 0
        self.num_for_low_qual = 0
        self.qual_counter_barcode = Counter()
        self.qual_counter_umi = Counter()
        self.qual_counter_read = Counter()

    def update(self, stat_info):
        assert isinstance(stat_info, StatInfo), "def assert"
        self.clean_num += stat_info.clean_num
        self.num_for_no_link += stat_info.num_for_no_link
        self.num_for_no_barcode += stat_info.num_for_no_barcode
        self.num_for_no_polyt += stat_info.num_for_no_polyt
        self.num_for_low_qual += stat_info.num_for_low_qual
        self.total_num += stat_info.total_num
        self.qual_counter_umi += stat_info.qual_counter_umi
        self.qual_counter_read += stat_info.qual_counter_read
        self.qual_counter_barcode += stat_info.qual_counter_barcode
