#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : keynames.py
@Time       : 2024/07/06
@Version    : 1.0
@Desc       : None
"""
from enum import Enum, unique


@unique
class KeyName(Enum):
    # name for output file,value for get data by index
    no_polyt = 1
    no_barcode = 2
    no_linker = 3
    low_qual = 4
    clean_num = 5
