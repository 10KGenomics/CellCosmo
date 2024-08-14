#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : get_sequence_by_pattern.py
@Time       : 2024/07/15
@Version    : 1.0
@Desc       : None
"""


def get_sequence_by_pattern(seq, pattern, merge=True):
    res = [seq[s:e] for s, e in pattern]
    if merge:
        res = ''.join(res)
    return res
