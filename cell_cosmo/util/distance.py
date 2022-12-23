#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : distance.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""


def hamming_correct(str1, str2):
    threshold = len(str1) / 10 + 1
    if hamming_distance(str1, str2) < threshold:
        return True
    return False


def hamming_distance(str1, str2):
    distance = 0
    str1_len = len(str1)
    str2_len = len(str2)
    if str1_len != str2_len:
        raise Exception(f"string1({str1_len}) and string2({str2_len}) do not have same length")
    for i in range(str1_len):
        if str1[i] != str2[i]:
            distance += 1
    return distance
