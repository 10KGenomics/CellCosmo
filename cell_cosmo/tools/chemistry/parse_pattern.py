#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : parse_pattern.py
@Time       : 2024/07/11
@Version    : 1.0
@Desc       : None
"""
import re
from collections import defaultdict


def check_is_equal(findall_pattern, pattern):
    # [('C','8'),...] => [('C8'),...]
    findall_pattern = list(map(lambda x: x[0] + x[1], findall_pattern))
    if "".join(findall_pattern) != pattern:
        start, error_code = 0, ""
        for current_pattern in findall_pattern:
            if pattern[start:start + len(current_pattern)] != current_pattern:
                error_code = pattern[start]
                break
            start += len(current_pattern)
        if error_code == "" and len(pattern) > start:
            # 处理前面都一样，但是pattern后有额外字符的情况
            error_code = pattern[start]
        raise Exception(f"Find invalid code `{error_code}` "
                        f"in pattern:`{pattern}`,please "
                        f"note that the pattern only support this "
                        f"code: [C, L, U, N, T]")


def parse_pattern(pattern: str):
    p = re.compile(r'([CLUNT])(\d+)')
    findall_pattern = p.findall(pattern)
    if not findall_pattern:
        raise Exception(f"无效的pattern:{pattern}")
    check_is_equal(findall_pattern, pattern)

    start, pattern_dict, pattern_list = 0, defaultdict(list), []
    for c, n in findall_pattern:
        end = start + int(n)
        pattern_dict[c].append((start, end))
        pattern_list.append((c, start, end))
        start = end
    return pattern_dict, pattern_list
