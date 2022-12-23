#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : parse_pattern.py
@Time       : 2022/05/26
@Version    : 1.0
@Desc       : None
"""
import re
from typing import List
from cell_cosmo.util.runtime import runtime

C, U, L, T, N = list("CULTN")
ALLOW_CODES = {C, U, L, T, N}


class PatternParser:
    def __init__(self, pattern_str):
        # C5U3L15U3C6U3C5L6T30
        self.pattern_str = pattern_str
        self.allow_codes = ALLOW_CODES
        self.pattern_dict = {'pattern_list': []}
        self._parse_pattern()

    @runtime(__name__)
    def _parse_pattern(self):
        compile_pattern = re.compile(f"([{''.join(self.allow_codes)}]\\d+)")
        findall_pattern = re.findall(compile_pattern, self.pattern_str)
        if not findall_pattern:
            raise Exception(f"无效的pattern:{self.pattern_str}")
        if "".join(findall_pattern) != self.pattern_str:
            start, error_code = 0, ""
            for current_pattern in findall_pattern:
                if self.pattern_str[start:start + len(current_pattern)] != current_pattern:
                    error_code = self.pattern_str[start]
                    break
                start += len(current_pattern)
            raise Exception(f"Find invalid code `{error_code}` "
                            f"in pattern:`{self.pattern_str}`,please "
                            f"note that the pattern only support this "
                            f"code: [`{'`,`'.join(self.allow_codes)}`]")
        start = 0
        for current_pattern in findall_pattern:
            code, length = current_pattern[0], int(current_pattern[1:])
            end = start + length
            if code not in self.pattern_dict:
                self.pattern_dict[code] = []
            self.pattern_dict[code].append((start, end))
            self.pattern_dict['pattern_list'].append((code, start, end))
            start = end

    def get_pattern(self, code: str or List[str]):
        if type(code) == str:
            if code not in self.allow_codes:
                raise Exception(f"unsupported code!`{code}`")
            return self.pattern_dict[code]
        elif type(code) == list:
            segments = []
            for c in code:
                segments.extend(self.get_pattern(c))
            return segments
        else:
            raise Exception(f"Unsupported type!{type(code)}")

    def get_pattern_len(self, code: str or List[str]):
        return sum([e - s for s, e in self.get_pattern(code)])

    def get_seq_list(self, seq, code: str or List[str]) -> List:
        return [seq[item[0]: item[1]] for item in self.get_pattern(code)]

    def get_seq_str(self, seq, code: str or list[str]) -> str:
        return ''.join(self.get_seq_list(seq, code))


if __name__ == '__main__':
    test_str = "C5U3L15U3C6U3C5L6T30"
    p = PatternParser(test_str)
    # get_pattern
    assert p.get_pattern(C) == [(0, 5), (26, 32), (35, 40)]
    assert p.get_pattern(U) == [(5, 8), (23, 26), (32, 35)]
    assert p.get_pattern(L) == [(8, 23), (40, 46)]
    assert p.get_pattern(T) == [(46, 76)]
    assert p.get_pattern([C, U]) == [(0, 5), (26, 32), (35, 40), (5, 8), (23, 26), (32, 35)]
    # get_pattern_len
    assert p.get_pattern_len(C) == 16
    assert p.get_pattern_len(U) == 9
    assert p.get_pattern_len(L) == 21
    assert p.get_pattern_len(T) == 30
    assert p.get_pattern_len([C, U]) == 25
    test_seq_str = "GAAGCAGGGTTACAACCGTCACGTAAGCAGTATTATGTGTGCTCACTTTTTTT" \
                   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTG" \
                   "TTTGTTTTTGTTGTTTTGTTTGTTTGTTTTTTTTTGTGTTTTG "
    # get_seq_list
    assert p.get_seq_list(test_seq_str, C) == ['GAAGC', 'GCAGTA', 'TGTGT']
    assert p.get_seq_list(test_seq_str, U) == ['AGG', 'TAA', 'TTA']
    assert p.get_seq_list(test_seq_str, L) == ['GTTACAACCGTCACG', 'GCTCAC']
    assert p.get_seq_list(test_seq_str, T) == ['TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT']
    assert p.get_seq_list(test_seq_str, [C, U]) == [
        'GAAGC', 'GCAGTA', 'TGTGT', 'AGG', 'TAA', 'TTA']

    # get_seq_str
    assert p.get_seq_str(test_seq_str, C) == 'GAAGCGCAGTATGTGT'
    assert p.get_seq_str(test_seq_str, U) == 'AGGTAATTA'
    assert p.get_seq_str(test_seq_str, L) == 'GTTACAACCGTCACGGCTCAC'
    assert p.get_seq_str(test_seq_str, T) == 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
    assert p.get_seq_str(test_seq_str, [C, U]) == 'GAAGCGCAGTATGTGTAGGTAATTA'
