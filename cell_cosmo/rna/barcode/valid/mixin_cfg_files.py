#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : mixin_cfg_files.py
@Time       : 2022/06/02
@Version    : 1.0
@Desc       : None
"""
import os


class CfgFilesMixin:
    def join_path(self, directory, filename):
        filepath = os.path.join(directory, filename.strip())
        if not os.path.exists(filepath):
            raise Exception(f"{self.__class__.__name__}: config file `{filename}` not exists in {directory}")
        return filepath

    def read_list(self, filepath, n):
        res = []
        with open(filepath, "r", encoding="utf8") as fh:
            for line in fh:
                line = line.strip()
                if line == "":
                    continue
                if len(line) != n:
                    # f'length of L in pattern ({length}) do not equal to length in {seq_file} ({len(seq)}) !'
                    raise Exception(f"{self.__class__.__name__}:{filepath} should have {n} char one line !")
                res.append(line)
            assert len(res) > 0, f"the file is empty! {filepath} "
            return res
