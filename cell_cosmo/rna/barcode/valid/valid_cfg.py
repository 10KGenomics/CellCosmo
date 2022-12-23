#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : valid_cfg.py
@Time       : 2022/06/01
@Version    : 1.0
@Desc       : None
"""
import os
from configparser import ConfigParser


class ValidCfg:
    def __init__(self, config_file):
        self._s = "chemistry"  # 约定 section
        self._o = ["barcode", "link"]  # 约定 options
        self.filepath = config_file
        self._valid()
        self.cfg = self._init_cfg()
        self.directory, _ = os.path.split(self.filepath)
        self.barcode_files = self._get(self._o[0])
        self.link_files = self._get(self._o[1])

    def _valid(self):
        if str(self.filepath).lower() in {"none", ""} or \
                not os.path.exists(self.filepath):
            # 约定 config 文件必须存在, TODO 需要交流规范后的配置文件名称,
            # current fmt: section=chemistry,options=barcode,link
            raise Exception("Can't find chemistry config file")
        if not os.path.isfile(self.filepath):
            # 这里要求: 指定到文件(非目录)
            raise Exception(f"{self.filepath} is not a file")
        # 将路径强制转换成绝对路径
        self.filepath = os.path.abspath(self.filepath)

    def _init_cfg(self):
        cfg = ConfigParser()
        cfg.read(self.filepath)
        if self._s not in cfg.sections():
            raise Exception(f"config file must contain `{self._s}` section")
        options = cfg.options(self._s)
        for o in self._o:
            if o not in options:
                raise Exception(f"the `{self._s}` section must contain "
                                f"`{'`, `'.join(self._o)}` option")
        return cfg

    def _get(self, val):
        files = self.cfg.get(self._s, val).strip()
        return files, self.directory
