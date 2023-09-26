#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : base_output.py
@Time       : 2022/06/02
@Version    : 1.0
@Desc       : None
"""
import os
from xopen import xopen


class BaseOutput:
    __filename__ = None

    def __init__(self, need_valid=False, need_output=False, **kwargs):
        """
        只有需要验证且需要输出时才会进行输出处理
        :param need_valid: 是否需要验证
        :param need_output: 是否输出
        :param kwargs:
        """
        self.num = 0
        self.need_valid = need_valid
        self.need_output = need_output
        self.outdir = kwargs.pop("outdir")
        self.sample = kwargs.pop("sample")
        self.suffix = ".gz" if kwargs.pop("gzip", False) else ""
        self.prefix = os.path.join(self.outdir, self.sample)
        self.fh1, self.fh2 = None, None
        if self.need_valid and self.need_output:
            if self.__filename__ is None:
                raise Exception(f"please set __filename__ in {self.__module__}.{self.__class__}")
            self.fh1 = self._get_fh(f"{self.__filename__}1")
            self.fh2 = self._get_fh(f"{self.__filename__}2")

    def _get_fh(self, name, mode="w"):
        """自动拼接文件全路径名,并返回处理文件的函数句柄"""
        return xopen(f'{self.prefix}_{name}.fq{self.suffix}', mode)

    def _close(self):
        # 只允许子类调用
        for fh in (self.fh1, self.fh2):
            if fh is not None:
                fh.close()

    def _write(self, r1, r2, not_valid):
        # 只允许子类调用
        if self.need_valid and not_valid:
            # 需要校验且无效，统计值+1,并输出
            # TODO need_output 目前还没使用
            self.num += 1
            self.fh1.write("\n".join(r1) + "\n")
            self.fh2.write("\n".join(r2) + "\n")
