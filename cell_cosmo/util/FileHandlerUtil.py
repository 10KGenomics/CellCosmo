#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : FileHandlerUtil.py
@Time       : 2022/06/28
@Version    : 1.0
@Desc       : None
"""


class Writer(object):
    def __init__(self, file_path: str, mode="w", encoding="UTF8", limit: int = 10000):
        self.file_path = file_path
        self._limit = limit
        self.encoding = encoding
        self.mode = mode
        self._cache = []
        self._fh = self._init_file_handler()

    def _init_file_handler(self):
        # TODO 根据输出文件名的后缀决定是否为gz格式
        return open(self.file_path, self.mode, encoding=self.encoding)

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()
        return type_ is None

    def _write(self):
        self._fh.write("\n".join(self._cache) + "\n")
        self._cache = []

    def write(self, data: str or list):
        """write single or multi data"""
        datatype = type(data)
        assert datatype in {str, list}, f"not support data type {datatype}!!!"
        if datatype == list:
            self._cache.extend(data)
        else:
            self._cache.append(data)
        if len(self._cache) > self._limit:
            self._write()

    def close(self):
        if self._cache:
            self._write()
        self._fh.close()
