#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : base_runner.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
from abc import ABCMeta, abstractmethod


class BaseRunner(metaclass=ABCMeta):
    def __init__(self):
        pass

    def get(self, var):
        assert hasattr(self, var), "正常情况一定存在,由调用方确认,这里只当做编程过程中的一个提醒"
        return getattr(self, var)

    @abstractmethod
    def run(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass
