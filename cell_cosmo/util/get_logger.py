#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : get_logger.py
@Time       : 2022/06/02
@Version    : 1.0
@Desc       : None
"""
import logging


def get_logger(name):
    fmt_str = '%(asctime)s,%(msecs)03d %(levelname)s %(name)s >> %(message)s'
    fmt = logging.Formatter(fmt=fmt_str, datefmt='%Y-%m-%dT%H:%M:%S')
    # 2022-06-01T15:20:03,846 [INFO] E:\workspace\genostack\cell_cosmo\cell_cosmo\util\runtime.py.__main__ start...
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    default_handler = logging.StreamHandler()
    default_handler.setFormatter(fmt)
    logger.addHandler(default_handler)
    return logger
