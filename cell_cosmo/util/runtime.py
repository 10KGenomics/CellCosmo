#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : runtime.py
@Time       : 2022/06/01
@Version    : 1.0
@Desc       : None
"""
import time
import logging
from functools import wraps
from datetime import timedelta


def runtime(logger_name=None):
    if logger_name is None:
        logger_name = __name__
    # fmt_str = '%(asctime)s,%(msecs)03d %(levelname)s %(name)s >> %(message)s'
    # fmt = logging.Formatter(fmt=fmt_str, datefmt='%Y-%m-%dT%H:%M:%S')
    # 2022-06-01T15:20:03,846 [INFO] E:\workspace\genostack\cell_cosmo\cell_cosmo\util\runtime.py.__main__ start...
    logger = logging.getLogger(logger_name)
    # logger.setLevel(logging.INFO)
    # default_handler = logging.StreamHandler()
    # default_handler.setFormatter(fmt)
    # logger.addHandler(default_handler)

    def wrapper(func):
        @wraps(func)
        def inner(*args, **kwargs):
            logger.info('start...')
            start = time.time()
            res = func(*args, **kwargs)
            end = time.time()
            used = timedelta(seconds=end - start)
            logger.info('done. time used: %s', used)
            return res

        return inner

    return wrapper
