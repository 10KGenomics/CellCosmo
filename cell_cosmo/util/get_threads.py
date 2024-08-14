#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : get_threads.py
@Time       : 2024/07/01
@Version    : 1.0
@Desc       : None
"""
import logging
import multiprocessing

logger = logging.getLogger(__name__)


def get_threads(thread: int, min_limit=None):
    """

    Args:
        thread:
        min_limit: 最小要求线程数，默认0

    Returns:

    """
    c = multiprocessing.cpu_count() * 2
    if min_limit is not None and thread < min_limit:
        if min_limit > c:
            logger.error(f"程序运行最少需要{min_limit}线程,当前机器不满足运行条件")
            raise Exception(f"程序运行最少需要{min_limit}线程,当前机器不满足运行条件")

        thread = min_limit

    if thread > c:
        logger.warning(f"设置的线程数({thread})大于机器最大线程数({c}),忽略线程配置，使用 thread={c}")
        thread = c
    return thread
