#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : archive_data.py
@Time       : 2022/07/01
@Version    : 1.0
@Desc       : None
"""
import os
from cell_cosmo.util.ScanpyUtil.constant import TopDirFmt
from cell_cosmo.util import PathUtil
from functools import wraps


def archive_data(subdir=None, write=None):
    def wrapper(func):
        @wraps(func)
        def inner(*args, **kwargs):
            res = func(*args, **kwargs)
            if write is not None \
                    and hasattr(write, '__call__') \
                    and hasattr(args[0], "outdir") \
                    and hasattr(args[0], "sample") \
                    and hasattr(args[0], "adata"):
                outdir = getattr(args[0], "outdir")
                sample = getattr(args[0], "sample")
                adata = getattr(args[0], "adata")
                outdir = os.path.join(outdir, TopDirFmt % sample)
                PathUtil.create_dir_if_not_exists(outdir)
                with PathUtil.chdir(outdir):
                    d = func.__name__ if subdir is None else subdir
                    if not os.path.exists(d):
                        os.mkdir(d)
                    with PathUtil.chdir(d):
                        write(adata)
            return res

        return inner

    return wrapper
