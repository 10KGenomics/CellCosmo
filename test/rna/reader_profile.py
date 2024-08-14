#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : reader_profile.py
@Time       : 2024/07/01
@Version    : 1.0
@Desc       : None
"""
import pstats
import cProfile
from cell_cosmo.rna.barcode_1.reader1 import reader

# snakeviz load_data.out
DO_PROF = True


def do_cprofile(filename):
    """
    Decorator for function profiling.
    """

    def wrapper(func):
        def profiled_func(*args, **kwargs):
            if DO_PROF:
                profile = cProfile.Profile()
                profile.enable()
                result = func(*args, **kwargs)
                profile.disable()
                # Sort stat by internal time.
                sortby = "tottime"
                ps = pstats.Stats(profile).sort_stats(sortby)
                ps.dump_stats(filename)
            else:
                result = func(*args, **kwargs)
            return result

        return profiled_func

    return wrapper


@do_cprofile('./load_data.out')
def test():
    file1 = "/mnt/data/万乘/0120W_0324-mb_L04_R1.fq.gz"
    file2 = "/mnt/data/万乘/0120W_0324-mb_L04_R2.fq.gz"

    reader(file1, file2)


test()
