#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : PathUtil.py
@Time       : 2022/06/20
@Version    : 1.0
@Desc       : None
"""
import os
import shutil
from contextlib import contextmanager


@contextmanager
def chdir(target):
    origin_dir = os.getcwd()
    os.chdir(target)
    try:
        yield
    finally:
        os.chdir(origin_dir)


def create_dir_if_not_exists(target):
    if os.path.exists(target):
        return
    try:
        os.mkdir(target)
    except Exception as e1:
        try:
            os.makedirs(target)
        except Exception as e2:
            print(f"创建目录失败:{target}")
            print(str(e1))
            print(str(e2))


def clean_dir(target):
    if os.path.exists(target):
        try:
            shutil.rmtree(target)
            os.mkdir(target)
        except Exception as e:
            if "[WinError 32]" in str(e):
                # 父级文件夹被占用,选择清理内部所有文件
                for fp in os.listdir(target):
                    shutil.rmtree(os.path.join(target, fp))
            else:
                print(f"rm dir fail :{target},{str(e)}")
                import sys
                sys.exit(1)
    else:
        # 目录不存在,则创建该目录
        create_dir_if_not_exists(target)
