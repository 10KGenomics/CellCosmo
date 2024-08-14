#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : get_config_path_by_chemistry_name.py
@Time       : 2024/07/10
@Version    : 1.0
@Desc       : None
"""
import os
from pathlib import Path
from cell_cosmo import ROOT_PATH

DEFAULT_CONFIG_FILE = "chemistry.ini"


def get_config_path_by_name(name) -> Path:
    """通过文库名称获取文库信息(pattern_str, config_file)"""
    # chemistry 位置位于项目根目录的 cc_data 下
    chemistry_root = ROOT_PATH / "cc_data" / "chemistry"
    names = {n for n in os.listdir(chemistry_root) if (Path(chemistry_root) / n).is_dir()}
    if name not in names:
        name_str = "`, `".join(names)
        raise Exception(f"不支持该名称的文库:`{name}`,当前仅支持:`{name_str}`"
                        f",或者联系管理员增加该文库配置")
    libname_path = chemistry_root / name
    libname_dirs = []
    for p in os.listdir(libname_path):
        pattern_path = os.path.join(libname_path, p)
        config_path = libname_path / p / DEFAULT_CONFIG_FILE
        if os.path.isdir(pattern_path) and config_path.exists():
            libname_dirs.append(config_path)
    if len(libname_dirs) == 0:
        raise Exception(f"没有识别名为`{name}`的文库，请检查文库是否满足要求")
    if len(libname_dirs) != 1:
        raise Exception(f"名为`{name}`的文库中包含多个文库配置目录,当前仅支持配置一个目录")
    return libname_dirs[0]  # pattern_str, config_file


if __name__ == '__main__':
    res = get_config_path_by_name("KitVersion1")
    print(res)
