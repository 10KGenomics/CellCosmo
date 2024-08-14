#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : chemistry.py
@Time       : 2024/07/10
@Version    : 1.0
@Desc       : None
"""
import os.path
import re
import logging
from pathlib import Path
from typing import List, Tuple, Set
from dataclasses import dataclass, field
from .chemistry_config_parser import ChemistryConfigParser
from .get_config_path_by_name import get_config_path_by_name
from .init_chemistry_db import init_chemistry_db
from .parse_pattern import parse_pattern

logger = logging.getLogger(__name__)

C, U, L, T, N = "CULTN"


def read_seq(files, pattern_list):
    seqs = []
    for filepath, (s, e) in zip(files, pattern_list):
        res, n = set(), e - s
        with open(filepath, "r", encoding="utf8") as fh:
            for line in fh:
                line = line.strip()
                if line == "":
                    continue
                if len(line) != n:
                    raise Exception(
                        f"{filepath} "
                        f"should have {n} char one line !")
                res.add(line)
            assert len(res) > 0, f"the file is empty! {filepath} "
            seqs.append(res)
    return seqs


@dataclass
class LibraryInfo:
    library_name: str | None = None  # 文库名称
    pattern_str: str | None = None  # 文库模式字符串，eg: C8L6C8L6C8U8T30
    config_file: str | None = None  # chemistry.ini 文件的位置

    barcode_n_mismatch: int = 1
    link_n_mismatch: int = 2
    auto_init_db: bool = True
    force_re_init_db: bool = False
    pattern_C: List[Tuple[int, int]] = field(init=False)
    pattern_L: List[Tuple[int, int]] = field(init=False)
    pattern_U: List[Tuple[int, int]] = field(init=False)
    pattern_T: List[Tuple[int, int]] = field(init=False)
    pattern_N: List[Tuple[int, int]] | None = field(init=False)
    pattern_list: List[Tuple[str, int, int]] = field(init=False)
    n_C: int = field(init=False)
    n_L: int = field(init=False)
    n_U: int = field(init=False)
    n_T: int = field(init=False)
    n_N: int | None = field(init=False)

    link_files: List[str] = field(init=False)
    link_library: List[Set[str]] = field(init=False)
    barcode_files: List[str] = field(init=False)
    barcode_library: List[Set[str]] = field(init=False)
    sqlite_file: str | None = None

    def __post_init__(self):
        self._check()
        self._init_pattern_field()
        self._init_config_field()

    def _check(self):
        if self.library_name is not None:
            # 优先以文库名识别相关信息
            if self.config_file is not None \
                    or self.pattern_str is not None:
                # todo 这里的验证放到命令行中
                logger.warning("当指定`--chemistry-name`时, "
                               "`--chemistry-config` 和 `--pattern` 无效，"
                               "已忽略这些参数值!")
            if self.library_name == "default":
                self.library_name = "KitVersion1"  # 默认使用该版本的库
            config_path = get_config_path_by_name(self.library_name)
            self.config_file = str(config_path)
            self.pattern_str = config_path.parent.name
        else:
            if self.config_file is None:
                raise Exception(f"未通过`--chemistry-name`指定文库名时,"
                                f"`--chemistry-config`必须指定")

            c = Path(self.config_file)
            if not c.exists():
                raise Exception(f"`--chemistry-config`指定的文件不存在！")

            if self.pattern_str is None:
                # todo 校验父文件夹是否满足模式字符串的要求
                self.pattern_str = c.parent.name
            self.library_name = self.pattern_str  # or = c.parent.parent.name?

    def _init_config_field(self):
        cfg = ChemistryConfigParser(self.config_file)
        self.link_files = cfg.get_link_files()
        self.barcode_files = cfg.get_barcode_files()

        self.link_library = read_seq(self.link_files, self.pattern_L)
        self.barcode_library = read_seq(self.barcode_files, self.pattern_C)

        library_dir = Path(self.config_file).absolute().parent
        self.sqlite_file = str(library_dir / ".mismatch_library.sqlite")
        if not os.path.exists(self.sqlite_file) and self.auto_init_db:
            # todo 初始化 sqlite
            init_chemistry_db(
                self.sqlite_file,
                (
                    self.barcode_files,
                    self.barcode_n_mismatch,
                    self.pattern_C
                ),
                (
                    self.link_files,
                    self.link_n_mismatch,
                    self.pattern_L
                )
            )

    def _init_pattern_field(self):
        pattern_dict, pattern_list = parse_pattern(self.pattern_str)

        self.pattern_list = pattern_list
        self.pattern_C = pattern_dict.get(C)
        self.pattern_L = pattern_dict.get(L)
        self.pattern_U = pattern_dict.get(U)
        self.pattern_T = pattern_dict.get(T)
        self.pattern_N = pattern_dict.get(N, None)
        self.n_C = sum([e - s for s, e in self.pattern_C])
        self.n_L = sum([e - s for s, e in self.pattern_L])
        self.n_U = sum([e - s for s, e in self.pattern_U])
        self.n_T = sum([e - s for s, e in self.pattern_T])
        self.n_N = None if self.pattern_N is None else sum([e - s for s, e in self.pattern_N])
