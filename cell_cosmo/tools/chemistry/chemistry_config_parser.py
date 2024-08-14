#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : chemistry_config_parse.py
@Time       : 2024/07/10
@Version    : 1.0
@Desc       : None
"""
import os
import logging
from pathlib import Path
from configparser import ConfigParser

logger = logging.getLogger(__name__)


def _check_config_file(config_file: str) -> Path:
    # 将路径强制转换成绝对路径
    config_file = Path(config_file).absolute()
    if not config_file.exists():
        # 约定 config 文件必须存在
        raise Exception("Can't find chemistry config file")
    if not os.path.isfile(config_file):
        # 这里要求: 指定到文件(非目录)
        raise Exception(f"{config_file} is not a file")
    return config_file


class ChemistryConfigParser(ConfigParser):
    SECTION = "chemistry"
    LINK = "link"
    BARCODE = "barcode"
    OPTIONS = [BARCODE, LINK]

    def __init__(self, config_path):
        """
        约定配置文件格式: section=chemistry,options=barcode,link;
        eg:
        [chemistry]
        barcode=Barcode1.list,Barcode2.list,Barcode3.list
        link=Link1.list,Link2.list
        """
        super(ChemistryConfigParser, self).__init__()
        self.config_path = _check_config_file(config_path)
        self._chemistry_root = self.config_path.parent

        self.read(self.config_path)

        self._check()

    def _check(self):
        options = self.options(self.SECTION)
        for o in self.OPTIONS:
            if o not in options:
                raise Exception(
                    f"the `{self.SECTION}` section must contain "
                    f"`{'`, `'.join(self.OPTIONS)}` option")

    def get_barcode_files(self) -> list:
        return self._get_files(self.BARCODE)

    def get_link_files(self) -> list:
        return self._get_files(self.LINK)

    def _get_files(self, option):
        filepaths = []
        files = self.get(self.SECTION, option).strip()
        for filename in files.split(","):
            filename = filename.strip()
            fp = os.path.join(self._chemistry_root, filename)
            if not os.path.exists(fp):
                raise Exception(f"config error!`{option}`:"
                                f"`{filename}` not exists in {self._chemistry_root}")
            filepaths.append(fp)
        return filepaths
