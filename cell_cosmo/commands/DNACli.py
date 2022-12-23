#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : DNACli.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.commands.BaseCli import BaseCli


class DNACli(BaseCli):
    __package_name__ = "dna"

    def __init__(self, *args, **kwargs):
        super(DNACli, self).__init__(*args, **kwargs)
