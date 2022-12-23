#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : rna.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
from cell_cosmo.commands import RNACli


@click.group(cls=RNACli)
def rna():
    """
    Single-cell rna.
    """
    pass
