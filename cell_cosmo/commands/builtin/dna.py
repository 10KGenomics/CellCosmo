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
from cell_cosmo.commands import DNACli


@click.group(cls=DNACli)
def dna():
    """
    Single-cell dna.
    """
    pass
