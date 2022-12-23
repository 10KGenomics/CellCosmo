#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : dna1.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""

import click


@click.command()
@click.option('--debug/--no-debug', default=False)
def dna1(debug):
    """
    xxxxxxxxxxxxx
    """
    print(f"dna1 {debug}")
