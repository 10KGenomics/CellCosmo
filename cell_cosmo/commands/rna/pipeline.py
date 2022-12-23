#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : pipeline.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
import os
import click
import shutil
from cell_cosmo import ROOT_PATH
from cell_cosmo.rna.Pipeline import Pipeline


def gen_cfg_tpl_file(ctx: click.core.Context, param: click.Option, value):
    if not value or ctx.resilient_parsing:
        return
    src = os.path.join(ROOT_PATH, "rna_pipeline.tpl.cfg")
    dst = os.path.join(os.getcwd(), "rna_pipeline.cfg")
    shutil.copyfile(src, dst)
    ctx.exit()


@click.command()
@click.option("-g", "--gen-config-tpl", is_flag=True, callback=gen_cfg_tpl_file, help="Generate config template file.")
@click.option("-c", "--config", help="Config file for the pipeline.")
def pipeline(config, **kwargs):
    Pipeline(config).run()
