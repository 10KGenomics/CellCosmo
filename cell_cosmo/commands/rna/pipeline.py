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
from cell_cosmo.rna.PipelineStarsolo import PipelineStarsolo


def gen_cfg_tpl_file(pipe_type):
    if pipe_type == "star":
        tpl = "rna_pipeline.tpl.cfg"
    else:
        tpl = "rna_starsolo_pipeline.tpl.cfg"
    src = os.path.join(ROOT_PATH, tpl)
    dst = os.path.join(os.getcwd(), tpl.replace(".tpl.", "."))
    shutil.copyfile(src, dst)


@click.command()
@click.option("-g", "--gen-config-tpl", is_flag=True, help="Generate config template file.")
@click.option("-t", "--pipe-type", default="star", help="Which type template to generate.")
@click.option("-c", "--config", default="", help="Config file for the pipeline.")
def pipeline(gen_config_tpl, pipe_type, config):
    if gen_config_tpl and config:
        raise Exception(f"不能同时指定-g和-c参数")
    assert pipe_type in {'starsolo', 'star'}, f"-t只能指定为 `starsolo` 或者 `star`"

    if gen_config_tpl:
        gen_cfg_tpl_file(pipe_type)
    else:
        if pipe_type == "star":
            Pipeline(config).run()
        else:
            PipelineStarsolo(config).run()
