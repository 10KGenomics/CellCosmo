#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : __main__.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import click
import logging
from cell_cosmo.commands import MyCli
from cell_cosmo import PACKAGE_NAME, __VERSION__


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(__VERSION__)
    ctx.exit()


def init_logger():
    fmt_str = '%(asctime)s,%(msecs)03d %(levelname)s %(name)s >> %(message)s'
    # fmt_str = '%(asctime)s,%(msecs)03d [%(process)d]%(filename)s %(levelname)s  >> %(message)s'
    fmt = logging.Formatter(fmt=fmt_str, datefmt='%Y-%m-%dT%H:%M:%S')
    # 2022-06-01T15:20:03,846 [INFO] E:\workspace\genostack\cell_cosmo\cell_cosmo\util\runtime.py.__main__ start...
    logger = logging.getLogger(PACKAGE_NAME)
    logger.setLevel(logging.INFO)
    default_handler = logging.StreamHandler()
    default_handler.setFormatter(fmt)
    logger.addHandler(default_handler)
    # logger.info("logger init complete! welcome to use cell_cosmo.")


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'], show_default=True)


@click.group(cls=MyCli, context_settings=CONTEXT_SETTINGS)
@click.option('-v', '--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True, help="show program's version and exit.")
def cli():
    """
    万乘软件
    """
    init_logger()


if __name__ == '__main__':
    cli()
