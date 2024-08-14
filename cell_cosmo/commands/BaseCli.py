#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : BaseCli.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import os
import click
import typing as t
import itertools
import importlib
from pathlib import Path
from click import Context, Command

COMMANDS_PACKAGE = Path(__file__).parent


class BaseCli(click.MultiCommand):
    __package_name__ = None

    def __init__(self, *args, **kwargs):
        super(BaseCli, self).__init__(*args, **kwargs)
        self._builtin_cmds = {}
        self._init_cmds()

    def list_commands(self, ctx: Context) -> t.List[str]:
        return list(itertools.chain(self._builtin_cmds))

    def get_command(self, ctx: Context, cmd_name: str) -> t.Optional[Command]:
        return self._builtin_cmds[cmd_name]

    def _init_cmds(self):
        if self.__package_name__ is None:
            raise Exception("the subclass must init `__package_name__`!")
        commands_path = COMMANDS_PACKAGE / self.__package_name__
        for filename in os.listdir(commands_path):
            if filename == "__init__.py" or not filename.endswith(".py"):
                continue
            filename = filename[:-len(".py")]
            try:
                module = importlib.import_module(
                    f'cell_cosmo.commands.{self.__package_name__}.{filename}')
                if not hasattr(module, filename):
                    raise Exception('')
                func = getattr(module, filename)
                self._builtin_cmds[filename] = func
            except Exception as e:
                print(str(e))
        # print(self.__package_name__, "====", self._builtin_cmds)
