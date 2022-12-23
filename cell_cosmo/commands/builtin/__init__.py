#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : __init__.py.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
# import os
# import click
# import typing as t
# from click import Context, Command
# from cell_cosmo import ROOT_PATH
# import itertools
# import importlib
#
#
# class BuiltinCli(click.MultiCommand):
#     def __init__(self, sub_package_name, *args, **kwargs):
#         super(BuiltinCli, self).__init__(*args, **kwargs)
#         self._sub_package_name = sub_package_name
#         self._commands = {}
#         self._init_cmds()
#
#     def list_commands(self, ctx: Context) -> t.List[str]:
#         return itertools.chain(self._commands)
#
#     def get_command(self, ctx: Context, cmd_name: str) -> t.Optional[Command]:
#         return self._commands[cmd_name]
#
#     def _init_cmds(self):
#         commands_pkg = os.path.join(ROOT_PATH, "commands")
#         commands_path = os.path.join(commands_pkg, self._sub_package_name)
#
#         for f in os.listdir(commands_path):
#             if f == "__init__.py" or not f.endswith(".py"):
#                 continue
#             f_name = f[:-len(".py")]
#             func = importlib.import_module(f'cell_cosmo.commands.{self._sub_package_name}.{f_name}')
#             try:
#                 self._commands[f_name] = func
#             except Exception as e:
#                 print(f_name)
#                 print(str(e))
