#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : CMDMeta.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
import os
import logging
from typing import List
from configparser import ConfigParser

logger = logging.getLogger(__name__)


class CMDMeta:
    def __init__(self, cmd_prefix, default: str = "",
                 required: bool = False, choice=None,
                 is_flag: bool = False,
                 type_=None):
        self.cmd_prefix = cmd_prefix
        self.option = str(cmd_prefix).lstrip("-").replace("-", "_")
        self.default = default
        self.required = required
        self.choice = choice
        self.is_flag = is_flag

    def __str__(self):
        fs = [f"{f}=`{getattr(self, f)}`" for f in dir(self) if not f.startswith("_")]
        return ", ".join(fs)

    def __repr__(self):
        return self.__str__()


class PipeLineConfigParser(ConfigParser):
    G = "global"

    def __init__(self, config_path):
        super(PipeLineConfigParser, self).__init__(inline_comment_prefixes=('#', ';'))
        self.read(config_path)
        # 获取全局参数
        self.fq1 = self.get_from_global("fq1", required=True)
        self.fq2 = self.get_from_global("fq2", required=True)
        self.genomeDir = self.get_from_global("genomeDir", default=os.getcwd())
        self.outdir = self.get_from_global("outdir", default=os.getcwd())
        self.sample = self.get_from_global("sample", default=os.path.basename(self.fq1).split("_")[0])
        self.thread = self.get_from_global("thread", default='4')
        self.gzip = self.get_from_global("gzip", default='0', isFlag=True)

        logger.info("\tRUN scRNA Pipeline,read global param is:")
        logger.info(f"\t\t fq1 = {self.fq1}")
        logger.info(f"\t\t fq2 = {self.fq2}")
        logger.info(f"\t\t sample = {self.sample}")
        logger.info(f"\t\t genomeDir = {self.genomeDir}")
        logger.info(f"\t\t outdir = {self.outdir}")
        logger.info(f"\t\t thread = {self.thread}")
        logger.info(f"\t\t gzip = {self.gzip}")

        if not os.path.exists(self.outdir):
            try:
                os.mkdir(self.outdir)
            except Exception as e:
                raise Exception(f"outdir path not exists({self.outdir}),"
                                f"and failed to create the folder!{str(e)}")

    def _clean_val(self, val: str):
        # 去除注释信息
        val = val.strip()
        if val and val[0] in self._inline_comment_prefixes:
            for c in self._inline_comment_prefixes:
                val = val.split(c, 1)[0]
                val = val.strip()
        return val

    def get_from_cmd_meta(self, sec, m: CMDMeta):
        """解析指定 section 中 CMDMeta 配置对应的命令行"""
        # CMDMeta中配置但是配置文件中不存在,该情况下先忽略,只提示
        # (thread和gzip两个参数例外,这两个参数可以在全局配置中设置)
        if sec not in self.sections():
            raise Exception(f"please do not modify section "
                            f"name,the `{sec}` unrecognized")
        if m.option not in self.options(sec):
            if m.option == "thread" and m.default != self.thread:
                return [m.cmd_prefix, self.thread]
            if m.option == "gzip" and self.gzip:
                # gzip is flag param
                return [m.cmd_prefix]
            if m.option in {"thread", "gzip"}:
                logger.info(f"{sec}.{m.option} use {self.G}.{m.option}:{getattr(self, m.option)}")
            else:
                logger.info(f"ignore: {sec}.{m.option}")
            return []

        val = self.get(sec, m.option)
        val = self._clean_val(val)
        if val == "" and m.required:
            raise Exception(f"{sec}.{m.option} is required,but get empty!")
        if m.is_flag:
            if val == "":
                # '--shout/--no-shout' case, 暂不考虑
                # '/debug;/no-debug' case, 暂不考虑
                # 没有设置参数,不进行拼接
                # 该判断 避免没有设置参数得到空字符串导致转换失败
                return []
            val = self._convert_to_boolean(val)
            if val and m.default == "False":
                # --gzip case, and default=False in option
                return [m.cmd_prefix]
            if not val and m.default == "True":
                # -x case, and default=True in option
                return [m.cmd_prefix]
        else:
            if val and val != m.default:
                # 该参数设置了,并且不是默认值,则拼接到命令行
                return [m.cmd_prefix, val]
        return []

    def get_from_global(self, opt, required=False, default=None, isFlag=False):
        """解析global中配置的全局参数"""
        if self.G not in self.sections():
            raise Exception(f"please do not modify section "
                            f"name `{self.G}`")
        if opt not in self.options(self.G):
            raise Exception(f"please do not modify `{self.G}` option "
                            f"name,the `{self.G}.{opt}` unrecognized")
        val = self.get(self.G, opt)
        val = self._clean_val(val)
        if val == "" and required:
            raise Exception(f"{self.G}.{opt} is required,but get empty!")
        if val == "" and default is not None:
            val = default
        if isFlag:
            # 避免没有设置参数得到空字符串导致转换失败, 默认为 false
            val = '0' if val == "" else val
            val = self._convert_to_boolean(val)
        return val

    def optionxform(self, optionStr):
        # 避免 genomeDir 配置读取错误,重载该方法避免忽略大小写
        return optionStr

    def _get_file_path(self, subName):
        return os.path.join(self.result_dir, f"{self.sample}{subName}")


class CMDBase:
    def __init__(self, section: str, meta: List[CMDMeta], cmd: List[str]):
        self.section = section
        self.meta = meta
        self.cmd = cmd
        self._gzip = None

    def getcmd(self, cfg: PipeLineConfigParser) -> str:
        """从配置文件中解析所需要的参数"""
        for m in self.meta:
            try:
                res = cfg.get_from_cmd_meta(self.section, m)
                if m.option == "gzip":
                    self._gzip = bool(res)
                self.cmd.extend(res)
            except Exception as e:
                print(str(e))
        return " ".join(self.cmd)

    @property
    def gzip(self):
        assert self._gzip is not None, \
            f"{self.__class__.__name__} gzip should " \
            f"not be none,please check it"
        return self._gzip
