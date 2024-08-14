#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : new_write.py
@Time       : 2024/07/08
@Version    : 1.0
@Desc       : None
"""
import time
import logging
from xopen import xopen
from threading import Thread
from .validators import Validators

logger = logging.getLogger(__name__)


def _get_current_chunk_id(temp_file):
    with open(temp_file, "r") as fh:
        try:
            x = int(fh.read().strip())
            return x
        except Exception as e:
            # 多进程对文件读的时候，其中一个进程写文件会导致读空
            str(e)
            return -1


def _update_chunk_id(temp_file, chunk_id: int):
    with open(temp_file, "w") as fh:
        fh.write(str(chunk_id))


class Write(Thread):
    def __init__(self, f1, f2, chunk_data, chunk_id):
        super(Write, self).__init__()
        self.f1 = f1
        self.f2 = f2
        self.chunk_data = chunk_data
        self.chunk_id = str(chunk_id).zfill(5)

    def run(self):
        # s = time.time()
        # batch = 1000

        with xopen(self.f1, "a") as fh1, xopen(self.f2, "a") as fh2:
            for r1_name, r1_seq, _, r1_qual, r2_name, r2_seq, _, r2_qual in self.chunk_data:
                fh1.write(f"{r1_name}\n{r1_seq}\n+\n{r1_qual}\n")
                fh2.write(f"{r2_name}\n{r2_seq}\n+\n{r2_qual}\n")

        # e = time.time()
        # logger.info(f"index={self.chunk_id} write in {e - s}s, {Path(self.f1).name}")


def reads_write(chunk_id: int, batch_data, validators: Validators):
    while chunk_id != _get_current_chunk_id(validators.temp_state_file):
        time.sleep(0.01)

    # do write
    ts = []
    # logger.info(f"chunk_id={chunk_id} 开始写入 。。。")
    for v in validators.list:
        if not v.is_valid:
            continue  # 无需校验，跳过
        this_batch_data = batch_data[v.name]
        if not this_batch_data:
            continue  # 没有满足的数据，跳过
        t = Write(v.f1, v.f2, this_batch_data, chunk_id)
        t.start()
        ts.append(t)
    _ = [t.join() for t in ts]
    # logger.info(f"chunk_id={chunk_id} 写入完成")

    # 更新 chunk_id 进入下一个chunk_id所在进程开始写数据
    _update_chunk_id(validators.temp_state_file, chunk_id + 1)
