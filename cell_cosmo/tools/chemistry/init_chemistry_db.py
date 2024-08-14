#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : init_chemistry_db.py
@Time       : 2024/07/10
@Version    : 1.0
@Desc       : None
"""
import os
import sqlite3
import logging
from itertools import combinations, product

if os.environ.get("CellCosmoRunMode", None) == "TEST":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(filename)s %(levelname)s %(message)s',
        datefmt='%a %d %b %Y %H:%M:%S',
    )

logger = logging.getLogger(__name__)


def yield_all_mismatch_seq(seq, n_mismatch=1, bases="A C G T N"):
    bases = [b for b in list(set(bases.upper())) if b.strip()]
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for pos in combinations(range(seq_len), n_mismatch):
        seq_pos = [[base] for base in seq]
        for loc in pos:
            seq_pos[loc] = bases
        for poss in product(*seq_pos):
            yield ''.join(poss)


def read_seq(files, pattern_list):
    seqs = []
    for filepath, (s, e) in zip(files, pattern_list):
        res, n = [], e - s
        with open(filepath, "r", encoding="utf8") as fh:
            for line in fh:
                line = line.strip()
                if line == "":
                    continue
                if len(line) != n:
                    raise Exception(
                        f"{filepath} "
                        f"should have {n} char one line !")
                res.append(line)
            assert len(res) > 0, f"the file is empty! {filepath} "
            seqs.append(list(set(res)))  # 文件内部序列去重一下
    return seqs


def build_library(conn, cursor, name, files, n_allow, pattern_list):
    seqs = read_seq(files, pattern_list)
    max_len = max([e - s for (s, e) in pattern_list])

    cursor.execute(f'create table {name}_set (seq varchar({max_len}) PRIMARY KEY);')
    cursor.execute(f'create table {name}_dict (mismatch_seq varchar({max_len}) PRIMARY KEY,seq varchar({max_len}) );')

    conn.commit()
    library_set, library_mismatch = set(), {}
    logger.info(f"init pattern {pattern_list} start ... ")
    # 将多个库按顺序进行组合
    for seq in product(*seqs):
        seq = ''.join(seq)
        library_set.add(str(seq))
        for cb in yield_all_mismatch_seq(seq, n_allow):
            library_mismatch[cb] = seq

    logger.info("init complete!!")

    sql1 = f"INSERT INTO {name}_set (seq) VALUES (?);"
    cursor.executemany(sql1, [[d] for d in list(library_set)])
    conn.commit()
    logger.info(f"insert {name}_set complete!")

    cache = []
    i, batch_size = 0, 500000
    logger.info(f"(total={len(library_mismatch)})insert {name}_dict start ...")
    for k, v in library_mismatch.items():
        cache.append([k, v])
        if len(cache) == batch_size:
            i += 1
            sql2 = f"INSERT INTO {name}_dict (mismatch_seq,seq) VALUES (?,?);"
            cursor.executemany(sql2, cache)
            conn.commit()
            logger.info(f"load {batch_size * i} compete!")
            cache = []
    if cache:
        sql2 = f"INSERT INTO {name}_dict (mismatch_seq,seq) VALUES (?,?);"
        cursor.executemany(sql2, cache)
        conn.commit()
        logger.info(f"load {batch_size * i + len(cache)} compete!")

    logger.info(f"insert {name}_dict complete!")


def init_chemistry_db(db, barcode_args: tuple, link_args: tuple):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    build_library(conn, cursor, "barcode", *barcode_args)
    build_library(conn, cursor, "link", *link_args)
    cursor.close()
    conn.close()
    # todo 通过barcodes 初始化 set,dict


if __name__ == '__main__':
    import time

    _db_file = "./t1.sqlite"
    _library = "KitVersion1"
    # _library = "DemoV1"
    _root1 = "/home/lx/Documents/gitlab/wancheng/cellcosmo/cell_cosmo"
    _root2 = "/ceph_disk3/lx/miniconda3/envs/cellcosmo_root/lib/python3.10/site-packages/cell_cosmo"
    _base = f'{_root2}/cc_data/chemistry/{_library}/C8L6C8L6C8U8T30'
    bs = [os.path.join(_base, f"Barcode{i}.list") for i in range(1, 4)]
    ls = [os.path.join(_base, f"Link{i}.list") for i in range(1, 3)]

    start = time.time()
    init_chemistry_db(
        _db_file,
        (
            bs, 1, [(0, 8), (14, 22), (28, 36)]
        ), (
            ls, 2, [(8, 14), (22, 28)]
        )
    )
    end = time.time()
    print(end - start)