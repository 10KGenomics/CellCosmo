#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : reads_processor.py
@Time       : 2024/07/05
@Version    : 1.0
@Desc       : None
"""
import time
import logging
import sqlite3
from datetime import timedelta
from cell_cosmo.tools.chemistry import LibraryInfo, get_sequence_by_pattern
from .stat_info import StatInfo
from .validators import Validators
from .reads_write import reads_write

logger = logging.getLogger(__name__)


def reads_processor(
        validators: Validators,
        library_info: LibraryInfo,
        one_batch,
):
    start = time.time()
    chunk_id, chunk_start, chunk = one_batch
    chunk_size = len(chunk)
    batch_data = {v.name: [] for v in validators.list}

    stat_info = StatInfo()
    stat_info.total_num += chunk_size

    conn = sqlite3.connect(library_info.sqlite_file)
    cursor = conn.cursor()

    for i, reads_pair in enumerate(chunk):
        r1_name, r1_seq, _, r1_qual, r2_name, r2_seq, _, r2_qual = reads_pair
        is_pass, barcode_replace, barcode_error, link_error = True, None, False, False
        reads_no = chunk_start + i + 1

        if validators.params_polyt.is_valid:
            polyt = get_sequence_by_pattern(r1_seq, library_info.pattern_T)
            if polyt.count("T") < validators.params_polyt.n_allow:
                is_pass = False
                batch_data[validators.params_polyt.name].append(reads_pair)
                stat_info.num_for_no_polyt += 1
        if validators.params_barcode.is_valid:
            seq = "".join([r1_seq[s:e] for s, e in library_info.pattern_C])
            r = cursor.execute(f"select * from barcode_set where seq = '{seq}';")
            r = r.fetchone()
            if r is None:
                r = cursor.execute(f"select * from barcode_dict where mismatch_seq = '{seq}';")
                r = r.fetchone()
                if r is None:
                    # not fond in library, this is a error reads
                    barcode_error = True
                else:
                    _, barcode_replace = r
            if barcode_error:
                is_pass = False
                batch_data[validators.params_barcode.name].append(reads_pair)
                stat_info.num_for_no_barcode += 1
        if validators.params_link.is_valid:
            seq = "".join([r1_seq[s:e] for s, e in library_info.pattern_L])
            r = cursor.execute(f"select seq from link_set where seq = '{seq}';")
            r = r.fetchone()
            if r is None:
                r = cursor.execute(f"select * from link_dict where mismatch_seq = '{seq}';").fetchone()
                if r is None:
                    link_error = True
                else:
                    _, link_replace = r
            if link_error:
                is_pass = False
                batch_data[validators.params_link.name].append(reads_pair)
                stat_info.num_for_no_link += 1
        if validators.params_qual.is_valid:
            qual_str = get_sequence_by_pattern(r1_seq, library_info.pattern_C + library_info.pattern_U)
            lows = [q for q in qual_str if ord(q) < validators.params_qual.low_qual]
            if len(lows) > validators.params_qual.n_allow:
                is_pass = False
                batch_data[validators.params_qual.name].append(reads_pair)
                stat_info.num_for_low_qual += 1

        if is_pass:
            # todo 处理reads
            umi = get_sequence_by_pattern(r1_seq, library_info.pattern_U)

            if barcode_replace is not None:
                # 对barcode进行校正
                barcode = barcode_replace
            else:
                barcode = get_sequence_by_pattern(r1_seq, library_info.pattern_C)

            qual_b = get_sequence_by_pattern(r1_qual, library_info.pattern_C)
            qual_umi = get_sequence_by_pattern(r1_qual, library_info.pattern_U)
            stat_info.qual_counter_barcode.update(qual_b)
            stat_info.qual_counter_umi.update(qual_umi)
            stat_info.qual_counter_read.update(r2_qual)

            name = f"@{barcode}_{umi}_{reads_no}"
            reads_pair = (name, r1_seq, '+', r1_qual, name, r2_seq, '+', r2_qual,)
            batch_data[validators.params_cleaned.name].append(reads_pair)
            stat_info.clean_num += 1

    cursor.close()
    conn.close()

    reads_write(
        chunk_id=chunk_id,
        batch_data=batch_data,
        validators=validators
    )
    end = time.time()
    logger.info(
        f"chunk_id={chunk_id}(s={chunk_start},n={chunk_size}) "
        f"completed. TimeUsage: {timedelta(seconds=end - start)}."
    )
    return stat_info
