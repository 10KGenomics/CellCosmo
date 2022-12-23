#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : ImageUtil.py
@Time       : 2022/08/29
@Version    : 1.0
@Desc       : None
"""
import base64


def get_base64(img_path):
    with open(img_path, 'rb') as f:
        image_data = f.read()
        base64_data = base64.b64encode(image_data)  # base64编码
        base64_data = base64_data.decode("utf-8")
        return base64_data


def get_img_src_base64_str(img_path):
    return f"data:image/png;base64,{get_base64(img_path)}"
