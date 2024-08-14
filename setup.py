#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : setup.py.py
@Time       : 2022/05/20
@Version    : 1.0
@Desc       : None
"""
from setuptools import find_packages, setup
from cell_cosmo import __VERSION__, PACKAGE_NAME

with open("README.md", "r", encoding="UTF8") as fh:
    LONG_DESCRIPTION = fh.read()
with open('requirements.txt', "r", encoding="UTF8") as fh:
    install_requires = fh.read()
setup(
    name=PACKAGE_NAME,
    version=__VERSION__,
    author="ice-melt",
    author_email="ice-melt@outlook.com",
    maintainer="ice-melt",
    maintainer_email="ice-melt@outlook.com",
    description='万乘单细胞分析流程软件',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url='http://console.genostack.com/',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Free for non-commercial use",
        "Operating System :: OS Independent",
    ],
    # 需要安装的依赖
    # install_requires=[],
    # package_dir={"": "src"},
    packages=find_packages(exclude=["update_package.py"]),
    # package_data={"": ["*.*"]},
    include_package_data=True,
    entry_points={'console_scripts': ['CellCosmo=cell_cosmo.__main__:cli', ]},
    install_requires=install_requires
)
