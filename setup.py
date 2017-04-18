#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name="MetaCov",
    use_scm_version=True,
    packages=find_packages(),
    setup_requires=[
        'setuptools_scm',
    ],
    install_requires=[
        'Click',
        'pysam'
    ],
    entry_points='''
        [console_scripts]
        metacov=metacov.cli:main
    '''
)
